"""Turn the 3D topology result into an organic, printable chassis (3 pieces).

density -> smooth field + fused bolt bosses -> isosurface -> exact flat
mating faces, holes, seam cuts. Writes stl/print + stl/assembly chassis_*.
"""
import json, math, os, sys
import numpy as np
import trimesh
from manifold3d import Manifold, Mesh, OpType
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from skimage.measure import marching_cubes
import params as P

d = np.load("topo3d.npz")
rho, X0, Y0, ZB, H = d["rho"], float(d["x0"]), float(d["y0"]), float(d["z0"]), float(d["h"])
nx, ny, nz = rho.shape
ZTOP = ZB + H * nz
Z0 = P.DECK_Z0
A_ = P.DECK_L / 2 - P.DECK_R

# ---------------- fine field ----------------
S = 1.25                                          # voxel size (mm)
cx = X0 + H * (np.arange(nx) + 0.5)
cy = Y0 + H * (np.arange(ny) + 0.5)
cz = ZB + H * (np.arange(nz) + 0.5)
fx = np.arange(X0 - 6, X0 + H * nx + 6, S)
fy = np.arange(Y0 - 6, Y0 + H * ny + 6, S)
fz = np.arange(ZB - 6, ZTOP + 6, S)
FX, FY, FZ = np.meshgrid(fx, fy, fz, indexing="ij")
interp = RegularGridInterpolator((cx, cy, cz), rho, bounds_error=False, fill_value=0.0)
f = interp(np.stack([FX.ravel(), FY.ravel(), FZ.ravel()], 1)).reshape(FX.shape)
f = gaussian_filter(f, 2.4)                        # bone-like smooth members
f = f / max(f.max(), 1e-9)

def soft(inside_dist, w=1.2):
    """0..1 indicator from signed distance (negative = inside)."""
    return 1 / (1 + np.exp(inside_dist / w))

def cyl_z(x, y, r, z0, z1):
    dr = np.hypot(FX - x, FY - y) - r
    dz = np.maximum(z0 - FZ, FZ - z1)
    return soft(np.maximum(dr, dz))

def box(x0, x1, y0, y1, z0, z1):
    dd = np.maximum.reduce([x0 - FX, FX - x1, y0 - FY, FY - y1, z0 - FZ, FZ - z1])
    return soft(dd)

bosses = []
pads = P.pod_bolts(1) + P.pod_bolts(-1)
for b in pads:
    bosses.append(cyl_z(*b, 8.0, Z0, Z0 + 11))
for xc in (P.CASTER_X, -P.CASTER_X):
    bosses.append(cyl_z(xc, 0, P.CASTER_BOLT_R + 7, Z0, Z0 + 7))
bosses.append(cyl_z(P.HUB_X, 0, P.HUB_R, ZTOP - 7, ZTOP))
for b in P.PI_HOLES:
    bosses.append(cyl_z(*b, 5.5, ZTOP - 12, ZTOP))
for sx in (1, -1):
    bosses.append(box(sx * P.STRAP_X - 8, sx * P.STRAP_X + 8, -17, 17, Z0, Z0 + 6))
SEAM_PTS = [(-76, Z0 + 4.5), (76, Z0 + 4.5), (-76, ZTOP - 4.5), (76, ZTOP - 4.5),
            (-38, ZTOP - 4.5), (38, ZTOP - 4.5), (0, ZTOP - 4.5),
            (-38, Z0 + 4.5), (38, Z0 + 4.5)]
for sx in (1, -1):                                 # seam flange rings + lugs
    xs = sx * P.SEAM_X
    bosses.append(box(xs - 6, xs + 6, -84, -70, Z0, ZTOP))
    bosses.append(box(xs - 6, xs + 6, 70, 84, Z0, ZTOP))
    bosses.append(box(xs - 6, xs + 6, -84, 84, ZTOP - 9, ZTOP))
    bosses.append(box(xs - 6, xs + 6, -84, 84, Z0, Z0 + 9))
fb = np.maximum.reduce(bosses)
f = np.maximum(f, fb)
f = gaussian_filter(f, 1.0)                        # fillet bosses into members

# keep-outs and envelope
inplan = np.hypot(np.maximum(np.abs(FX) - A_, 0), FY) - P.DECK_R
keep = (np.abs(FX) < 64) & (np.abs(FY) < 44) & (FZ < Z0)
keep |= (FX > -72) & (FX < 50) & (np.abs(FY) > 40) & (FZ < Z0)
keep |= (np.abs(FY) > 78) & (np.hypot(FX, FZ - P.WHEEL_R) < 60)
for xc in (P.CASTER_X, -P.CASTER_X):
    keep |= (np.hypot(FX - xc, FY) < 26) & (FZ < Z0)
f[keep] = 0
f[inplan > -0.5] = 0
f[FZ > ZTOP] = 0
f[FZ < ZB] = 0
LEVEL = 0.42

verts, faces, _, _ = marching_cubes(f, LEVEL, spacing=(S, S, S))
verts += np.array([fx[0], fy[0], fz[0]])
tm = trimesh.Trimesh(verts, faces, process=True)
trimesh.smoothing.filter_taubin(tm, iterations=8)
if tm.volume < 0:
    tm.invert()
parts = tm.split(only_watertight=True)
parts = [p for p in parts if p.volume > 2000]      # drop floating crumbs
print("bodies kept:", len(parts), [round(p.volume) for p in parts])
tm = trimesh.util.concatenate(parts)
mesh = Mesh(vert_properties=np.asarray(tm.vertices, np.float32),
            tri_verts=np.asarray(tm.faces, np.uint32))
chassis = Manifold(mesh)
print("manifold status:", chassis.status(), "vol cm3:", round(chassis.volume() / 1000))

# ---------------- exact features ----------------
def cylz(r, z0, z1, x, y, seg=32):
    return Manifold.cylinder(z1 - z0, r, r, seg).translate((x, y, z0))

def cylx(r, x0, x1, y, z, seg=24):
    return Manifold.cylinder(x1 - x0, r, r, seg).rotate((0, 90, 0)).translate((x0, y, z))

def cube(x0, x1, y0, y1, z0, z1):
    return Manifold.cube((x1 - x0, y1 - y0, z1 - z0)).translate((x0, y0, z0))

def U(ms):
    return Manifold.batch_boolean(ms, OpType.Add)

cut = []
cut.append(cube(-300, 300, -300, 300, ZTOP, 300))               # flat top
cut.append(cube(-64, 64, -44, 44, -100, Z0))                    # battery bay, flat ceiling
for s_ in (1, -1):                                              # pod mounting faces
    cut.append(cube(-72, 50, 40 if s_ > 0 else -300, 300 if s_ > 0 else -40, -100, Z0))
for xc in (P.CASTER_X, -P.CASTER_X):
    cut.append(cylz(26, -100, Z0, xc, 0, 64))
for b in pads:                                                 # M3 heat-set inserts
    cut.append(cylz(2.0, Z0 - 1, Z0 + 6, *b))
for xc in (P.CASTER_X, -P.CASTER_X):
    for a in (0, 120, 240):
        x = xc + P.CASTER_BOLT_R * math.cos(math.radians(a))
        y = P.CASTER_BOLT_R * math.sin(math.radians(a))
        cut.append(cylz(2.0, Z0 - 1, Z0 + 6, x, y))
for (hx, hy) in P.hub_holes():                                 # arm pad holes
    cut.append(cylz(P.M3 / 2, ZTOP - 30, ZTOP + 1, hx, hy))
cut.append(cylz(13, ZB - 10, ZTOP + 1, P.HUB_X, 0, 64))        # cable pass-through
for b in P.PI_HOLES:
    cut.append(cylz(P.M25 / 2, ZTOP - 9, ZTOP + 1, *b))
for sx in (1, -1):                                             # battery strap slots
    for dx in (-4.5, 4.5):
        cut.append(cube(sx * P.STRAP_X + dx - 1.6, sx * P.STRAP_X + dx + 1.6, -11, 11, Z0 - 1, Z0 + 7))
for sx in (1, -1):                                             # seam bolts along x
    xs = sx * P.SEAM_X
    for (y, z) in SEAM_PTS:
        cut.append(cylx(P.M3 / 2, xs - 7, xs + 7, y, z))
        cut.append(cylx(3.6, xs + 6, xs + 22, y, z))           # head / nut access
        cut.append(cylx(3.6, xs - 22, xs - 6, y, z))
chassis = chassis - U(cut)

pieces = {
    "chassis_center": chassis ^ cube(-P.SEAM_X + 0.1, P.SEAM_X - 0.1, -200, 200, -100, 300),
    "chassis_front": chassis ^ cube(P.SEAM_X + 0.1, 300, -200, 200, -100, 300),
    "chassis_rear": chassis ^ cube(-300, -P.SEAM_X - 0.1, -200, 200, -100, 300),
}
os.makedirs("stl/print", exist_ok=True); os.makedirs("stl/assembly", exist_ok=True)
report = []
for name, m in pieces.items():
    out = m.to_mesh()
    t = trimesh.Trimesh(np.asarray(out.vert_properties)[:, :3], np.asarray(out.tri_verts))
    t.export(f"stl/assembly/{name}.stl")
    lo, hi = t.bounds
    t2 = t.copy(); t2.apply_translation([-(lo[0] + hi[0]) / 2, -(lo[1] + hi[1]) / 2, -lo[2]])
    t2.export(f"stl/print/{name}.stl")
    ext = (hi - lo).round(1).tolist()
    fits = all(e <= P.BED for e in ext)
    report.append(dict(name=name, qty=1, material="PETG / PLA+", size_mm=ext, fits_a1_mini=fits,
                       volume_cm3=round(m.volume() / 1000, 1), est_grams=round(m.volume() / 1000 * 1.27 * 0.6),
                       note="Print upright, tree supports (auto). 3 walls, 15% gyroid. Heat-set M3 inserts in the bottom pads.",
                       watertight=bool(t.is_watertight), bodies=len(t.split(only_watertight=False))))
    print(f"{name:16s} {ext} fits={fits} {m.volume()/1000:.0f} cm3 watertight={t.is_watertight}")
json.dump(report, open("stl/chassis.json", "w"), indent=1)

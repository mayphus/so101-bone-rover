"""Optimized density -> one smooth, bone-like structure -> 4 printable pieces.

1. upsample + blur the 3D density into smooth members
2. union the exact mounting features (flat plates, pads, bosses)
3. exact cuts: motor bores, plug space, wheel wells, battery bay, all holes
4. split for the A1 mini: rear body | front body at x = SPLIT_X, motor pods at |y| = SPLIT_Y.
   Wherever a member crosses a cut it gets a collar (its own section grown 9 mm,
   10 mm long) with M3 screws into heat-set inserts. No walls anywhere else.
usage: structure.py <npz> <variant> [level]
"""
import json, os, sys
import numpy as np
import trimesh
from manifold3d import Manifold, Mesh, OpType, CrossSection, FillRule, JoinType
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from skimage.measure import marching_cubes
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
import params as P
import layout as L

SRC = sys.argv[1] if len(sys.argv) > 1 else "topo_print.npz"
VAR = sys.argv[2] if len(sys.argv) > 2 else "print"
LEVEL = float(sys.argv[3]) if len(sys.argv) > 3 else 0.40
OUTD = f"stl/{VAR}"

d = np.load(SRC)
rho, X0, Y0, ZB, H = d["rho"], float(d["x0"]), float(d["y0"]), float(d["z0"]), float(d["h"])
nx, ny, nz = rho.shape

S = float(os.environ.get("ISO_S", "1.0"))
cx = X0 + H * (np.arange(nx) + 0.5)
cy = Y0 + H * (np.arange(ny) + 0.5)
cz = ZB + H * (np.arange(nz) + 0.5)
YM = P.F + P.PLATE_T + 6
fx = np.arange(P.X_TAIL - 6, P.X_NOSE + 6 + S, S)
fy = np.arange(-YM, YM + S, S)
fz = np.arange(P.Z_FLOOR - 6, P.Z_TOP + 6, S)
FX, FY, FZ = np.meshgrid(fx, fy, fz, indexing="ij")
interp = RegularGridInterpolator((cx, cy, cz), rho, bounds_error=False, fill_value=0.0)
f = interp(np.stack([FX.ravel(), FY.ravel(), FZ.ravel()], 1)).reshape(FX.shape)
# The light mast and board loads leave half-dense struts (rho ~0.5-0.6) behind the motors.
# Firm up grey material there so those struts survive the isosurface instead of vanishing.
alpha = np.clip((-14.0 - FX) / 20.0, 0, 1)
f = np.clip(f, 0, 1) ** (1 - 0.55 * alpha)
f = gaussian_filter(f, 2.4)
f /= max(f.max(), 1e-9)
feat = gaussian_filter(L.solid(FX, FY, FZ).astype(float), 1.2)
f = gaussian_filter(np.maximum(f, feat), 1.0)          # fillet features into members

# ---- FEA-guided clean-up (needs stl/<variant>/stress_field.npz from verify.py) ----
# Material thinner than 2*R_MIN is either a dead-end fin (low stress: the 45-degree overhang rule
# grows these to avoid supports) or a thin load-carrying bridge. Fins go, bridges thicken.
R_MIN = 2.5
SF_PATH = f"{OUTD}/stress_field.npz"
if os.environ.get("REFINE") and os.path.exists(SF_PATH):
    from scipy.ndimage import distance_transform_edt as _edt0
    sfd = np.load(SF_PATH)
    env, sx0, sh = sfd["env"], sfd["x0"], float(sfd["h"])
    ii = [np.clip(np.round((g - sx0[k]) / sh - 0.5).astype(int), 0, env.shape[k] - 1) for k, g in enumerate((FX, FY, FZ))]
    st = env[ii[0], ii[1], ii[2]]
    from scipy.ndimage import maximum_filter
    st = maximum_filter(st, size=7)                       # stress of the nearby material (~7 mm)
    solid_now = f > LEVEL
    core = _edt0(solid_now) * S > R_MIN
    thick = (_edt0(~core) * S) <= R_MIN + S
    thin = solid_now & ~thick
    featm = L.solid(FX, FY, FZ)
    near_feat = (_edt0(~featm) * S) < 1.5          # features come back exactly afterwards
    joint_zone = (np.abs(FX - P.SPLIT_X) < P.COLLAR_T + 8) | \
                 ((np.abs(np.abs(FY) - P.SPLIT_Y) < P.COLLAR_T + 8) & (np.abs(FX) < P.POD_X + 8))
    # dead material: under 4% of the frame's peak stress in every load case (fins, tails, stubs)
    core7 = _edt0(solid_now) * S > 3.5
    blade = solid_now & ~((_edt0(~core7) * S) <= 3.5 + S)          # thinner than ~7 mm
    fins = ((solid_now & (st < 0.04)) | (blade & (st < 0.15))) & ~near_feat & ~joint_zone
    bridges = thin & (st >= 0.12)
    print(f"refine: thin {thin.sum()*S**3/1000:.1f} cm3 -> fins removed {fins.sum()*S**3/1000:.1f} cm3, bridges thickened {bridges.sum()*S**3/1000:.1f} cm3")
    f = np.where(fins, np.minimum(f, 0.25 * LEVEL), f)
    grow = gaussian_filter(((_edt0(~bridges) * S) <= R_MIN).astype(float), 1.0)
    f = np.maximum(f, grow)
    f = gaussian_filter(f, 0.8)

# ---- joints grown into the bone (not bolted-on blocks) ----
# Where a member crosses a print cut, the member swells into a smooth knuckle: full width
# (section + COLLAR_GROW) for COLLAR_T each side of the cut, then a rounded taper back into
# the member. It is part of the same field, so it fillets in like the rest of the structure.
from scipy.ndimage import distance_transform_edt
TAPER = 14.0
def grow_at(dw):
    t = np.clip((dw - P.COLLAR_T) / TAPER, 0, 1)
    return P.COLLAR_GROW * np.cos(np.pi / 2 * t)
knuckle = np.zeros(f.shape, bool)
for axis, val in (("x", P.SPLIT_X), ("y", P.SPLIT_Y), ("y", -P.SPLIT_Y)):
    ax_c = fx if axis == "x" else fy
    i0 = int(np.argmin(np.abs(ax_c - val)))
    sec2 = (f[i0] if axis == "x" else f[:, i0, :]) > LEVEL
    if axis == "y":                                        # only the motor pod members
        sec2 &= (np.abs(fx) < P.POD_X + 6)[:, None]
    dist = distance_transform_edt(~sec2) * S
    for i, w in enumerate(ax_c):
        dw = abs(w - val)
        if dw > P.COLLAR_T + TAPER:
            continue
        reg = dist <= grow_at(dw) + 1e-6
        if axis == "x":
            knuckle[i] |= reg
        else:
            knuckle[:, i, :] |= reg
f = np.maximum(f, gaussian_filter(knuckle.astype(float), 1.6))
# Soft trim: fade the field out over ~1.2 mm before every keep-out, so members end in a rounded
# edge just short of a cut instead of being sliced into a thin, feathered lip by the exact cut.
from scipy.ndimage import distance_transform_edt as _edt
KO = L.keepout(FX, FY, FZ)
dko = _edt(~KO) * S
f = np.minimum(f, LEVEL * (0.55 + dko / 1.2))
f[KO] = 0.0

verts, faces, _, _ = marching_cubes(f, LEVEL, spacing=(S, S, S))
verts += np.array([fx[0], fy[0], fz[0]])
tm = trimesh.Trimesh(verts, faces, process=True)
trimesh.smoothing.filter_taubin(tm, iterations=30, lamb=0.5, nu=-0.53)
if tm.volume < 0:
    tm.invert()
bodies = [b for b in tm.split(only_watertight=True) if abs(b.volume) > 3000]
print("bodies kept:", [round(b.volume) for b in bodies])
N_ISLANDS = len(bodies)
tm = trimesh.util.concatenate(bodies)
body = Manifold(Mesh(vert_properties=np.asarray(tm.vertices, np.float32),
                     tri_verts=np.asarray(tm.faces, np.uint32)))
print("isosurface:", body.status(), round(body.volume() / 1000), "cm3")

# ---------------- exact geometry helpers ----------------
def U(ms):
    ms = [m for m in ms if not m.is_empty()]
    return Manifold.batch_boolean(ms, OpType.Add) if ms else Manifold()

def cube(x0, x1, y0, y1, z0, z1):
    return Manifold.cube((x1 - x0, y1 - y0, z1 - z0)).translate((x0, y0, z0))

def cylz(r, z0, z1, x, y, seg=32):
    return Manifold.cylinder(z1 - z0, r, r, seg).translate((x, y, z0))

def cylx(r, x0, x1, y, z, seg=24):
    return Manifold.cylinder(x1 - x0, r, r, seg).rotate((0, 90, 0)).translate((x0, y, z))

def cyly(r, y0, y1, x, z, seg=64):
    return Manifold.cylinder(y1 - y0, r, r, seg).rotate((-90, 0, 0)).translate((x, y0, z))

def extrude_poly(poly, z0, z1):
    rings = [np.asarray(poly.exterior.coords)[:-1]] + [np.asarray(r.coords)[:-1] for r in poly.interiors]
    return Manifold.extrude(CrossSection(rings, FillRule.EvenOdd), z1 - z0).translate((0, 0, z0))

def rrect(x0, x1, hw, r):
    return box(x0 + r, -hw + r, x1 - r, hw - r).buffer(r, 16)

plan = unary_union([rrect(P.X_TAIL, P.X_NOSE, P.HALF_W, P.CORNER_R),
                    rrect(-P.POD_X, P.POD_X, P.F + P.PLATE_T, 8.0)]).buffer(-0.4)
plan_solid = extrude_poly(plan, P.Z_FLOOR - 1, P.Z_TOP + 30)
m_in = P.F - P.MOTOR_TOTAL_L

# ---- preserved geometry: exact functional features ----
keep = [extrude_poly(L.ARM_PAD, P.Z_TOP - 6, P.Z_TOP)]
for hx, hy in P.ARM_HOLES:
    keep.append(cylz(7, P.Z_TOP - 16, P.Z_TOP, hx, hy, 40))
for hx, hy in P.DRV_HOLES:                                       # driver bosses; over the caster pad they
    keep.append(cylz(5.5, P.Z_TOP - 10, P.Z_TOP, hx, hy, 32))
for s_ in (1, -1):
    y0, y1 = sorted((s_ * P.F, s_ * (P.F + P.PLATE_T)))
    keep.append(cyly(24, y0, y1, 0, P.MOTOR_AXIS_Z, 96))                     # flat motor face plate
    y0, y1 = sorted((s_ * (m_in + 4), s_ * (m_in + 14)))
    keep.append(cyly(P.MOTOR_CLEAR_R + 5, y0, y1, 0, P.MOTOR_AXIS_Z, 96))     # rear cradle ring
for hx, hy in P.CASTER_HOLES:                                     # flat-bottomed caster bolt bosses
    keep.append(cylz(P.CASTER_BOSS, P.CASTER_H, P.CASTER_H + 10, P.CASTER_X + hx, hy, 40))
for sx in P.STRAP_X:
    keep.append(cube(sx - P.STRAP_HALF, sx + P.STRAP_HALF, -18, 18, P.DECK_Z, P.DECK_Z + 6))
keep.append(cylz(P.MAST_D / 2 + 6, P.Z_TOP - P.MAST_SOCKET, P.Z_TOP, P.MAST_X, 0, 64))
body = body + (U(keep) ^ plan_solid)

cut = [cube(-400, 400, -400, 400, P.Z_TOP, 500), cube(-400, 400, -400, 400, -100, P.Z_FLOOR - 0.01)]
for s in (1, -1):
    ys = lambda a, b: (min(s * a, s * b), max(s * a, s * b))
    y0, y1 = ys(m_in - 1, P.F)                                    # motor bore, flat face seat
    cut.append(cyly(P.MOTOR_CLEAR_R, y0, y1, 0, P.MOTOR_AXIS_Z))
    y0, y1 = ys(0, m_in + 0.5)                                    # encoder plug + cable
    cut.append(cyly(13, y0, y1, P.PLUG_X, P.MOTOR_AXIS_Z + P.PLUG_Z, 48))
    y0, y1 = ys(P.F + P.PLATE_T, 300)                             # wheel well
    cut.append(cyly(P.WHEEL_R + 6, y0, y1, 0, P.AXLE_Z))
    y0, y1 = ys(P.F - 2, 300)                                     # boss + coupling
    cut.append(cyly(P.BOSS_D / 2 + 0.6, y0, y1, 0, P.AXLE_Z, 48))
    for hx, hz in P.face_holes():                                 # 6 x M3 face screws
        y0, y1 = ys(P.F - 2, P.F + P.PLATE_T + 10)
        cut.append(cyly(P.M3 / 2, y0, y1, hx, hz, 20))
b = P.BAY
cut.append(cube(b["x0"], b["x1"], -b["y"], b["y"], -50, b["z1"]))
cut.append(cube(-4, b["x0"] + 1, -9, 9, 52, 76))                  # motor cables -> bay
cut.append(cylz(10, P.Z_TOP - 40, P.Z_TOP + 1, *P.ARM_CABLE, 64))
for hx, hy in P.ARM_HOLES:                                        # M5 heat-set inserts
    cut.append(cylz(P.M5_INSERT_D / 2, P.Z_TOP - P.M5_INSERT_L, P.Z_TOP + 1, hx, hy, 32))
    cut.append(cylz(2.7, P.Z_TOP - 20, P.Z_TOP + 1, hx, hy, 24))
for hx, hy in P.DRV_HOLES:                                        # M3 inserts for driver standoffs
    cut.append(cylz(P.M3_INSERT_D / 2, P.Z_TOP - P.M3_INSERT_L, P.Z_TOP + 1, hx, hy, 24))
for sx in P.STRAP_X:
    for dx in (-4.5, 4.5):
        cut.append(cube(sx + dx - 1.6, sx + dx + 1.6, -11, 11, P.DECK_Z - 2, P.DECK_Z + 8))
cut.append(cylz(P.CASTER_SWING, -50, P.CASTER_H, P.CASTER_X, 0, 96))
for hx, hy in P.CASTER_HOLES:                                     # M4 inserts from below
    cut.append(cylz(P.M4_INSERT_D / 2, P.CASTER_H - 1, P.CASTER_H + P.M4_INSERT_L, P.CASTER_X + hx, hy, 24))
mast_bore = cylz(P.MAST_D / 2 + 0.15, P.Z_TOP - P.MAST_SOCKET + 4, P.Z_TOP + 40, P.MAST_X, 0, 48)
cut.append(mast_bore)
CUT = U(cut)
body = body - CUT

# ---------------- collars where members cross a print cut ----------------
# rotate so the cut plane's normal is +z, slice there, grow each section into a collar
PERM = {"x": (np.array([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]], float),
              np.array([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0]], float)),
        "y": (np.array([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0]], float),
              np.array([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]], float))}

def to_world(axis, u, v, w):
    return (w, u, v) if axis == "x" else (v, w, u)

collars, holes, joints = [], [], []
for axis, val, head_sign in (("x", P.SPLIT_X, 1), ("y", P.SPLIT_Y, 1), ("y", -P.SPLIT_Y, -1)):
    fwd, back = PERM[axis]
    sec = body.transform(fwd).slice(val)
    for comp in sec.decompose():
        if comp.area() < 25:
            continue
        rings = comp.to_polygons()
        outer = Polygon(max(rings, key=lambda r: Polygon(r).area))
        inner = outer.buffer(-(P.COLLAR_GROW - 4.7), 16)    # knuckle edge -> screw circle
        if inner.is_empty:
            continue
        if inner.geom_type == "MultiPolygon":
            inner = max(inner.geoms, key=lambda g_: g_.area)
        ring = inner.exterior
        want = int(np.clip(round(ring.length / 40), 2, 6))
        nc = 36
        cand = [ring.interpolate(t, normalized=True) for t in (np.arange(nc) + 0.5) / nc]
        n_ok, last = 0, -99
        for ci, p in enumerate(cand):
            if n_ok >= want or ci - last < nc // (want + 1):
                continue
            u, v = p.x, p.y
            wx, wy, wz = to_world(axis, u, v, val)
            if axis == "y" and np.hypot(wx, wz - P.MOTOR_AXIS_Z) < 27:     # stay off the motor face plate
                continue
            if not (P.Z_FLOOR + 3 < wz < P.Z_TOP - 3) or L.plan_sd(np.array([wx]), np.array([wy]))[0] > -3:
                continue
            hs = head_sign
            parts = [Manifold.cylinder(2 * P.COLLAR_T + 2, P.M3 / 2, P.M3 / 2, 20).translate((u, v, val - P.COLLAR_T - 1)),
                     Manifold.cylinder(14, 3.1, 3.1, 24).translate((u, v, val + P.COLLAR_T if hs > 0 else val - P.COLLAR_T - 14)),
                     Manifold.cylinder(P.M3_INSERT_L + 0.5, P.M3_INSERT_D / 2, P.M3_INSERT_D / 2, 24)
                     .translate((u, v, val - P.M3_INSERT_L - 0.5 if hs > 0 else val))]
            holes.append(U(parts).transform(back))
            n_ok += 1; last = ci
        joints.append(dict(plane=f"{axis}={val:+.0f}", section_mm2=round(comp.area()), screws=n_ok))
print("collars:", joints)
body = (body - CUT) - U(holes)
# mast socket rim above the deck with an M3 grub screw
rim = cylz(P.MAST_D / 2 + 5, P.Z_TOP - 0.5, P.Z_TOP + 14, P.MAST_X, 0, 64) - mast_bore
rim -= cylx(P.M3_TAP / 2, P.MAST_X - 20, P.MAST_X, 0, P.Z_TOP + 7, 16)
body = body + rim

g = 0.1
pieces = {
    "frame_rear": body ^ cube(-400, P.SPLIT_X - g, -P.SPLIT_Y + g, P.SPLIT_Y - g, -100, 400),
    "frame_front": body ^ cube(P.SPLIT_X + g, 400, -400, 400, -100, 400),
    "frame_pod_left": body ^ cube(-400, P.SPLIT_X - g, P.SPLIT_Y + g, 400, -100, 400),
    "frame_pod_right": body ^ cube(-400, P.SPLIT_X - g, -400, -P.SPLIT_Y - g, -100, 400),
}
NOTES = {
    "frame_rear": "Upright as on the robot. Tree supports (auto). M5 inserts in the arm pad, M3 inserts in the knuckle faces.",
    "frame_front": "Stand it on its flat cut face. Tree supports (auto). M4 inserts in the caster bosses, M3 inserts for the driver.",
    "frame_pod_left": "Lay the motor face plate flat on the bed. M3 inserts in the knuckle faces.",
    "frame_pod_right": "Mirror of the left pod. Lay the face plate flat on the bed.",
}
os.makedirs(f"{OUTD}/print", exist_ok=True)
os.makedirs(f"{OUTD}/assembly", exist_ok=True)
report = []
for name, m in pieces.items():
    out = m.to_mesh()
    # export manifold3d's mesh as is: it is closed by construction (merging/cleanup in trimesh
    # can open tiny gaps); closure is checked with manifold3d itself below
    t = trimesh.Trimesh(np.asarray(out.vert_properties)[:, :3], np.asarray(out.tri_verts), process=False)
    t.export(f"{OUTD}/assembly/{name}.stl")
    t2 = t.copy()
    if name == "frame_front":                           # flat mating face down: half the supports
        t2.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [0, 1, 0]))
    if "pod" in name:                                   # face plate down on the bed
        t2.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2 * (1 if "left" in name else -1), [1, 0, 0]))
    lo, hi = t2.bounds
    t2.apply_translation([-(lo[0] + hi[0]) / 2, -(lo[1] + hi[1]) / 2, -lo[2]])
    t2.export(f"{OUTD}/print/{name}.stl")
    fz_ = t2.face_normals[:, 2]; cz_ = t2.triangles_center[:, 2]
    over = float(t2.area_faces[(fz_ < -0.7071) & (cz_ > 0.6)].sum())
    ext = (t2.bounds[1] - t2.bounds[0]).round(1).tolist()
    fits = all(e <= P.BED for e in ext)
    real = [p for p in t.split(only_watertight=False) if abs(p.volume) > 50]
    report.append(dict(name=name, qty=1, material="PETG (or PLA+)", size_mm=ext, fits_a1_mini=fits,
                       volume_cm3=round(m.volume() / 1000, 1),
                       est_grams=round(m.volume() / 1000 * 1.27 * 0.55),
                       note=NOTES[name], watertight=bool(m.status().name == 'NoError' and not m.is_empty() and m.genus() >= 0), bodies=len(real),
                       overhang_cm2=round(over / 100, 1)))
    print(f"{name:16s} {ext} fits={fits} {m.volume()/1000:.0f} cm3 watertight={t.is_watertight} bodies={len(real)}")
json.dump(report, open(f"{OUTD}/frame.json", "w"), indent=1)
summary = dict(variant=VAR, compliance=float(d["c"]) if "c" in d else None,
               grams=sum(r["est_grams"] for r in report), volume_cm3=round(sum(r["volume_cm3"] for r in report), 1),
               overhang_cm2=round(sum(r["overhang_cm2"] for r in report), 1),
               fits=all(r["fits_a1_mini"] for r in report), bodies=[r["bodies"] for r in report],
               islands=N_ISLANDS, joints=joints, screws=sum(j["screws"] for j in joints))
json.dump(summary, open(f"{OUTD}/summary.json", "w"), indent=1)
print(summary)

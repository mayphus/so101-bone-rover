"""FEA check of a finished frame variant with real loads and PETG.

The smoothed frame (4 pieces, bolted collars treated as bonded) is voxelized
at 2.5 mm and solved as linear-elastic hex8.
Supports: each wheel is a hinge about its axle; the swivel caster pad stands on the floor.
Writes stl/<variant>/fea.json and per-vertex stress for the viewer.
usage: verify.py <variant>
"""
import json, sys, time
import numpy as np
import pyamg
import shapely
from shapely.geometry import Polygon
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg
import trimesh
from manifold3d import Manifold, Mesh, OpType
import params as P
import layout as L

VAR = sys.argv[1]
H = 2.5                                     # voxel size, mm
E, NU = 2000.0, 0.35                        # PETG, MPa
ALLOW = 25.0                                # MPa, conservative for FDM layer bonding
G = 9.81
ARM_KG, PAYLOAD_KG, BATT_KG, MOTOR_KG, PI_KG = 0.75, 0.2, 0.5, 0.35, 0.15
MAST_KG, MAST_ARM = 0.12, 380.0          # tube + head + camera, lever to the socket
DYN = 2.0                                   # dynamic factor on arm cases

# ---------------- voxelize ----------------
pieces = []
for n in P.PIECES:
    t = trimesh.load(f"stl/{VAR}/assembly/{n}.stl", process=True)
    pieces.append(Manifold(Mesh(vert_properties=np.asarray(t.vertices, np.float32), tri_verts=np.asarray(t.faces, np.uint32))))
frame = Manifold.batch_boolean(pieces, OpType.Add)
lo, hi = np.array(frame.bounding_box()[:3]), np.array(frame.bounding_box()[3:])
X0 = np.floor(lo / H) * H - H * 0.37        # offset so voxel centres avoid the 0.2 mm cut gaps
n = np.ceil((hi - X0) / H).astype(int) + 1
xc = X0[0] + H * (np.arange(n[0]) + 0.5)
yc = X0[1] + H * (np.arange(n[1]) + 0.5)
zc = X0[2] + H * (np.arange(n[2]) + 0.5)
GX, GY = np.meshgrid(xc, yc, indexing="ij")
vox = np.zeros(n, bool)
for k, z in enumerate(zc):
    cs = frame.slice(z)
    polys = cs.to_polygons()
    if not polys:
        continue
    # even-odd fill of all contours
    inside = np.zeros(GX.shape, bool)
    for ring in polys:
        if len(ring) >= 3:
            inside ^= shapely.contains_xy(Polygon(ring), GX, GY)
    vox[:, :, k] = inside
from scipy.ndimage import label as _label
_lab, _nl = _label(vox, structure=np.ones((3, 3, 3)))
_sz = np.bincount(_lab.ravel()); _sz[0] = 0
DROPPED = int(vox.sum() - _sz.max())
vox = _lab == _sz.argmax()                  # drop voxel specks thinner than the 2.5 mm grid can connect
print("voxels", vox.sum(), "dropped specks", DROPPED, "grid", n.tolist(), flush=True)
if DROPPED > 40:
    raise SystemExit(f"FRAME NOT CONNECTED: {DROPPED} voxels ({DROPPED * H**3 / 1000:.1f} cm3) are separate from the main body")

# ---------------- FE model ----------------
def hex8_ke(nu):
    C = 1 / ((1 + nu) * (1 - 2 * nu)) * np.array([
        [1-nu, nu, nu, 0, 0, 0], [nu, 1-nu, nu, 0, 0, 0], [nu, nu, 1-nu, 0, 0, 0],
        [0, 0, 0, (1-2*nu)/2, 0, 0], [0, 0, 0, 0, (1-2*nu)/2, 0], [0, 0, 0, 0, 0, (1-2*nu)/2]])
    nodes = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]], float)
    def Bmat(a, b, c):
        dN = np.zeros((8, 3))
        for i, (nx_, ny_, nz_) in enumerate(nodes):
            sx = 1 if nx_ else -1; sy = 1 if ny_ else -1; sz = 1 if nz_ else -1
            fx = a if nx_ else 1 - a; fy = b if ny_ else 1 - b; fz = c if nz_ else 1 - c
            dN[i] = [sx * fy * fz, sy * fx * fz, sz * fx * fy]
        B = np.zeros((6, 24))
        for i in range(8):
            B[0, 3*i] = dN[i, 0]; B[1, 3*i+1] = dN[i, 1]; B[2, 3*i+2] = dN[i, 2]
            B[3, 3*i] = dN[i, 1]; B[3, 3*i+1] = dN[i, 0]
            B[4, 3*i+1] = dN[i, 2]; B[4, 3*i+2] = dN[i, 1]
            B[5, 3*i] = dN[i, 2]; B[5, 3*i+2] = dN[i, 0]
        return B
    g = [0.5 - 0.5/np.sqrt(3), 0.5 + 0.5/np.sqrt(3)]
    K = sum(Bmat(a, b, c).T @ C @ Bmat(a, b, c) / 8.0 for a in g for b in g for c in g)
    return K, C, Bmat(0.5, 0.5, 0.5), nodes.astype(int)

KEu, Cu, Bc, corners = hex8_ke(NU)
KE = E * H * KEu                             # scale unit cube to size H, modulus E
ei, ej, ek = np.nonzero(vox)
def nid(i, j, k):
    return (i * (n[1] + 1) + j) * (n[2] + 1) + k
enodes = np.stack([nid(ei + c[0], ej + c[1], ek + c[2]) for c in corners], 1)
used = np.unique(enodes)
remap = np.full((n[0] + 1) * (n[1] + 1) * (n[2] + 1), -1); remap[used] = np.arange(len(used))
en = remap[enodes]
edof = np.stack([3 * en + d for d in range(3)], 2).reshape(len(ei), 24)
ndof = 3 * len(used)
K = coo_matrix((np.tile(KE.ravel(), len(ei)), (np.repeat(edof, 24, 1).ravel(), np.tile(edof, (1, 24)).ravel())),
               shape=(ndof, ndof)).tocsr()
ii, jj, kk = np.unravel_index(used, n + 1)
NX, NY, NZ = X0[0] + H * ii, X0[1] + H * jj, X0[2] + H * kk
print("elements", len(ei), "dofs", ndof, flush=True)

def sel(mask):
    return np.nonzero(mask)[0]

# supports: wheel hinge (nodes on each face plate around the shaft hole), casters on the floor
plate = (np.abs(NY) >= P.F - 0.5) & (np.abs(NY) <= P.F + P.PLATE_T + 0.5)
dz, dxn = NZ - P.AXLE_Z, NX
ring = plate & (np.hypot(dxn, dz) > P.BOSS_D / 2 + 0.5) & (np.hypot(dxn, dz) < 12)
nA = sel(ring & (np.abs(dxn) < H * 0.6))          # above / below shaft: hold y, z
nB = sel(ring & (np.abs(dz) < H * 0.6))           # beside shaft: hold x
hinge = [3 * nA + 1, 3 * nA + 2, 3 * nB]
cnodes = sel((np.min([np.hypot(NX - P.CASTER_X - hx, NY - hy) for hx, hy in P.CASTER_HOLES], 0) < P.CASTER_BOSS)
              & (np.abs(NZ - P.CASTER_H) < H * 0.6))
print("hinge nodes", len(nA), len(nB), "caster nodes", len(cnodes), flush=True)
FREE = np.setdiff1d(np.arange(ndof), np.unique(np.concatenate(hinge + [3 * cnodes + 2])))

# ---------------- load cases (N, mm) ----------------
pad = sel(L.in_pad(NX, NY) & (NZ > P.Z_TOP - H * 0.6))
dxp, dyp = NX[pad] - P.HUB_X, NY[pad]
def arm_case(W, Mx=0.0, My=0.0):
    f = np.zeros(ndof)
    f[3 * pad + 2] += -W / len(pad)
    if My: f[3 * pad + 2] += -My * dxp / (dxp ** 2).sum()     # nose-down pitch
    if Mx: f[3 * pad + 2] += -Mx * dyp / (dyp ** 2).sum()     # roll
    return f
W_arm = (ARM_KG + PAYLOAD_KG) * G
M_reach = (ARM_KG * G * 150 + PAYLOAD_KG * G * 400)            # N*mm at full reach
cases = {
    "Arm reaching forward, 200 g payload, 2x dynamic": (arm_case(DYN * W_arm, My=DYN * M_reach), "fixed"),
    "Arm reaching sideways, 200 g payload, 2x dynamic": (arm_case(DYN * W_arm, Mx=DYN * M_reach), "fixed"),
}
f = arm_case(3 * W_arm)
strap = sel((np.abs(NZ - P.DECK_Z) < H * 0.6) & np.any([np.abs(NX - sx) < 7 for sx in P.STRAP_X], 0) & (np.abs(NY) < 18))
f[3 * strap + 2] -= 3 * BATT_KG * G / max(len(strap), 1)
m_in = P.F - P.MOTOR_TOTAL_L
crad = sel((np.abs(NY) > m_in + 4) & (np.abs(NY) < m_in + 14) & (np.abs(np.hypot(NX, NZ - P.MOTOR_AXIS_Z) - P.MOTOR_CLEAR_R - 1) < 2))
f[3 * crad + 2] -= 3 * 2 * MOTOR_KG * G / max(len(crad), 1)
pi = sel(np.any([np.hypot(NX - hx, NY - hy) < 6 for hx, hy in P.DRV_HOLES], 0) & (NZ > P.Z_TOP - H * 0.6))
f[3 * pi + 2] -= 3 * PI_KG * G / max(len(pi), 1)
mast = sel((np.hypot(NX - P.MAST_X, NY) < P.MAST_D / 2 + 3) & (NZ > P.Z_TOP - P.MAST_SOCKET) & (NZ < P.Z_TOP + 14))
f[3 * mast + 2] -= 3 * MAST_KG * G / max(len(mast), 1)
cases["3 g bump: arm, battery, motors, boards, mast"] = (f, "fixed")
T_STALL = 1000.0                                                # N*mm, JGB37-520 1:90 stall ~ 10 kg*cm
f = np.zeros(ndof)
for s in (1, -1):
    hr = sel(plate & (np.sign(NY) == s) & (np.abs(np.hypot(NX, NZ - P.MOTOR_AXIS_Z) - P.FACE_PCD / 2) < 2.5))
    rx, rz = NX[hr], NZ[hr] - P.MOTOR_AXIS_Z
    r2 = (rx ** 2 + rz ** 2).sum()
    f[3 * hr] += -T_STALL * rz / r2
    f[3 * hr + 2] += T_STALL * rx / r2
cases["Both motors at stall torque"] = (f, "fixed")
# camera mast hit sideways: 3 g on the camera at the top of the tube -> couple on the socket
f = np.zeros(ndof)
top_r = sel((np.hypot(NX - P.MAST_X, NY) < P.MAST_D / 2 + 3) & (NZ > P.Z_TOP + 8))
bot_r = sel((np.hypot(NX - P.MAST_X, NY) < P.MAST_D / 2 + 3) & (NZ < P.Z_TOP - P.MAST_SOCKET + 10) & (NZ > P.Z_TOP - P.MAST_SOCKET))
Fm = 3 * MAST_KG * G
lever = (P.Z_TOP + 11) - (P.Z_TOP - P.MAST_SOCKET + 6)
Ftop = Fm * (MAST_ARM + lever) / lever
f[3 * top_r + 1] += Ftop / max(len(top_r), 1)
f[3 * bot_r + 1] -= (Ftop - Fm) / max(len(bot_r), 1)
cases["Camera mast knocked sideways, 3 g"] = (f, "fixed")
print("pad nodes", len(pad), "strap", len(strap), "cradle", len(crad), "boards", len(pi), "mast", len(mast), len(top_r), len(bot_r), flush=True)

def solver(free):
  Kf = K[free][:, free]
  cx_, cy_, cz_ = NX[free // 3] - NX.mean(), NY[free // 3] - NY.mean(), NZ[free // 3] - NZ.mean()
  d = free % 3
  B = np.zeros((len(free), 6))
  for k in range(3): B[d == k, k] = 1
  B[d == 0, 3] = -cy_[d == 0]; B[d == 1, 3] = cx_[d == 1]
  B[d == 1, 4] = -cz_[d == 1]; B[d == 2, 4] = cy_[d == 2]
  B[d == 0, 5] = cz_[d == 0]; B[d == 2, 5] = -cx_[d == 2]
  ml = pyamg.smoothed_aggregation_solver(Kf, B=B / 100, max_coarse=1000, strength=("symmetric", {"theta": 0.0}))
  return Kf, ml.aspreconditioner()
SOLVERS = {}

CB = E * Cu @ Bc / H                          # stress at element centre per element dof vector
results, envelope = [], np.zeros(len(ei))
for name, (f, sup_key) in cases.items():
    t0 = time.time()
    free = FREE
    if sup_key not in SOLVERS:
        SOLVERS[sup_key] = solver(free)
    Kf, M = SOLVERS[sup_key]
    u = np.zeros(ndof)
    u[free], info = cg(Kf, f[free], rtol=1e-7, maxiter=2000, M=M)
    s = (CB @ u[edof].T).T
    vm = np.sqrt(0.5 * ((s[:, 0]-s[:, 1])**2 + (s[:, 1]-s[:, 2])**2 + (s[:, 2]-s[:, 0])**2) + 3 * (s[:, 3]**2 + s[:, 4]**2 + s[:, 5]**2))
    envelope = np.maximum(envelope, vm)
    disp = np.sqrt(u[0::3]**2 + u[1::3]**2 + u[2::3]**2)
    # pad tilt: best-fit plane through pad node z-displacements
    A_ = np.c_[dxp, dyp, np.ones(len(pad))]
    gx, gy, _ = np.linalg.lstsq(A_, u[3 * pad + 2], rcond=None)[0]
    tilt = np.degrees(np.hypot(gx, gy))
    p995 = float(np.percentile(vm, 99.5))
    r = dict(case=name, max_disp_mm=round(float(disp.max()), 3), pad_tilt_deg=round(float(tilt), 3),
             gripper_error_mm=round(float(np.radians(tilt) * 400), 2),
             stress_p995_mpa=round(p995, 2), stress_max_mpa=round(float(vm.max()), 2),
             safety_factor=round(ALLOW / p995, 1), cg_info=int(info))
    results.append(r)
    print(r, f"{time.time()-t0:.0f}s", flush=True)

summary = dict(variant=VAR, voxel_mm=H, material="PETG E=2000 MPa, allowable 25 MPa", cases=results,
               worst_safety_factor=min(r["safety_factor"] for r in results),
               worst_gripper_error_mm=max(r["gripper_error_mm"] for r in results))
json.dump(summary, open(f"stl/{VAR}/fea.json", "w"), indent=1)

# ---------------- per-vertex stress for the viewer ----------------
SCALE = float(np.percentile(envelope, 99.5))
# stress field on the voxel grid, for the FEA-guided clean-up pass in structure.py
_env = np.zeros(n, np.float32); _env[ei, ej, ek] = envelope / SCALE
np.savez_compressed(f"stl/{VAR}/stress_field.npz", env=_env, x0=X0, h=H)
summary["stress_view_scale_mpa"] = round(SCALE, 3)
json.dump(summary, open(f"stl/{VAR}/fea.json", "w"), indent=1)
from scipy.spatial import cKDTree
cent = np.c_[X0[0] + H * (ei + 0.5), X0[1] + H * (ej + 0.5), X0[2] + H * (ek + 0.5)]
tree = cKDTree(cent)
for nm in P.PIECES:
    j = json.load(open(f"web/{VAR}/assembly/{nm}.json"))
    v = np.asarray(j["v"], float).reshape(-1, 3) * j["s"]
    _, idx = tree.query(v, k=4)
    val = envelope[idx].max(1) / SCALE                  # relative to this frame's own peak (p99.5)
    json.dump({"c": np.round(np.clip(val, 0, 1.0), 3).tolist()}, open(f"web/{VAR}/assembly/{nm}_stress.json", "w"), separators=(",", ":"))
print("worst safety factor", summary["worst_safety_factor"], "worst gripper error mm", summary["worst_gripper_error_mm"])

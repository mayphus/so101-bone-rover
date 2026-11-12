"""3D SIMP topology optimization of the chassis (hex8, multi-load, y-mirrored).

Design space: the stadium plan (330 x 170) from z=30 to z=ZTOP, minus
keep-outs (battery bay, motor pods + wheels, casters). Supports: motor pod
bolt pads (wheels, fully fixed) and caster pads (vertical only).
Loads: SO-101 arm on the top pad (weight, pitch / roll moment at full reach,
yaw twist, side push), battery hanging from the bay ceiling, Pi on top.
Writes topo3d.npz used by build3d.py.
"""
import os, sys, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu, cg
import pyamg
import params as P
import layout as L

H = float(os.environ.get("GRID_H", "6.0"))
ZBOT, ZTOP = P.Z_FLOOR - 2.0, P.Z_TOP
nelx = int(np.ceil((P.X_NOSE - P.X_TAIL + 8) / H))
X0 = P.X_TAIL - 4.0
nely = 2 * int(np.ceil((P.F + P.PLATE_T + 1) / H))
Y0 = -H * nely / 2
nelz = int(round((ZTOP - ZBOT) / H))
xs = X0 + H * (np.arange(nelx) + 0.5)
ys = Y0 + H * (np.arange(nely) + 0.5)
zs = ZBOT + H * (np.arange(nelz) + 0.5)
X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
active = ~L.keepout(X, Y, Z)
solid = L.solid(X, Y, Z) & active
free = active & ~solid
print("elements: active", active.sum(), "solid", solid.sum(), "free", free.sum(), flush=True)

# ---------------- hex8 stiffness (unit cube, E=1, nu=0.3) ----------------
def hex8_ke(nu=0.3):
    E = 1.0
    C = E / ((1 + nu) * (1 - 2 * nu)) * np.array([
        [1-nu, nu, nu, 0, 0, 0], [nu, 1-nu, nu, 0, 0, 0], [nu, nu, 1-nu, 0, 0, 0],
        [0, 0, 0, (1-2*nu)/2, 0, 0], [0, 0, 0, 0, (1-2*nu)/2, 0], [0, 0, 0, 0, 0, (1-2*nu)/2]])
    nodes = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]], float)
    g = [0.5 - 0.5/np.sqrt(3), 0.5 + 0.5/np.sqrt(3)]
    K = np.zeros((24, 24))
    for a in g:
        for b in g:
            for c in g:
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
                K += B.T @ C @ B / 8.0
    return K, nodes.astype(int)

KE, corners = hex8_ke()

def nid(i, j, k):
    return (i * (nely + 1) + j) * (nelz + 1) + k

ei, ej, ek = np.nonzero(active)
eidx_full = np.full(active.shape, -1)
eidx_full[ei, ej, ek] = np.arange(len(ei))
enodes = np.stack([nid(ei + c[0], ej + c[1], ek + c[2]) for c in corners], 1)
used = np.unique(enodes)
remap = np.full((nelx + 1) * (nely + 1) * (nelz + 1), -1)
remap[used] = np.arange(len(used))
en = remap[enodes]
edof = np.stack([3 * en + d for d in range(3)], 2).reshape(len(ei), 24)
ndof = 3 * len(used)
iK = np.repeat(edof, 24, axis=1).ravel()
jK = np.tile(edof, (1, 24)).ravel()

ii, jj, kk = np.unravel_index(used, (nelx + 1, nely + 1, nelz + 1))
NX, NY, NZ = X0 + H * ii, Y0 + H * jj, ZBOT + H * kk
print("dofs", ndof, flush=True)

def near(cx, cy, cz, r, rz=0.1):
    return np.where((np.hypot(NX - cx, NY - cy) < r) & (np.abs(NZ - cz) < rz + 1e-6))[0]

def near3(mask_fn):
    return np.where(mask_fn(NX, NY, NZ))[0]

aNY = np.abs(NY)
# wheels: shaft region of each face plate carries the wheel (fixed)
# wheels spin on their shafts: a hinge about the axle, so pitch must go to the casters.
# Nodes above/below the shaft hole hold y and z; nodes beside it hold x.
zax = ZBOT + H * np.round((P.AXLE_Z - ZBOT) / H)
inplate = lambda y: (np.abs(y) > P.F - 4) & (np.abs(y) < P.F + P.PLATE_T + 4)
nA = near3(lambda x, y, z: inplate(y) & (np.abs(x) < 1) & (np.abs(z - zax) > H * 0.5) & (np.abs(z - zax) < H * 1.5))
nB = near3(lambda x, y, z: inplate(y) & (np.abs(np.abs(x) - H) < 1) & (np.abs(z - zax) < 1))
sup = np.concatenate([nA, nB])
fixed = list(3 * nA + 1) + list(3 * nA + 2) + list(3 * nB)
# the caster's steel top plate presses on every bolt boss: one vertical spring per boss,
# so each boss must be tied into the structure (a rigid support would let the optimizer skip some)
K_BOSS = float(os.environ.get("K_BOSS", "4.0"))
sd_, sv_, ncast = [], [], []
for hx, hy in P.CASTER_HOLES:
    nb_ = near3(lambda x, y, z: (np.hypot(x - P.CASTER_X - hx, y - hy) < P.CASTER_BOSS + 2) & (np.abs(z - P.CASTER_H) < H * 0.75))
    ncast += list(nb_); sd_.append(3 * nb_ + 2); sv_.append(np.full(len(nb_), K_BOSS / max(len(nb_), 1)))
ncast = np.unique(ncast)
fixed = np.unique(fixed)
print("hinge nodes", len(nA), len(nB), "support nodes", len(sup), "caster", len(ncast), "fixed dofs", len(fixed), flush=True)
freedofs = np.setdiff1d(np.arange(ndof), fixed)

top = near3(lambda x, y, z: L.in_pad(x, y) & (z > ZTOP - 0.5))
NC = 13
F = np.zeros((ndof, NC))
dx, dy = NX[top] - P.HUB_X, NY[top]
F[3*top+2, 0] = -1.0                                   # arm weight
F[3*top+2, 1] = -dx / 40.0                             # pitch (reach forward), about the pan axis
F[3*top+2, 2] = -dy / 40.0 - 0.3                       # roll (reach sideways)
rr = np.hypot(dx, dy) + 1e-9
F[3*top, 3] = -dy / rr; F[3*top+1, 3] = dx / rr        # yaw twist
F[3*top, 4] = 1.0                                      # push / brake
F[3*top+1, 5] = 1.0                                    # side push
for sx in P.STRAP_X:                                   # battery hangs from straps
    n = near3(lambda x, y, z: (np.abs(x - sx) < 8) & (np.abs(y) < 16) & (np.abs(z - P.DECK_Z) < 3.6))
    F[3*n+2, 6] -= 1.0
n = near3(lambda x, y, z: (np.min([np.hypot(x - hx, y - hy) for hx, hy in P.DRV_HOLES], 0) < 6) & (z > ZTOP - 0.5))
F[3*n+2, 7] -= 1.0                                     # driver board (+ Pi)
n = near3(lambda x, y, z: (np.abs(y) > P.F - 4) & (np.abs(y) < P.F + P.PLATE_T + 4) & (np.hypot(x, z - P.AXLE_Z) > 13) & (np.hypot(x, z - P.MOTOR_AXIS_Z) < 26))
ddx, ddz = NX[n], NZ[n] - P.AXLE_Z
rn = np.hypot(ddx, ddz) + 1e-9
F[3*n, 8] = -ddz / rn; F[3*n+2, 8] = ddx / rn          # motor torque reaction
m_in = P.F - P.MOTOR_TOTAL_L
n = near3(lambda x, y, z: (np.abs(y) > m_in) & (np.abs(y) < m_in + 18) & (np.abs(np.hypot(x, z - P.MOTOR_AXIS_Z) - P.MOTOR_CLEAR_R - 2) < 5))
F[3*n+2, 9] -= 1.0                                     # motor weight on the rear cradles
F[3*top+2, 10] = dx / 40.0                             # pitch (reach backward)
rm_ = np.hypot(NX - P.MAST_X, NY)                      # camera mast whip: couple on the socket
n1 = near3(lambda x, y, z: (np.hypot(x - P.MAST_X, y) < 12) & (z > ZTOP - 4))
n2 = near3(lambda x, y, z: (np.hypot(x - P.MAST_X, y) < 12) & (np.abs(z - (ZTOP - P.MAST_SOCKET + 4)) < 4))
F[3*n1+1, 11] = 1.0 / max(len(n1), 1); F[3*n2+1, 11] = -1.0 / max(len(n2), 1)
F[3*n1, 11] = 0.5 / max(len(n1), 1); F[3*n2, 11] = -0.5 / max(len(n2), 1)
# caster plate pushes up equally on every bolt boss (a steel plate shares its load), wheels held:
# every boss needs its own load path, the old solid pad is not needed
for k_, sd in enumerate(sd_):
    F[sd, 12] += 1.0 / len(sd_) / max(len(sd), 1)
F /= np.abs(F).sum(0, keepdims=True)
W = np.array([1.0, 1.4, 1.4, 0.8, 0.5, 0.5, 0.8, 0.6, 0.8, 0.3, 0.8, 0.5, 1.0])
CASE_SET = ["caster"] * 12 + ["clamp"]
clamp_d = np.concatenate([3 * sup + d for d in range(3)])
SPRINGS = {"caster": (np.concatenate(sd_), np.concatenate(sv_)),
           "clamp": (clamp_d, np.full(len(clamp_d), 20.0))}

# density filter on active elements
rmin = 1.6
R = int(np.ceil(rmin))
rows, cols, vals = [], [], []
for a in range(-R + 1, R):
    for b in range(-R + 1, R):
        for c in range(-R + 1, R):
            w = rmin - np.sqrt(a*a + b*b + c*c)
            if w <= 0:
                continue
            i2, j2, k2 = ei + a, ej + b, ek + c
            ok = (i2 >= 0) & (i2 < nelx) & (j2 >= 0) & (j2 < nely) & (k2 >= 0) & (k2 < nelz)
            tgt = np.full(len(ei), -1)
            tgt[ok] = eidx_full[i2[ok], j2[ok], k2[ok]]
            m = tgt >= 0
            rows.append(np.nonzero(m)[0]); cols.append(tgt[m]); vals.append(np.full(m.sum(), w))
Hf = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                shape=(len(ei),) * 2).tocsr()
Hs = np.asarray(Hf.sum(1)).ravel()

is_solid = solid[ei, ej, ek]
is_free = free[ei, ej, ek]
mirror = eidx_full[ei, nely - 1 - ej, ek]             # y-mirror partner
volfrac = float(sys.argv[1]) if len(sys.argv) > 1 else 0.12
target = volfrac * is_free.sum()

x = np.where(is_free, volfrac, 1.0)
if os.environ.get("INIT"):                       # warm start from another design
    r0 = np.load(os.environ["INIT"])["rho"][ei, ej, ek]
    r0 = np.clip(r0 * volfrac / max(r0[is_free].mean(), 1e-6), 0.001, 1.0)
    x[is_free] = r0[is_free]
    print("warm start from", os.environ["INIT"], "mean", x[is_free].mean(), flush=True)
xPhys = x.copy()
NIT = int(sys.argv[2]) if len(sys.argv) > 2 else 90
OUT = sys.argv[3] if len(sys.argv) > 3 else "topo3d.npz"
AM = len(sys.argv) > 4 and sys.argv[4] == "am"

# ---------- overhang (AM) filter, Langelaar 2017, build direction +z ----------
# An element may only be as dense as its supporting elements one layer down
# (itself + 4 side neighbours = 45 deg). Cells kept clear inside the body
# (bay, motors, wheel wells) count as support material; outside the body
# counts as air. The first layer sits on the bed.
inplan_full = L.plan_sd(X, Y) < -0.5
aug = np.where(~active & inplan_full, 1.0, 0.0)
solid_full = solid
P_AM, EPS = 40.0, 1e-4
Q_AM = P_AM + np.log(5) / np.log(0.5)

def am_forward(xt):
    xp = xt.copy(); cache = []
    for k in range(1, nelz):
        below = np.where(active[:, :, k-1], xp[:, :, k-1], aug[:, :, k-1])
        bp = np.pad(below, 1)
        Ss = [bp[1:-1, 1:-1], bp[2:, 1:-1], bp[:-2, 1:-1], bp[1:-1, 2:], bp[1:-1, :-2]]
        sumP = sum(np.maximum(S_, 1e-3) ** P_AM for S_ in Ss)
        smax = sumP ** (1 / Q_AM)
        dsmax = [(1 / Q_AM) * sumP ** (1 / Q_AM - 1) * P_AM * np.maximum(S_, 1e-3) ** (P_AM - 1) for S_ in Ss]
        a, b_ = xt[:, :, k], smax
        r = np.sqrt((a - b_) ** 2 + EPS)
        val = 0.5 * (a + b_ - r + np.sqrt(EPS))
        da, db = 0.5 * (1 - (a - b_) / r), 0.5 * (1 + (a - b_) / r)
        sol = solid_full[:, :, k]
        xp[:, :, k] = np.where(sol, 1.0, val)
        cache.append((np.where(sol, 0.0, da), np.where(sol, 0.0, db), dsmax))
    return xp, cache

def am_backward(g, cache):
    lam = g.copy(); out = np.zeros_like(g)
    for k in range(nelz - 1, 0, -1):
        da, db, dsmax = cache[k - 1]
        out[:, :, k] = lam[:, :, k] * da
        t = lam[:, :, k] * db
        acc = np.zeros((nelx + 2, nely + 2))
        acc[1:-1, 1:-1] += t * dsmax[0]; acc[2:, 1:-1] += t * dsmax[1]; acc[:-2, 1:-1] += t * dsmax[2]
        acc[1:-1, 2:] += t * dsmax[3]; acc[1:-1, :-2] += t * dsmax[4]
        lam[:, :, k - 1] += acc[1:-1, 1:-1] * active[:, :, k - 1]
    out[:, :, 0] = lam[:, :, 0]
    return out

am_cache = None
Emin, penal = 1e-6, 3.0
U = np.zeros((ndof, F.shape[1])); Uprev = U.copy()
# rigid body modes for AMG
cx_, cy_, cz_ = NX - NX.mean(), NY - NY.mean(), NZ - NZ.mean()
rbm_full = np.zeros((ndof, 6))
rbm_full[0::3, 0] = 1; rbm_full[1::3, 1] = 1; rbm_full[2::3, 2] = 1
rbm_full[0::3, 3] = -cy_; rbm_full[1::3, 3] = cx_
rbm_full[1::3, 4] = -cz_; rbm_full[2::3, 4] = cy_
rbm_full[0::3, 5] = cz_; rbm_full[2::3, 5] = -cx_
rbm = rbm_full[freedofs] / 100.0
for it in range(NIT):
    t0 = time.time()
    Ev = Emin + xPhys ** penal * (1 - Emin)
    sK = (KE.ravel()[None, :] * Ev[:, None]).ravel()
    K0 = coo_matrix((sK, (iK, jK)), shape=(ndof, ndof))
    for sname, (sd, sv) in SPRINGS.items():
        K = (K0 + coo_matrix((sv, (sd, sd)), shape=(ndof, ndof))).tocsc()
        Kf = K[freedofs][:, freedofs].tocsr()
        Kf = 0.5 * (Kf + Kf.T)
        ml = pyamg.smoothed_aggregation_solver(Kf, B=rbm, max_coarse=800,
                                               strength=("symmetric", {"theta": 0.0}))
        M = ml.aspreconditioner(cycle="V")
        for lc in [i for i, c in enumerate(CASE_SET) if c == sname]:
            u, info = cg(Kf, F[freedofs, lc], x0=Uprev[freedofs, lc], rtol=1e-5, maxiter=300, M=M)
            U[freedofs, lc] = u
    Uprev = U.copy()
    dc = np.zeros(len(ei)); c = 0.0
    for lc in range(F.shape[1]):
        Ue = U[edof, lc]
        ce = np.einsum("ij,jk,ik->i", Ue, KE, Ue)
        c += W[lc] * (Ev * ce).sum()
        dc -= W[lc] * penal * xPhys ** (penal - 1) * (1 - Emin) * ce
    if AM and am_cache is not None:
        g = np.zeros(active.shape); g[ei, ej, ek] = dc
        dc = am_backward(g, am_cache)[ei, ej, ek]
    dc = np.asarray(Hf @ (dc / Hs)).ravel()
    dc = 0.5 * (dc + dc[mirror])
    move = 0.2 if it < 40 else 0.08
    l1, l2 = 0.0, 1e12
    xf, dcf = x[is_free], np.minimum(dc[is_free], -1e-30)
    while (l2 - l1) / (l1 + l2 + 1e-30) > 1e-4:
        lm = 0.5 * (l1 + l2)
        xn = np.clip(xf * np.sqrt(-dcf / lm), np.maximum(0, xf - move), np.minimum(1, xf + move))
        l1, l2 = (lm, l2) if xn.sum() > target else (l1, lm)
    change = np.abs(xn - xf).max()
    x[is_free] = xn
    xPhys = np.asarray(Hf @ x).ravel() / Hs
    if it >= NIT // 2:
        beta = 1 + 7 * (it - NIT // 2) / (NIT - NIT // 2)
        xPhys = (np.tanh(beta * 0.5) + np.tanh(beta * (xPhys - 0.5))) / (2 * np.tanh(beta * 0.5))
    xPhys[is_solid] = 1.0
    if AM:
        xt = np.zeros(active.shape); xt[ei, ej, ek] = xPhys
        xp, am_cache = am_forward(xt)
        xPhys = xp[ei, ej, ek]
        xPhys[is_solid] = 1.0
    print(f"it {it:3d} c={c:.4e} ch={change:.3f} vol={xPhys[is_free].mean():.3f} {time.time()-t0:.1f}s", flush=True)

rho = np.zeros(active.shape)
rho[ei, ej, ek] = xPhys
np.savez(OUT, rho=rho, c=c, solid=solid, active=active, x0=X0, y0=Y0, z0=ZBOT, h=H)
print("saved", OUT, "c", c, "vol", xPhys[is_free].mean())

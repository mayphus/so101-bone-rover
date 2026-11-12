"""3D SIMP topology optimization of the chassis (hex8, multi-load, y-mirrored).

Design space: the stadium plan (330 x 170) from z=30 to z=ZTOP, minus
keep-outs (battery bay, motor pods + wheels, casters). Supports: motor pod
bolt pads (wheels, fully fixed) and caster pads (vertical only).
Loads: SO-101 arm on the top pad (weight, pitch / roll moment at full reach,
yaw twist, side push), battery hanging from the bay ceiling, Pi on top.
Writes topo3d.npz used by build3d.py.
"""
import sys, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu, cg
import pyamg
import params as P

H = 6.0
ZBOT, ZTOP = 30.0, 96.0
X0, Y0 = -168.0, -87.0
nelx, nely, nelz = 56, 29, int((ZTOP - ZBOT) / H)
xs = X0 + H * (np.arange(nelx) + 0.5)
ys = Y0 + H * (np.arange(nely) + 0.5)
zs = ZBOT + H * (np.arange(nelz) + 0.5)
X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
Z0 = P.DECK_Z0

A_ = P.DECK_L / 2 - P.DECK_R
def stadium_sd(x, y):
    return np.hypot(np.maximum(np.abs(x) - A_, 0), y) - P.DECK_R

inplan = stadium_sd(X, Y) < -1.0
keep = np.zeros_like(inplan)
keep |= (np.abs(X) < 64) & (np.abs(Y) < 44) & (Z < Z0)               # battery bay
keep |= (X > -72) & (X < 50) & (np.abs(Y) > 40) & (Z < Z0)           # pods
keep |= (np.abs(Y) > 78) & (np.hypot(X, Z - P.WHEEL_R) < 60)         # wheels
for cx in (P.CASTER_X, -P.CASTER_X):
    keep |= (np.hypot(X - cx, Y) < 26) & (Z < Z0)                    # casters
active = inplan & ~keep

solid = np.zeros_like(active)
solid |= (np.hypot(X - P.HUB_X, Y) < P.HUB_R) & (Z > ZTOP - H)        # arm pad
pads = []
for s in (1, -1):
    pads += P.pod_bolts(s)
for b in pads:
    solid |= (np.hypot(X - b[0], Y - b[1]) < 9) & (Z > Z0) & (Z < Z0 + H)
for cx in (P.CASTER_X, -P.CASTER_X):
    solid |= (np.hypot(X - cx, Y) < 22) & (Z > Z0) & (Z < Z0 + H)
for b in P.PI_HOLES:
    solid |= (np.hypot(X - b[0], Y - b[1]) < 7) & (Z > ZTOP - H)
for sx in (1, -1):
    solid |= (np.abs(X - sx * P.STRAP_X) < 7) & (np.abs(Y) < 18) & (Z > Z0) & (Z < Z0 + H)
for sx in (1, -1):                                                    # seam flanges
    band = np.abs(X - sx * P.SEAM_X) < 6
    halfw = np.sqrt(np.maximum(P.DECK_R**2 - np.maximum(np.abs(X) - A_, 0)**2, 0))
    ring = (np.abs(Y) > halfw - 13) | (Z > ZTOP - H) | (Z < Z0 + H)
    solid |= band & ring & (Z > Z0)
solid &= active
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

zpad = ZBOT + H * np.round((Z0 - ZBOT) / H)        # node plane at bottom of pads
fixed = []
for b in pads:
    n = near(b[0], b[1], zpad, 7)
    fixed += list(3 * n) + list(3 * n + 1) + list(3 * n + 2)
for cx in (P.CASTER_X, -P.CASTER_X):
    n = near(cx, 0, zpad, 16)
    fixed += list(3 * n + 2)
fixed = np.unique(fixed)
freedofs = np.setdiff1d(np.arange(ndof), fixed)

top = near(P.HUB_X, 0, ZTOP, 40)
F = np.zeros((ndof, 8))
dx, dy = NX[top] - P.HUB_X, NY[top]
F[3*top+2, 0] = -1.0                                   # arm weight
F[3*top+2, 1] = -dx / 40.0                             # pitch (reach forward)
F[3*top+2, 2] = -dy / 40.0 - 0.3                       # roll (reach sideways)
rr = np.hypot(dx, dy) + 1e-9
F[3*top, 3] = -dy / rr; F[3*top+1, 3] = dx / rr        # yaw twist
F[3*top, 4] = 1.0                                      # push / brake
F[3*top+1, 5] = 1.0                                    # side push
for sx in (1, -1):                                     # battery hangs from straps
    n = near(sx * P.STRAP_X, 0, zpad, 14, rz=0.1)
    F[3*n+2, 6] -= 1.0
for b in P.PI_HOLES:
    n = near(b[0], b[1], ZTOP, 5)
    F[3*n+2, 7] -= 1.0
F /= np.abs(F).sum(0, keepdims=True)
W = np.array([1.0, 1.4, 1.4, 0.8, 0.5, 0.5, 0.8, 0.25])

# density filter on active elements
rmin = 1.7
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
xPhys = x.copy()
NIT = int(sys.argv[2]) if len(sys.argv) > 2 else 90
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
    K = coo_matrix((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()
    Kf = K[freedofs][:, freedofs].tocsr()
    Kf = 0.5 * (Kf + Kf.T)
    if True:
        ml = pyamg.smoothed_aggregation_solver(Kf, B=rbm, max_coarse=800,
                                               strength=("symmetric", {"theta": 0.0}))
        M = ml.aspreconditioner(cycle="V")
    for lc in range(F.shape[1]):
        u, info = cg(Kf, F[freedofs, lc], x0=Uprev[freedofs, lc], rtol=1e-5, maxiter=300, M=M)
        U[freedofs, lc] = u
    Uprev = U.copy()
    dc = np.zeros(len(ei)); c = 0.0
    for lc in range(F.shape[1]):
        Ue = U[edof, lc]
        ce = np.einsum("ij,jk,ik->i", Ue, KE, Ue)
        c += W[lc] * (Ev * ce).sum()
        dc -= W[lc] * penal * xPhys ** (penal - 1) * (1 - Emin) * ce
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
    print(f"it {it:3d} c={c:.4e} ch={change:.3f} vol={xPhys[is_free].mean():.3f} {time.time()-t0:.1f}s", flush=True)

rho = np.zeros(active.shape)
rho[ei, ej, ek] = xPhys
np.savez("topo3d.npz", rho=rho, solid=solid, active=active, x0=X0, y0=Y0, z0=ZBOT, h=H)
print("saved topo3d.npz")

"""2D SIMP topology optimization of the deck truss layout.

Multi-load-case minimum compliance (top88 style, density filter, OC update),
mirrored about y=0. Loads come from the arm hub (push/pull, side, twist) and
the battery straps; supports are the motor-pod and caster bolts.
Writes topo_density.npz used by build.py.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.ndimage import distance_transform_edt
import params as P

H = 2.0                                   # element size (mm)
nelx, nely = int(P.DECK_L / H) + 2, int(P.DECK_W / H) + 2
x0, y0 = -nelx * H / 2, -nely * H / 2
xc = x0 + H * (np.arange(nelx) + 0.5)
yc = y0 + H * (np.arange(nely) + 0.5)
X, Y = np.meshgrid(xc, yc, indexing="xy")   # shape (nely, nelx)

def stadium_sd(x, y):
    """Signed distance (negative inside) to the stadium outline."""
    a = P.DECK_L / 2 - P.DECK_R
    dx = np.maximum(np.abs(x) - a, 0)
    return np.hypot(dx, y) - P.DECK_R

inside = stadium_sd(X, Y) < 0
solid = np.zeros_like(inside)
void = ~inside

def disc(cx, cy, r):
    return np.hypot(X - cx, Y - cy) < r

solid |= inside & (stadium_sd(X, Y) > -4.0)            # perimeter rim
solid |= disc(P.HUB_X, 0, P.HUB_R) & ~disc(P.HUB_X, 0, 15)   # arm hub ring
void |= disc(P.HUB_X, 0, 13)                                  # cable hole
for s in (1, -1):
    for b in P.pod_bolts(s):
        solid |= disc(*b, 7)
for xc_ in (P.CASTER_X, -P.CASTER_X):
    solid |= disc(xc_, 0, P.CASTER_BOLT_R + 7)
for xs in (P.SEAM_X, -P.SEAM_X):
    solid |= inside & (np.abs(X - xs) < 4.0)            # joint face strip
    for b in P.splice_bolts(xs):
        solid |= disc(*b, 6.5)
for b in P.PI_HOLES:
    solid |= disc(*b, 5)
for sx in (1, -1):                                      # strap slot frames
    solid |= (np.abs(X - sx * P.STRAP_X) < 6) & (np.abs(Y) < 16)
for s in (1, -1):                                       # cable pass-through rings
    solid |= disc(-45, s * 38, 11)
free = inside & ~solid

# ---------- FE setup (bilinear quad, plane stress) ----------
E0, Emin, nu, penal = 1.0, 1e-9, 0.3, 3.0
k = np.array([1/2-nu/6, 1/8+nu/8, -1/4-nu/12, -1/8+3*nu/8,
              -1/4+nu/12, -1/8-nu/8, nu/6, 1/8-3*nu/8])
KE = E0/(1-nu**2)*np.array([
    [k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]],
    [k[1],k[0],k[7],k[6],k[5],k[4],k[3],k[2]],
    [k[2],k[7],k[0],k[5],k[6],k[3],k[4],k[1]],
    [k[3],k[6],k[5],k[0],k[7],k[2],k[1],k[4]],
    [k[4],k[5],k[6],k[7],k[0],k[1],k[2],k[3]],
    [k[5],k[4],k[3],k[2],k[1],k[0],k[7],k[6]],
    [k[6],k[3],k[4],k[1],k[2],k[7],k[0],k[5]],
    [k[7],k[2],k[1],k[4],k[3],k[6],k[5],k[0]]])

# node (i=col 0..nelx, j=row 0..nely) -> id = i*(nely+1)+j
def nid(i, j):
    return i * (nely + 1) + j

ndof = 2 * (nelx + 1) * (nely + 1)
ii, jj = np.meshgrid(np.arange(nelx), np.arange(nely), indexing="xy")
n1 = nid(ii, jj + 1).ravel()   # lower-left in top88 orientation
n2 = nid(ii + 1, jj + 1).ravel()
edof = np.stack([2*n1, 2*n1+1, 2*n2, 2*n2+1,
                 2*n2-2+0, 2*n2-1, 2*n1-2, 2*n1-1], 1)
# simpler explicit: nodes (i,j+1),(i+1,j+1),(i+1,j),(i,j)
na, nb, nc, nd = nid(ii, jj+1).ravel(), nid(ii+1, jj+1).ravel(), nid(ii+1, jj).ravel(), nid(ii, jj).ravel()
edof = np.stack([2*na, 2*na+1, 2*nb, 2*nb+1, 2*nc, 2*nc+1, 2*nd, 2*nd+1], 1)
iK = np.kron(edof, np.ones((8, 1), int)).ravel()
jK = np.kron(edof, np.ones((1, 8), int)).ravel()

nx_ = x0 + H * np.arange(nelx + 1)
ny_ = y0 + H * np.arange(nely + 1)
NX, NY = np.meshgrid(nx_, ny_, indexing="ij")     # (nelx+1, nely+1) matches nid
NXf, NYf = NX.ravel(), NY.ravel()

def nodes_near(cx, cy, r):
    return np.where(np.hypot(NXf - cx, NYf - cy) < r)[0]

fixed = []
for s in (1, -1):
    for b in P.pod_bolts(s):
        fixed += list(nodes_near(*b, 4.0))
fixed = np.unique(fixed)
fixed_dofs = np.concatenate([2*fixed, 2*fixed+1])
free_dofs = np.setdiff1d(np.arange(ndof), fixed_dofs)

# load cases
ring = np.where(np.abs(np.hypot(NXf - P.HUB_X, NYf) - P.HUB_R + 2) < 2.2)[0]
F = np.zeros((ndof, 20))
F[2*ring, 0] = 1.0                               # arm pushes forward
F[2*ring+1, 1] = 1.0                             # arm pushes sideways
dx, dy = NXf[ring] - P.HUB_X, NYf[ring]
rr = np.hypot(dx, dy)
F[2*ring, 2] = -dy / rr; F[2*ring+1, 2] = dx / rr   # arm yaw torque
F[2*ring, 3] = -1.0                              # braking
for sx in (1, -1):                               # battery straps, lateral
    nn = nodes_near(sx * P.STRAP_X, 0, 8)
    F[2*nn+1, 4] += 0.4
for b in P.PI_HOLES:                             # electronics, small
    nn = nodes_near(*b, 3)
    F[2*nn, 5] += 0.15
# bumper impacts at 8 rim points, pointing inward
a_ = P.DECK_L / 2 - P.DECK_R
for n, ang in enumerate(np.linspace(0, 2*np.pi, 8, endpoint=False)):
    px = np.cos(ang) * (P.DECK_R - 3) + np.sign(np.cos(ang)) * a_ * (abs(np.cos(ang)) > 0.1)
    py = np.sin(ang) * (P.DECK_R - 3)
    nn = nodes_near(px, py, 5)
    F[2*nn, 6+n] = -np.cos(ang); F[2*nn+1, 6+n] = -np.sin(ang)
col = 14
for xc_ in (P.CASTER_X, -P.CASTER_X):            # caster reactions
    nn = np.concatenate([nodes_near(*b, 3.5) for b in P.caster_bolts(xc_)])
    F[2*nn, col] = 1.0; F[2*nn+1, col+1] = 1.0; col += 2
nn = np.concatenate([nodes_near(*b, 3) for b in P.PI_HOLES])
F[2*nn, col] = 1.0; F[2*nn+1, col+1] = 1.0
F /= np.abs(F).sum(0, keepdims=True)
weights = np.array([1.0, 1.0, 1.0, 0.6, 0.3, 0.2] + [0.12]*8 + [0.8]*4 + [0.3]*2)

# density filter
rmin = 2.2
dist_idx = []
rows, cols, vals = [], [], []
R = int(np.ceil(rmin))
for e_j in range(nely):
    pass
ex = ii.ravel(); ey = jj.ravel()
eidx = np.arange(nelx * nely).reshape(nely, nelx)
for di in range(-R + 1, R):
    for dj in range(-R + 1, R):
        w = rmin - np.hypot(di, dj)
        if w <= 0:
            continue
        i2, j2 = ex + di, ey + dj
        ok = (i2 >= 0) & (i2 < nelx) & (j2 >= 0) & (j2 < nely)
        rows.append(eidx[ey[ok], ex[ok]])
        cols.append(eidx[j2[ok], i2[ok]])
        vals.append(np.full(ok.sum(), w))
Hf = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                shape=(nelx*nely,)*2).tocsr()
Hs = np.asarray(Hf.sum(1)).ravel()

solid_f, void_f, free_f = solid.ravel(), void.ravel(), free.ravel()
volfrac = 0.15
target = volfrac * free_f.sum()

def mirror_y(a):
    a2 = a.reshape(nely, nelx)
    return (0.5 * (a2 + a2[::-1, :])).ravel()

x = np.where(free_f, volfrac, 0.0)
x[solid_f] = 1.0
xPhys = x.copy()
for it in range(260):
    sK = ((KE.ravel()[None, :]) * (Emin + xPhys[:, None] ** penal * (E0 - Emin))).ravel()
    K = coo_matrix((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()
    Kf = K[free_dofs][:, free_dofs]
    lu = splu(Kf.tocsc())
    U = np.zeros((ndof, F.shape[1]))
    U[free_dofs] = lu.solve(F[free_dofs])
    dc = np.zeros(nelx * nely); c = 0.0
    for lc in range(F.shape[1]):
        Ue = U[edof, lc]
        ce = np.einsum("ij,jk,ik->i", Ue, KE, Ue)
        c += weights[lc] * ((Emin + xPhys**penal*(E0-Emin)) * ce).sum()
        dc += weights[lc] * (-penal * xPhys**(penal-1) * (E0-Emin) * ce)
    dc = np.asarray(Hf @ (dc / Hs)).ravel()
    dc = mirror_y(dc)
    # OC on free elements only
    l1, l2, move = 0.0, 1e9, (0.15 if it < 80 else 0.05)
    xf = x[free_f]; dcf = np.minimum(dc[free_f], -1e-30)
    while (l2 - l1) / (l1 + l2 + 1e-30) > 1e-4:
        lm = 0.5 * (l1 + l2)
        xn = np.clip(xf * np.sqrt(-dcf / lm), np.maximum(0, xf - move), np.minimum(1, xf + move))
        if xn.sum() > target:
            l1 = lm
        else:
            l2 = lm
    change = np.abs(xn - xf).max()
    x[free_f] = xn
    xPhys = np.asarray(Hf @ x).ravel() / Hs
    xPhys[solid_f] = 1.0; xPhys[void_f] = 0.0
    # mild projection late in the run for crisp members
    if it > 100:
        beta = 1 + (it - 100) / 20
        xPhys[free_f] = (np.tanh(beta*0.5) + np.tanh(beta*(xPhys[free_f]-0.5))) / (2*np.tanh(beta*0.5))
    if it % 10 == 0 or it == 259:
        print(f"it {it:3d}  compliance {c:.4e}  change {change:.3f}", flush=True)

np.savez("topo_density.npz", rho=xPhys.reshape(nely, nelx), x0=x0, y0=y0, h=H)
print("saved topo_density.npz", xPhys.reshape(nely, nelx).shape)

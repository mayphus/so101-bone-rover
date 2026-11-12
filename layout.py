"""Where the structure may not go (keep-outs) and where it must be solid
(mounting features). Shared by the optimizer (coarse grid) and the final
smooth build (fine grid). All functions take same-shape arrays X, Y, Z."""
import numpy as np
import shapely
from shapely.geometry import Polygon, box
import params as P

ARM_PAD = Polygon(P.ARM_OUTLINE).buffer(3.0, 16)
DRV_PAD = box(P.DRV_X - 34, -30, P.DRV_X + 34, 30)

def in_pad(X, Y):
    return shapely.contains_xy(ARM_PAD, X, Y)

def _rrect_sd(X, Y, x0, x1, hw, r):
    cx, hx = (x0 + x1) / 2, (x1 - x0) / 2
    qx = np.abs(X - cx) - (hx - r)
    qy = np.abs(Y) - (hw - r)
    return np.hypot(np.maximum(qx, 0), np.maximum(qy, 0)) + np.minimum(np.maximum(qx, qy), 0) - r

def plan_sd(X, Y):
    """Signed distance to the body plan (negative inside): rounded body + motor pods."""
    body = _rrect_sd(X, Y, P.X_TAIL, P.X_NOSE, P.HALF_W, P.CORNER_R)
    pods = _rrect_sd(X, Y, -P.POD_X, P.POD_X, P.F + P.PLATE_T, 8.0)
    return np.minimum(body, pods)

def r_motor(X, Z):
    return np.hypot(X, Z - P.MOTOR_AXIS_Z)

def r_axle(X, Z):
    return np.hypot(X, Z - P.AXLE_Z)

def keepout(X, Y, Z):
    aY = np.abs(Y)
    k = plan_sd(X, Y) > -0.5
    k |= Z > P.Z_TOP
    k |= Z < P.Z_FLOOR
    rc = np.hypot(X - P.CASTER_X, Y)
    k |= (rc < P.CASTER_SWING) & (Z < P.CASTER_H)                         # swivel caster
    m_in = P.F - P.MOTOR_TOTAL_L
    k |= (aY < P.F) & (aY > m_in - 0.5) & (r_motor(X, Z) < P.MOTOR_CLEAR_R)   # motor body
    k |= (aY <= m_in) & (np.hypot(X - P.PLUG_X, Z - P.MOTOR_AXIS_Z - P.PLUG_Z) < 13)  # plugs + cable
    k |= (aY > P.F + P.PLATE_T + 0.4) & (r_axle(X, Z) < P.WHEEL_R + 6)   # wheels
    k |= (aY > P.F - 0.1) & (r_axle(X, Z) < 7.5)                         # shaft boss
    b = P.BAY
    k |= (X > b["x0"]) & (X < b["x1"]) & (aY < b["y"]) & (Z > b["z0"]) & (Z < b["z1"])
    k |= (X > -4) & (X < b["x0"] + 1) & (aY < 9) & (Z > 52) & (Z < 76)   # motor cables -> bay
    k |= (np.hypot(X - P.ARM_CABLE[0], Y - P.ARM_CABLE[1]) < 10) & (Z > P.Z_TOP - 14)
    k |= (np.hypot(X - P.MAST_X, Y) < P.MAST_D / 2 + 0.3) & (Z > P.Z_TOP - P.MAST_SOCKET + 4)
    return k

def solid(X, Y, Z):
    aY = np.abs(Y)
    s = in_pad(X, Y) & (Z > P.Z_TOP - 7)                                 # arm pad = SO-101 footprint
    for hx, hy in P.ARM_HOLES:                                           # M5 insert bosses
        s |= (np.hypot(X - hx, Y - hy) < 7) & (Z > P.Z_TOP - 16)
    for hx, hy in P.DRV_HOLES:                                           # driver board bosses
        s |= (np.hypot(X - hx, Y - hy) < 5.5) & (Z > P.Z_TOP - 10)
    s |= (aY >= P.F) & (aY <= P.F + P.PLATE_T) & (r_motor(X, Z) < 24)    # face plates
    rm = r_motor(X, Z)
    m_in = P.F - P.MOTOR_TOTAL_L
    s |= (aY > m_in + 4) & (aY < m_in + 14) & (rm > P.MOTOR_CLEAR_R) & (rm < P.MOTOR_CLEAR_R + 5)  # rear cradle
    for hx, hy in P.CASTER_HOLES:                                        # caster bolt bosses only;
        s |= (np.hypot(X - P.CASTER_X - hx, Y - hy) < P.CASTER_BOSS) & (Z >= P.CASTER_H) & (Z < P.CASTER_H + 10)  # optimizer fills between
    for sx in P.STRAP_X:                                                 # battery strap bars
        s |= (np.abs(X - sx) < P.STRAP_HALF) & (aY < 18) & (Z > P.DECK_Z) & (Z < P.DECK_Z + 6)
    s |= (np.hypot(X - P.MAST_X, Y) < P.MAST_D / 2 + 6) & (Z > P.Z_TOP - P.MAST_SOCKET)   # mast socket
    return s & ~keepout(X, Y, Z)

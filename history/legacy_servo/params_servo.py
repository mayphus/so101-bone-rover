"""Shared dimensions for the SO-101 2-motor mobile base (all mm).

Frame: x = forward, y = left, z = up. Origin = wheel axle center line at
ground level (x=0, y=0, z=0). Axle height = WHEEL_R.
"""

# ---------- Printer ----------
BED = 180.0                     # Bambu Lab A1 mini: 180 x 180 x 180

# ---------- Drive servo: Feetech STS3215 (same as SO-101 arm) ----------
# VERIFY with calipers before printing the motor pods.
SERVO_L = 45.2                  # long side (lies along x)
SERVO_W = 24.7                  # short side (lies along z)
SERVO_H = 35.0                  # body depth along output shaft (lies along y)
SERVO_SHAFT_FROM_END = 12.0     # shaft center to nearest short end
SERVO_CLEAR = 0.25              # pocket clearance per side
HORN_PROTRUDE = 5.0             # horn face distance beyond servo face
HORN_D = 22.5                   # horn disc diameter (hub recess)
HORN_DEPTH = 2.5                # recess depth in wheel hub
HORN_BOLT_CIRCLE = 14.0         # horn screw circle diameter (VERIFY)
HORN_BOLT_D = 2.2               # M2 clearance

# ---------- Wheels ----------
WHEEL_R = 50.0                  # with TPU tire -> 100 mm wheel
WHEEL_W = 24.0
RIM_R = 44.0                    # rim outer radius (tire seat)
TIRE_T = WHEEL_R - RIM_R        # 6 mm TPU tire

# ---------- Deck ----------
DECK_L = 330.0                  # overall length (stadium outline)
DECK_W = 170.0                  # overall width
DECK_R = DECK_W / 2             # end radius
SEAM_X = 80.0                   # center piece = |x| < SEAM_X
SEAM_BAND = 16.0                # solid band each side of a seam
DECK_T_MIN = 7.0                # thinnest truss height
DECK_T_HUB = 15.0               # height at the arm hub
SERVO_FACE_Y = 84.0             # servo output face (y), wheel outboard of deck

# deck bottom height above ground: servo top + pod top plate
POD_TOP_T = 3.0
DECK_Z0 = WHEEL_R + SERVO_W / 2 + SERVO_CLEAR + POD_TOP_T   # ~65.6

# ---------- Arm hub ----------
HUB_X = 15.0                    # arm hub center (slightly ahead of axle)
HUB_R = 46.0                    # flat solid pad radius
HUB_HOLE_RINGS = (22.0, 34.0)   # M3 polar hole rings (VERIFY vs your arm base)
HUB_HOLE_N = 12

# ---------- Casters ----------
CASTER_X = 140.0                # front and rear
BALL_D = 25.4                   # 1" steel ball (or 25 mm marble)
CASTER_LIFT = 0.5               # ball sits this much above floor (tire squish)
CASTER_BOLT_R = 15.0

# ---------- Electronics ----------
PI_X = -122.0                   # Raspberry Pi 5 on rear piece (58 x 49 holes)
PI_HOLES = [(PI_X - 24.5, -29.0), (PI_X - 24.5, 29.0),
            (PI_X + 24.5, -29.0), (PI_X + 24.5, 29.0)]
STRAP_X = 52.0                  # battery strap slots at x = +-STRAP_X

# ---------- Hardware ----------
M3 = 3.4                        # M3 clearance
M3_NUT_AF = 5.7                 # nut across flats (+clearance)
M3_TAP = 2.8                    # self-tap into plastic
M25 = 2.3                       # M2.5 tap for Pi standoffs

# motor pod geometry (shared by deck + pod)
def servo_x_span():
    """Servo body x range: shaft at x=0, long side toward the rear."""
    return (-(SERVO_L - SERVO_SHAFT_FROM_END), SERVO_SHAFT_FROM_END)

POD_BOLTS_X = (-48.0, 26.0)
POD_BOLTS_Y = (54.0, 79.0)

def pod_bolts(side=1):
    return [(x, side * y) for x in POD_BOLTS_X for y in POD_BOLTS_Y]

def caster_bolts(xc):
    import math
    return [(xc + CASTER_BOLT_R * math.cos(math.radians(a)),
             CASTER_BOLT_R * math.sin(math.radians(a))) for a in (0, 120, 240)]

SPLICE_Y = (-60.0, -22.0, 22.0, 60.0)
SPLICE_DX = 9.0

def splice_bolts(xs):
    return [(xs + s * SPLICE_DX, y) for s in (-1, 1) for y in SPLICE_Y]

def hub_holes():
    import math
    out = []
    for i, r in enumerate(HUB_HOLE_RINGS):
        for k in range(HUB_HOLE_N):
            a = 2 * math.pi * (k + 0.5 * i) / HUB_HOLE_N
            out.append((HUB_X + r * math.cos(a), r * math.sin(a)))
    return out

"""Bone Rover v4: SO-101 mobile base, 2 x JGB37-520 (encoder) rear drive + 1 swivel caster (all mm).

Frame: x forward, y left, z up. Origin on the floor under the drive axle.
The left drivetrain is listed; the right one is its mirror (y -> -y).
Values marked VERIFY depend on the exact part you buy.
"""
import math, json as _json

BED = 180.0                       # Bambu Lab A1 mini

# ---------- JGB37-520 gearmotor with Hall encoder (VERIFY against your motor) ----------
GEARBOX_D = 37.0
GEARBOX_L = 26.5                  # 19-29 depending on ratio; 1:90 ~ 26.5
CAN_D = 35.5
CAN_L = 31.0
ENC_D, ENC_L = 34.0, 12.0         # round encoder PCB + magnet disc behind the can
MOTOR_TOTAL_L = GEARBOX_L + CAN_L + ENC_L
PLUG_L = 22.0                     # 6-pin PH2.0 header + plug + cable bend, behind the encoder
PLUG_X, PLUG_Z = 0.0, 9.0         # header position on the PCB (from the motor axis)
SHAFT_OFFSET = 7.0                # output shaft is 7 mm off the body axis
SHAFT_D = 6.0
SHAFT_L = 15.5
BOSS_D, BOSS_H = 12.0, 6.0        # locating boss around the shaft
FACE_PCD = 31.0                   # 6 x M3 tapped holes on this circle
FACE_HOLE_ANGLE0 = 30.0           # angle of first hole (VERIFY)

# ---------- bought wheel on the 6 mm D-shaft (hex coupling) ----------
WHEEL_D = 100.0                   # 100 mm rubber wheel for 6 mm shaft (VERIFY)
WHEEL_R = WHEEL_D / 2
WHEEL_W = 40.0
WHEEL_GAP = 8.0                   # motor face to wheel inner face

AXLE_Z = WHEEL_R
MOTOR_AXIS_Z = AXLE_Z + SHAFT_OFFSET
PLATE_T = 4.0
F = MOTOR_TOTAL_L + PLUG_L + 1.5  # motor face plane |y|: the two plugs meet at y = 0
WHEEL_Y0 = F + PLATE_T + WHEEL_GAP - PLATE_T
WHEEL_Y1 = WHEEL_Y0 + WHEEL_W
MOTOR_CLEAR_R = GEARBOX_D / 2 + 1.5

# ---------- swivel caster (2" rubber, top plate; VERIFY) ----------
CASTER_X = 195.0
CASTER_H = 64.0                   # floor to top of plate
CASTER_SWING = 50.0               # keep-clear radius for the swivelling wheel
CASTER_PAD = 76.0                 # square flat pad under the frame
CASTER_HOLES = [(sx * 27.0, sy * 12.0) for sx in (1, -1) for sy in (1, -1)] + \
               [(sx * 19.0, sy * 19.0) for sx in (1, -1) for sy in (1, -1)]   # 54x24 and 38x38 patterns
M4_INSERT_D, M4_INSERT_L = 5.6, 8.0
CASTER_BOSS = 6.5                 # boss radius around each caster bolt (no solid plate)

# ---------- body envelope ----------
X_TAIL, X_NOSE = -80.0, 235.0
HALF_W = 80.0                     # body half width; motor pods stick out to the face plates
CORNER_R = 34.0
POD_X = 42.0                      # motor pod region |x| < POD_X reaches out to the plates (room for collars)
Z_FLOOR = 25.0                    # ground clearance
Z_TOP = 100.0
DECK_Z = 68.0                     # battery bay ceiling (straps)
BODY_L = X_NOSE - X_TAIL

# ---------- print split: 4 pieces, joined by bolted collars where members cross ----------
SPLIT_X = 104.0                   # rear body | front body: ahead of the arm pad, where few members cross
SPLIT_Y = 80.0                    # motor pods outside |y| > SPLIT_Y
COLLAR_GROW, COLLAR_T = 9.0, 5.0  # collar = member section grown 9 mm, 5 mm each side of the cut
PIECES = ["frame_rear", "frame_front", "frame_pod_left", "frame_pod_right"]

# ---------- arm ----------
HUB_X = 50.0                      # SO-101 pan axis, just ahead of the drive axle
_ARM = _json.load(open(__file__.rsplit("/", 1)[0] + "/so101_base.json" if "/" in __file__ else "so101_base.json"))
ARM_OUTLINE = [(HUB_X + x, y) for x, y in _ARM["outline"]]
ARM_HOLES = [(HUB_X + x, y) for x, y in _ARM["holes"]]
M5_INSERT_D, M5_INSERT_L = 6.4, 9.0
ARM_CABLE = (HUB_X - 40.0, 0.0)

# ---------- battery, electronics, camera ----------
BAY = dict(x0=30.0, x1=120.0, y=46.0, z0=Z_FLOOR - 4, z1=DECK_Z)   # 2S/3S pack up to ~80 x 80 x 40
STRAP_X = (50.0, 108.0)
STRAP_HALF = 9.0                  # strap bar half width: ~3 mm wall beside each strap slot
DRV_X = 160.0                     # Waveshare General Driver for Robots, 65 x 65, M3 at 58 x 49
DRV_HOLES = [(DRV_X + sx * 29.0, sy * 24.5) for sx in (1, -1) for sy in (1, -1)]
DRV_STANDOFF = 6.0                # Pi 5 (same 58 x 49 pattern) can stack above on 20 mm standoffs
MAST_X = -58.0                    # 10 mm carbon / aluminium tube socket
MAST_D = 10.0
MAST_L = 380.0                    # tube length (cut to taste)
MAST_SOCKET = 34.0
CAM_TILT = 45.0                   # degrees down: view centre lands ~420 mm ahead, the arm workspace

# ---------- hardware ----------
M3, M3_TAP, M3_INSERT_D, M3_INSERT_L, M25_TAP = 3.4, 2.8, 4.2, 6.0, 2.3

def face_holes():
    """(x, z) of the 6 motor face holes, around the gearbox axis."""
    return [(FACE_PCD / 2 * math.cos(math.radians(FACE_HOLE_ANGLE0 + 60 * k)),
             MOTOR_AXIS_Z + FACE_PCD / 2 * math.sin(math.radians(FACE_HOLE_ANGLE0 + 60 * k)))
            for k in range(6)]

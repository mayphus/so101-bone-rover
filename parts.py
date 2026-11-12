"""Printed camera head + reference models of bought parts for the viewer."""
import json, math, os
import numpy as np
import trimesh
from manifold3d import Manifold, OpType
import params as P

os.makedirs("stl/parts", exist_ok=True)
os.makedirs("stl/assembly", exist_ok=True)
for f_ in os.listdir("stl/assembly"):                      # drop v3 leftovers
    if f_.startswith("_ref_") or f_.startswith("caster_cap") or f_.startswith("camera_head"):
        os.remove(f"stl/assembly/{f_}")

def U(ms):
    return Manifold.batch_boolean([m for m in ms if not m.is_empty()], OpType.Add)

def tmesh(m):
    o = m.to_mesh()
    return trimesh.Trimesh(np.asarray(o.vert_properties)[:, :3], np.asarray(o.tri_verts), process=False)

def cyly(r, y0, y1, x, z, seg=64):
    return Manifold.cylinder(y1 - y0, r, r, seg).rotate((-90, 0, 0)).translate((x, y0, z))

def cylz(r, z0, z1, x, y, seg=48):
    return Manifold.cylinder(z1 - z0, r, r, seg).translate((x, y, z0))

def cube(x0, x1, y0, y1, z0, z1):
    return Manifold.cube((x1 - x0, y1 - y0, z1 - z0)).translate((x0, y0, z0))

def both_sides(m):
    return U([m, m.mirror((0, 1, 0))])

def save(m, name):
    tmesh(m).export(f"stl/assembly/{name}.stl")

# ---------------- JGB37-520 with Hall encoder (left side, mirrored) ----------------
F, ZA, ZM = P.F, P.AXLE_Z, P.MOTOR_AXIS_Z
g0 = F - P.GEARBOX_L
gear = cyly(P.GEARBOX_D / 2, g0, F, 0, ZM, 96)
for hx, hz in P.face_holes():
    gear -= cyly(1.3, F - 5, F + 0.1, hx, hz, 16)
boss = cyly(P.BOSS_D / 2, F, F + P.BOSS_H, 0, ZA, 48)
shaft = cyly(P.SHAFT_D / 2, F, F + P.BOSS_H + P.SHAFT_L, 0, ZA, 32) - cube(-5, 5, F + P.BOSS_H + 3, F + 40, ZA + 2.5, ZA + 5)
c0 = g0 - P.CAN_L
can = cyly(P.CAN_D / 2, c0, g0, 0, ZM, 96) - cube(-30, 30, c0 - 1, g0 + 1, ZM + P.CAN_D / 2 - 2.5, ZM + 30)  # flat on top
e0 = c0 - P.ENC_L
magnet = cyly(11, c0 - 4, c0, 0, ZM, 48)
pcb = cyly(P.ENC_D / 2, e0 + 2, e0 + 3.6, 0, ZM, 96)
hdr = cube(P.PLUG_X - 7.5, P.PLUG_X + 7.5, e0 - 4.5, e0 + 2, ZM + P.PLUG_Z - 3, ZM + P.PLUG_Z + 3)
plug = cube(P.PLUG_X - 8, P.PLUG_X + 8, e0 - 11, e0 - 4.5, ZM + P.PLUG_Z - 3.5, ZM + P.PLUG_Z + 3.5)
hall = U([cube(sx - 1.5, sx + 1.5, e0 + 3.6, e0 + 5, ZM - 13, ZM - 9) for sx in (-4, 4)])
save(both_sides(U([gear, boss, shaft])), "_ref_gearbox")
save(both_sides(U([can, magnet])), "_ref_motor")
save(both_sides(U([pcb, hall])), "_ref_encoder")
save(both_sides(U([hdr, plug])), "_ref_plug")

# ---------------- 100 mm wheel ----------------
y0, y1 = P.WHEEL_Y0, P.WHEEL_Y1
R = P.WHEEL_R
tire = cyly(R - 3, y0, y1, 0, ZA, 96) - cyly(R - 12, y0 - 1, y1 + 1, 0, ZA, 96)
lugs = []
for k in range(20):                                   # offset chevron lugs for grass and gravel
    a = 360 * k / 20
    off = 5 if k % 2 else -5
    lug = cube(R - 4, R, y0 + 3 + (off > 0) * 14, y0 + 23 + (off > 0) * 14, -4, 4)
    lugs.append(lug.rotate((0, a, 0)).translate((0, 0, ZA)))
tire = U([tire] + lugs)
hub = cyly(R - 12, y0 + 4, y1 - 4, 0, ZA, 64)
for k in range(5):
    a = math.radians(72 * k + 18)
    hub -= cyly(9, y0, y1, 24 * math.cos(a), ZA + 24 * math.sin(a), 32)
hub = U([hub, cyly(7.5, F + P.PLATE_T + 2.5, y0 + 6, 0, ZA, 6)])
save(both_sides(tire), "_ref_tire")
save(both_sides(hub), "_ref_wheelhub")

# ---------------- 2" swivel caster (trailing) ----------------
cx, ch = P.CASTER_X, P.CASTER_H
plate = cube(cx - 32.5, cx + 32.5, -24, 24, ch - 3, ch)
for hx, hy in P.CASTER_HOLES[:4]:
    plate -= cylz(2.2, ch - 4, ch + 1, cx + hx, hy, 16)
swivel = cylz(17, ch - 10, ch - 3, cx, 0, 48)
wx, wz = cx - 24, 25.4
fork = U([cube(wx - 6, cx + 12, sy * 13 - 1.5, sy * 13 + 1.5, wz - 4, ch - 10) for sy in (1, -1)]
         + [cube(cx - 14, cx + 14, -14.5, 14.5, ch - 13, ch - 10)])
cwheel = cyly(25.4, -10, 10, wx, wz, 64) - cyly(9, -11, 11, wx, wz, 32)
save(U([plate, swivel, fork]), "_ref_caster_fork")
save(cwheel, "_ref_caster_wheel")

# ---------------- battery, driver (+ optional Pi), mast, camera ----------------
b = P.BAY
save(cube(b["x0"] + 5, b["x1"] - 5, -b["y"] + 6, b["y"] - 6, b["z0"] + 8, b["z1"] - 2), "_ref_battery")
zs = P.Z_TOP + P.DRV_STANDOFF
drv = cube(P.DRV_X - 32.5, P.DRV_X + 32.5, -32.5, 32.5, zs, zs + 1.6)
drv = U([drv, cube(P.DRV_X - 32, P.DRV_X - 6, -12, 12, zs + 1.6, zs + 5.5),          # ESP32 module
         cube(P.DRV_X + 14, P.DRV_X + 30, -30, -10, zs + 1.6, zs + 8),               # connectors
         cube(P.DRV_X + 14, P.DRV_X + 30, 10, 30, zs + 1.6, zs + 8)])
for hx, hy in P.DRV_HOLES:
    drv = U([drv, cylz(2.5, P.Z_TOP, zs, hx, hy, 6)])
save(drv, "_ref_driver")
zp = zs + 22
pi = U([cube(P.DRV_X - 42.5, P.DRV_X + 42.5, -28, 28, zp, zp + 1.6),
        cube(P.DRV_X - 30, P.DRV_X + 12, -24, 24, zp + 1.6, zp + 16)]
       + [cylz(2.5, zs + 1.6, zp, hx, hy, 6) for hx, hy in P.DRV_HOLES])
save(pi, "_ref_pi5")

zb = P.Z_TOP - P.MAST_SOCKET + 4
zt = zb + P.MAST_L
save(cylz(P.MAST_D / 2, zb, zt, P.MAST_X, 0, 32), "_ref_mast")

# printed camera head: tube clamp + arm + camera plate tilted CAM_TILT down
clamp = cylz(P.MAST_D / 2 + 5, zt - 24, zt + 3, P.MAST_X, 0, 64)
clamp = U([clamp, cube(P.MAST_X - 22, P.MAST_X - 6, -6, 6, zt - 24, zt - 6)])         # pinch ears
clamp -= cylz(P.MAST_D / 2 + 0.2, zt - 30, zt, P.MAST_X, 0, 48)
clamp -= cube(P.MAST_X - 30, P.MAST_X - 3, -0.8, 0.8, zt - 30, zt - 4)                # slit
clamp -= Manifold.cylinder(30, P.M3 / 2, P.M3 / 2, 16).rotate((90, 0, 0)).translate((P.MAST_X - 15, 15, zt - 15))
t = math.radians(P.CAM_TILT)
pc = np.array([P.MAST_X + 30, 0, zt - 2])                # camera plate centre
plate = Manifold.cube((4, 40, 40), center=True)
plate -= Manifold.cylinder(10, 3.3, 3.3, 32).rotate((0, 90, 0)).translate((-5, 0, 0))    # 1/4"-20 screw
for a_, b_ in ((14, 14), (14, -14), (-14, 14), (-14, -14)):                            # 4 x M2, 28 x 28
    plate -= Manifold.cylinder(10, 1.1, 1.1, 12).rotate((0, 90, 0)).translate((-5, a_, b_))
plate = plate.rotate((0, P.CAM_TILT, 0)).translate(tuple(pc))
arm = U([cube(P.MAST_X + 4, pc[0] - 2, -8, 8, zt - 6, zt + 3)])
head = U([clamp, arm, plate])
save(head, "camera_head")
th = tmesh(head); lo, hi = th.bounds
th.apply_translation([-(lo[0] + hi[0]) / 2, -(lo[1] + hi[1]) / 2, -lo[2]]); th.export("stl/parts/camera_head.stl")
n = np.array([math.cos(t), 0, -math.sin(t)])          # camera looks forward and down
cam = Manifold.cube((8, 32, 32), center=True)
cam = U([cam, Manifold.cylinder(10, 7, 6, 32).rotate((0, 90, 0)).translate((4, 0, 0))])
cam = cam.rotate((0, P.CAM_TILT, 0))
save(cam.translate(tuple(pc + n * 8)), "_ref_camera")

wm = trimesh.load("vendor/so101_wristcam/SO-ARM101_camera_wrist_mount.stl")
json.dump([dict(name="so101_wrist_camera_mount (official)", qty=1, material="PETG (or PLA+)", size_mm=wm.extents.round(1).tolist(),
                fits_a1_mini=True, volume_cm3=round(wm.volume / 1000, 1), est_grams=round(wm.volume / 1000 * 1.27 * 0.8),
                note="TheRobotStudio SO101_Wrist_Cam_Hex-Nut_Mount_32x32_UVC_Module. Two M3 x 8 + hex nuts into the wrist roll part; camera on 4 x M2."),
           dict(name="camera_head", qty=1, material="PETG (or PLA+)", size_mm=(hi - lo).round(1).tolist(),
                fits_a1_mini=True, volume_cm3=round(head.volume() / 1000, 1),
                est_grams=round(head.volume() / 1000 * 1.27 * 0.8),
                note="Clamps on the 10 mm mast with one M3 pinch screw. Camera on a 1/4\"-20 screw or 4 x M2 (32 x 32 USB board). Tree supports.")],
          open("stl/parts.json", "w"), indent=1)
print("parts ok; camera at", (pc + n * 8).round(0), "head", (hi - lo).round(1).tolist())

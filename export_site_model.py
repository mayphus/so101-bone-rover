"""Export the rover for the mayphus.org SO-101 viewer (scripts/offline/prepare-bone-rover.py there).

Writes stl/site/*.stl and stl/site/manifest.json. Rover parts stay in rover coordinates (mm,
origin on the floor under the drive axle). The SO-101 itself is not exported: the website
converts it from the same pinned URDF and meshes and mounts it with `arm_mount`. The wrist
camera mount and camera are exported in the gripper_link frame so they move with the wrist.
usage: export_site_model.py [variant]   (default: print)
"""
import json, os, sys
import numpy as np
import trimesh
import params as P

VAR = sys.argv[1] if len(sys.argv) > 1 else "print"
OUT = "stl/site"
os.makedirs(OUT, exist_ok=True)
base = json.load(open("so101_base.json"))
arm_mount = trimesh.transformations.translation_matrix(
    [P.HUB_X - base["pan_xy"][0], -base["pan_xy"][1], P.Z_TOP - base["base_bottom_z"]])

FRAME = (0.81, 0.83, 0.84)
PARTS = [(f"stl/{VAR}/assembly/{n}.stl", n, FRAME) for n in P.PIECES] + [
    ("stl/assembly/_ref_gearbox.stl", "motor_gearbox", (0.91, 0.40, 0.11)),
    ("stl/assembly/_ref_motor.stl", "motor_can", (0.73, 0.75, 0.76)),
    ("stl/assembly/_ref_encoder.stl", "motor_encoder", (0.18, 0.56, 0.31)),
    ("stl/assembly/_ref_plug.stl", "encoder_plug", (0.95, 0.95, 0.94)),
    ("stl/assembly/_ref_tire.stl", "wheel_tire", (0.14, 0.16, 0.17)),
    ("stl/assembly/_ref_wheelhub.stl", "wheel_hub", (0.24, 0.29, 0.32)),
    ("stl/assembly/_ref_caster_fork.stl", "caster_fork", (0.54, 0.58, 0.60)),
    ("stl/assembly/_ref_caster_wheel.stl", "caster_wheel", (0.14, 0.16, 0.17)),
    ("stl/assembly/_ref_battery.stl", "battery", (0.50, 0.65, 0.79)),
    ("stl/assembly/_ref_driver.stl", "driver_board", (0.18, 0.44, 0.69)),
    ("stl/assembly/_ref_pi5.stl", "raspberry_pi", (0.24, 0.60, 0.37)),
    ("stl/assembly/_ref_mast.stl", "camera_mast", (0.16, 0.18, 0.20)),
    ("stl/assembly/camera_head.stl", "camera_head", FRAME),
    ("stl/assembly/_ref_camera.stl", "mast_camera", (0.16, 0.18, 0.20)),
]
visuals = []
for src, name, color in PARTS:
    m = trimesh.load(src)
    m.export(f"{OUT}/{name}.stl")
    visuals.append(dict(mesh=name, link="rover_base", transform_mm=np.eye(4).tolist(), color=color))

# wrist camera (same placement as arm_export.py), in the gripper_link frame
M = np.eye(4)
M[:3, :3] = [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]
M[:3, 3] = [-11.3, 26.0, -41.8]
mount = trimesh.load("vendor/so101_wristcam/SO-ARM101_camera_wrist_mount.stl")
mount.export(f"{OUT}/wrist_camera_mount.stl")
visuals.append(dict(mesh="wrist_camera_mount", link="gripper_link", transform_mm=M.tolist(), color=(0.88, 0.35, 0.16)))
pc = np.array([-48.5, -38.0, 17.5]); n = np.array([0.46, 0.89, 0.0]); n /= np.linalg.norm(n)
Rb = np.eye(4); Rb[:3, :3] = np.c_[[n[1], -n[0], 0.0], [0, 0, 1.0], n]; Rb[:3, 3] = pc
board = trimesh.creation.box((32, 32, 1.6)); board.apply_translation([0, 0, -3.5])
lens = trimesh.creation.cylinder(radius=7, height=14, sections=32); lens.apply_translation([0, 0, 4.0])
cam = trimesh.util.concatenate([board, lens]); cam.apply_transform(Rb)
cam.export(f"{OUT}/wrist_camera.stl")
visuals.append(dict(mesh="wrist_camera", link="gripper_link", transform_mm=M.tolist(), color=(0.16, 0.18, 0.20)))

json.dump(dict(variant=VAR, units="mm", arm_mount_mm=arm_mount.tolist(), visuals=visuals,
               wrist_mount_source="TheRobotStudio/SO-ARM100 Optional/SO101_Wrist_Cam_Hex-Nut_Mount_32x32_UVC_Module/stl/SO-ARM101_camera_wrist_mount.stl (main, Apache-2.0)"),
          open(f"{OUT}/manifest.json", "w"), indent=1)
print(len(visuals), "visuals ->", OUT)

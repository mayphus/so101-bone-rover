"""Place the official SO-101 on the frame's top pad, with the official wrist camera
mount (SO101_Wrist_Cam_Hex-Nut_Mount_32x32_UVC_Module) and a 32 x 32 USB camera."""
import json
import numpy as np
import trimesh
import params as P
import so101_arm as A

POSE = {"shoulder_pan": 0.0, "shoulder_lift": -0.35, "elbow_flex": 0.55, "wrist_flex": 0.6, "wrist_roll": 0.0, "gripper": 0.15}
base = json.load(open("so101_base.json"))
T = trimesh.transformations.translation_matrix([P.HUB_X - base["pan_xy"][0], -base["pan_xy"][1], P.Z_TOP - base["base_bottom_z"]])
asm = A.assemble(POSE)
TW = asm.pop("__T__")
printed, servos = [], []
for name, ms in asm.items():
    for servo, m in ms:
        m = m.copy(); m.apply_transform(T)
        (servos if servo else printed).append(m)

# wrist camera mount: its tab sits flat on the wrist roll part's flat side, its M3 hole on the
# part's M3 hole (gripper_link frame, mm); the camera plate trails behind, looking at the jaw tips.
M = np.eye(4)
M[:3, :3] = [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]       # mount x -> -y_g, y -> -z_g, z -> +x_g
M[:3, 3] = [-11.3, 26.0, -41.8]
G = T @ TW["gripper_link"].copy()
G[:3, 3] = G[:3, 3]                                    # URDF transforms are already in mm
mount = trimesh.load("vendor/so101_wristcam/SO-ARM101_camera_wrist_mount.stl")
mount.apply_transform(G @ M)
# 32 x 32 UVC board behind the angled camera plate, lens through it
pc = np.array([-48.5, -38.0, 17.5]); n = np.array([0.46, 0.89, 0.0]); n /= np.linalg.norm(n)
t1 = np.array([n[1], -n[0], 0.0]); t2 = np.array([0, 0, 1.0])
Rb = np.eye(4); Rb[:3, :3] = np.c_[t1, t2, n]
board = trimesh.creation.box((32, 32, 1.6)); board.apply_translation([0, 0, -3.5])
lens = trimesh.creation.cylinder(radius=7, height=14, sections=32); lens.apply_translation([0, 0, 4.0])
cam = trimesh.util.concatenate([board, lens])
Rb[:3, 3] = pc
cam.apply_transform(G @ M @ Rb)
mount.export("stl/assembly/_arm_wristmount.stl")
cam.export("stl/assembly/_arm_wristcam.stl")
for nm, group in (("_arm_printed", printed), ("_arm_servos", servos)):
    trimesh.util.concatenate(group).export(f"stl/assembly/{nm}.stl")
allm = trimesh.util.concatenate(printed + servos)
print("arm bounds", allm.bounds.round(0).tolist(), "mount", mount.bounds.round(0).tolist(), "cam", cam.bounds.round(0).tolist())

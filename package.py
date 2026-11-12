"""Collect the print-ready STLs of one frame variant into print_package/ and zip it.
usage: package.py [variant]   (default: print)"""
import os, shutil, sys
import trimesh

VAR = sys.argv[1] if len(sys.argv) > 1 else "print"
D = "print_package"
FILES = [("01_frame_rear.stl", f"stl/{VAR}/print/frame_rear.stl"),
         ("02_frame_front.stl", f"stl/{VAR}/print/frame_front.stl"),
         ("03_motor_pod_left.stl", f"stl/{VAR}/print/frame_pod_left.stl"),
         ("04_motor_pod_right.stl", f"stl/{VAR}/print/frame_pod_right.stl"),
         ("05_camera_head.stl", "stl/parts/camera_head.stl"),
         ("06_so101_wrist_camera_mount_OFFICIAL.stl", "vendor/so101_wristcam/SO-ARM101_camera_wrist_mount.stl")]
os.makedirs(D, exist_ok=True)
for dst, src in FILES:
    m = trimesh.load(src)
    assert m.is_watertight, f"{src} is not watertight"
    assert all(e <= 180 for e in m.extents), f"{src} does not fit the A1 mini: {m.extents.round(1)}"
    shutil.copy(src, f"{D}/{dst}")
    print(f"{dst:44s} {m.extents.round(0).tolist()}")
shutil.copy("vendor/so101/LICENSE", f"{D}/LICENSE_SO-ARM100_Apache-2.0.txt")
shutil.make_archive("bone_rover_print_package", "zip", ".", D)
print("-> bone_rover_print_package.zip (README.md in print_package/ is edited by hand)")

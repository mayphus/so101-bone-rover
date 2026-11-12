# Bone Rover: Print-friendly frame, print package

SO-101 mobile base for a Bambu Lab A1 mini (180 × 180 × 180 mm). Every file below fits the bed and
is already rotated the way it should be printed. Import into Bambu Studio and **do not rotate**.

## Prints

| # | File | Size (mm) | Orientation (already set) | Supports |
|---|------|-----------|---------------------------|----------|
| 01 | 01_frame_rear.stl | 173 × 160 × 89 | Upright, as it sits on the robot (it was optimized for this) | Tree (auto) |
| 02 | 02_frame_front.stl | 75 × 141 × 124 | Standing on its flat cut face | Tree (auto) |
| 03 | 03_motor_pod_left.stl | 83 × 74 × 17 | Motor face plate flat on the bed | Tree (auto), few |
| 04 | 04_motor_pod_right.stl | 81 × 74 × 17 | Motor face plate flat on the bed | Tree (auto), few |
| 05 | 05_camera_head.stl | 68 × 40 × 38 | As is | Tree (auto) |
| 06 | 06_so101_wrist_camera_mount_OFFICIAL.stl | 66 × 43 × 35 | As is | None |

File 06 is the official TheRobotStudio part (SO-ARM100 repo, `Optional/SO101_Wrist_Cam_Hex-Nut_Mount_32x32_UVC_Module`), Apache-2.0, license included.

## Slicer settings (Bambu Studio)

- **Material:** PETG (best) or PLA+. PETG handles heat from motors and sun outdoors.
- **Nozzle / layer:** 0.4 mm, 0.2 mm layers.
- **Walls:** 3. The heat-set insert holes need solid plastic around them.
- **Infill:** 15% gyroid. The FEA shows at least 40× strength margin; more infill only adds weight.
- **Supports:** tree (auto), threshold 45°. Paint-block supports inside the insert holes and motor bores if the slicer puts any there.
- **Brim:** 5 mm on 01 (rear body). Its floor contact is small for its height.
- **Plastic:** about 0.43 kg for the 4 frame pieces, plus about 30 g for 05 and 06.

## Hardware

| Qty | Part | Where |
|-----|------|-------|
| 22 | M3 heat-set insert (about 4.2 mm hole, 5.7 mm long) | 18 at the knuckle joints (smaller holes on the mating faces), 4 in the driver bosses |
| 18 | M3 × 10 socket head | Knuckle joints. Not longer: the hole is 11.5 mm deep. |
| 12 | M3 × 6 | Motor gearbox faces (6 per motor) |
| 4 | M5 heat-set insert (6.4 mm hole) + 4 M5 × 16 countersunk | SO-101 base onto the arm pad, through the base's own holes |
| 4 | M4 heat-set insert (5.6 mm hole) + 4 M4 × 8 | Caster plate. Pad takes 54 × 24 or 38 × 38 hole patterns. |
| 4 | M3 × 6 standoff (or 6 mm spacer) + 4 M3 × 6 | Waveshare driver on the front deck |
| 1 | M3 × 6 grub screw | Locks the mast tube in its socket |
| 1 | M3 × 16 + nut | Camera head pinch clamp |
| 2 | M3 × 8 + 2 M3 hex nut | Wrist camera mount (official) |
| 4 | M2 × 6 | Cameras onto the wrist mount and the camera head (spares from the Feetech servos fit) |

## Bought parts

- 2 × JGB37-520 12 V with Hall encoder, 1:90
- 2 × 100 mm rubber wheel with 6 mm hex coupling
- 1 × 2" (50 mm) rubber swivel caster, top plate, about 64 mm tall
- 1 × Waveshare General Driver for Robots (ESP32)
- 2 × 32 × 32 mm USB (UVC) camera (wide-angle one for the mast)
- 1 × 10 mm carbon or aluminium tube, 380 mm
- 1 × 3S pack up to 80 × 78 × 38 mm + 2 Velcro straps
- Optional: Raspberry Pi 5 on M2.5 standoffs above the driver

## Before printing, measure your parts

These depend on the seller. Change them in `params.py` and rebuild if yours differ:
- Motor: gearbox length (`GEARBOX_L`, 26.5 mm for 1:90), face hole angle (`FACE_HOLE_ANGLE0`), encoder plug position.
- Wheel: 100 × 40 mm (`WHEEL_D`, `WHEEL_W`).
- Caster height: 64 mm (`CASTER_H`). A taller caster tilts the base back.

## Print order

1. Print 03 (one motor pod) first. It is small and quick, and it checks the motor face holes and gearbox fit before you spend about a day on 01.
2. Then 01, 02, 04, 05, 06.

## Assembly (short)

1. Press in the inserts with a soldering iron: M3 in the knuckle faces and driver bosses, M5 in the arm pad, M4 in the caster bosses.
2. Motors into the pods (plug the encoder cables in first), 6 × M3 × 6 each.
3. Pods and front body onto the rear body: 18 × M3 × 10 through the knuckles.
4. Wheels on the D-shafts. Caster on with 4 × M4 × 8.
5. Battery strapped under the arm. Driver on its standoffs.
6. SO-101 onto the pad with 4 × M5 countersunk from inside its base.
7. Mast into the tail socket (grub screw), camera head on top.
8. Wrist camera: remove gripper servo (motor 6), drop 2 M3 nuts into the wrist roll part, refit, screw on mount 06, then camera with 4 × M2.

Keep the arm working in front of the base. Reaching fully backward with a load can tip it over the rear wheels.

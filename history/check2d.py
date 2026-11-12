import runpy, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import trimesh
g = runpy.run_path("build.py")
d2 = g["deck2d"]; polys = g["polys"]
fig, ax = plt.subplots(figsize=(12, 6.5))
for p in polys(d2):
    ax.fill(*p.exterior.xy, color="#334155")
    for r in p.interiors: ax.fill(*r.xy, color="white")
for xs in (80, -80): ax.axvline(xs, color="r", lw=0.6, ls="--")
ax.set_aspect("equal"); ax.set_title("deck truss (red = print seams)")
plt.savefig("renders/deck_2d.png", dpi=90, bbox_inches="tight")
for n in ["deck_center","deck_front","deck_rear","wheel_hub","motor_pod_left","caster_body","tire_tpu"]:
    m = trimesh.load(f"stl/print/{n}.stl"); print(n, "bodies:", len(m.split(only_watertight=False)))

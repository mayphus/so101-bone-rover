import glob, sys, numpy as np, trimesh, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
COL = {"deck":"#e8eef5","pod":"#3b82f6","wheel":"#f97316","tire":"#222","caster":"#3b82f6",
       "splice":"#3b82f6","servo_strap":"#3b82f6","_proxy_servos":"#555","_proxy_battery":"#9ca3af",
       "_proxy_balls":"#c0c0c0","_proxy_pi5":"#16a34a"}
def color(n):
    for k,v in COL.items():
        if n.startswith(k) or k in n: return v
    return "#ccc"
def draw(files, out, elev, azim, lim=None):
    fig = plt.figure(figsize=(12,8)); ax = fig.add_subplot(111, projection="3d")
    L = np.array([0.3,-0.4,0.85]); L/=np.linalg.norm(L)
    allv=[]
    for f in files:
        n = f.split("/")[-1][:-4]; m = trimesh.load(f)
        if len(m.faces) > 60000:
            m = m.simplify_quadric_decimation(face_count=60000) if hasattr(m,'simplify_quadric_decimation') else m
        tri = m.vertices[m.faces]; nrm = m.face_normals
        sh = 0.45 + 0.55*np.clip(nrm@L,0,1)
        base = np.array(matplotlib.colors.to_rgb(color(n)))
        fc = np.clip(base[None,:]*sh[:,None],0,1)
        ax.add_collection3d(Poly3DCollection(tri, facecolors=fc, edgecolors="none", linewidths=0))
        allv.append(m.vertices)
    v=np.vstack(allv); c=(v.max(0)+v.min(0))/2; r=(v.max(0)-v.min(0)).max()/2 if lim is None else lim
    ax.set_xlim(c[0]-r,c[0]+r); ax.set_ylim(c[1]-r,c[1]+r); ax.set_zlim(c[2]-r,c[2]+r)
    ax.set_box_aspect((1,1,1)); ax.view_init(elev, azim); ax.axis("off")
    plt.savefig(out, dpi=90, bbox_inches="tight", pad_inches=0); plt.close()
files = sorted(glob.glob("stl/assembly/*.stl"))
mode = sys.argv[1]
if mode=="iso": draw(files, "renders/assembly_iso.png", 22, -55)
if mode=="top": draw([f for f in files if "deck" in f], "renders/deck_top.png", 90, -90)
if mode=="under": draw(files, "renders/assembly_under.png", -25, -120)
if mode=="parts":
    for n in ["wheel_hub","motor_pod_left","caster_body"]:
        draw([f"stl/print/{n}.stl"], f"renders/{n}.png", 35, -50)

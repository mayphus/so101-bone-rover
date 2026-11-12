"""Convert the viewer's assembly STLs to compact indexed-mesh JSON (web/...json).
Published pages can't serve .stl; the full-resolution print STLs stay in stl/*/print."""
import glob, json, os
import numpy as np
import trimesh
from manifold3d import Manifold, Mesh

for src in sorted(glob.glob("stl/assembly/*.stl") + glob.glob("stl/*/assembly/*.stl")):
    dst = "web/" + src[len("stl/"):].replace(".stl", ".json")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    t = trimesh.load(src, process=True)
    try:
        m = Manifold(Mesh(vert_properties=np.asarray(t.vertices, np.float32),
                          tri_verts=np.asarray(t.faces, np.uint32))).simplify(0.06)
        o = m.to_mesh()
        v, f = np.asarray(o.vert_properties)[:, :3], np.asarray(o.tri_verts)
        if len(f) == 0:
            raise ValueError
    except Exception:
        v, f = t.vertices, t.faces
    q = np.round(v * 100).astype(np.int32)           # 0.01 mm grid
    json.dump({"s": 0.01, "v": q.ravel().tolist(), "i": f.astype(np.int64).ravel().tolist()},
              open(dst, "w"), separators=(",", ":"))
    print(f"{dst:42s} tris {len(t.faces):7d} -> {len(f):7d}  {os.path.getsize(dst)/1e6:.1f} MB")

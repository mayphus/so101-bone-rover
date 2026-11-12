"""Assemble the official SO-101 follower (TheRobotStudio/SO-ARM100, Apache-2.0)
from its URDF in a chosen pose. Returns one trimesh per link group, in mm,
in the URDF base_link frame."""
import math
import xml.etree.ElementTree as ET
import numpy as np
import trimesh

ROOT = "vendor/so101"
URDF = f"{ROOT}/so101_new_calib.urdf"

def rpy_matrix(r, p, y):
    cr, sr, cp, sp, cy, sy = math.cos(r), math.sin(r), math.cos(p), math.sin(p), math.cos(y), math.sin(y)
    Rx = np.array([[1, 0, 0], [0, cr, -sr], [0, sr, cr]])
    Ry = np.array([[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]])
    Rz = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx

def tf(origin):
    T = np.eye(4)
    if origin is not None:
        xyz = [float(v) for v in origin.get("xyz", "0 0 0").split()]
        rpy = [float(v) for v in origin.get("rpy", "0 0 0").split()]
        T[:3, :3] = rpy_matrix(*rpy)
        T[:3, 3] = np.array(xyz) * 1000.0
    return T

def axis_rot(axis, q):
    a = np.array(axis, float); a /= np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    T = np.eye(4)
    T[:3, :3] = np.eye(3) + math.sin(q) * K + (1 - math.cos(q)) * K @ K
    return T

def assemble(pose):
    root = ET.parse(URDF).getroot()
    links = {l.get("name"): l for l in root.findall("link")}
    joints = root.findall("joint")
    child_of = {j.find("child").get("link"): j for j in joints}
    world = {}

    def T_of(link):
        if link in world:
            return world[link]
        j = child_of.get(link)
        if j is None:
            world[link] = np.eye(4)
        else:
            T = T_of(j.find("parent").get("link")) @ tf(j.find("origin"))
            if j.get("type") in ("revolute", "continuous"):
                ax = [float(v) for v in j.find("axis").get("xyz").split()]
                T = T @ axis_rot(ax, pose.get(j.get("name"), 0.0))
            world[link] = T
        return world[link]

    out = {}
    for name, l in links.items():
        meshes = []
        for vis in l.findall("visual"):
            mesh = vis.find("geometry/mesh")
            if mesh is None:
                continue
            m = trimesh.load(f"{ROOT}/{mesh.get('filename')}", force="mesh")
            m.apply_scale(1000.0)                       # URDF meshes are in metres
            m.apply_transform(T_of(name) @ tf(vis.find("origin")))
            servo = "sts3215" in mesh.get("filename")
            meshes.append((servo, m))
        if meshes:
            out[name] = meshes
    out["__T__"] = {n: T_of(n) for n in links}
    return out

"""Build the drive, wheel and caster parts of the SO-101 2-motor base and export STLs.

stl/print/*.stl     parts laid flat in print orientation (Bambu A1 mini)
stl/assembly/*.stl  parts placed in the robot frame (for the 3D viewer)
"""
import json, math, os
import numpy as np
import trimesh
from manifold3d import Manifold, CrossSection, FillRule
from shapely.geometry import Polygon, Point, LineString, MultiPolygon, box
from shapely.geometry.polygon import orient
from shapely.ops import unary_union
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from skimage.measure import find_contours
import params as P

os.makedirs("stl/print", exist_ok=True)
os.makedirs("stl/assembly", exist_ok=True)
SEG = 48

# ---------------------------------------------------------------- helpers
def polys(g):
    if g.is_empty:
        return []
    return list(g.geoms) if hasattr(g, "geoms") else [g]

def to_cs(g):
    rings = []
    for p in polys(g):
        p = orient(p, 1.0)
        rings.append(np.asarray(p.exterior.coords)[:-1])
        rings += [np.asarray(r.coords)[:-1] for r in p.interiors]
    return CrossSection(rings, FillRule.EvenOdd)

def extrude(g, h, z0=0.0):
    return Manifold.extrude(to_cs(g), h).translate((0, 0, z0))

def cyl(r, h, x=0, y=0, z=0, r2=None, seg=SEG):
    return Manifold.cylinder(h, r, r if r2 is None else r2, seg).translate((x, y, z))

def cube(x0, x1, y0, y1, z0, z1):
    return Manifold.cube((x1 - x0, y1 - y0, z1 - z0)).translate((x0, y0, z0))

def hexprism(af, h, x, y, z):
    r = af / math.sqrt(3)
    return Manifold.cylinder(h, r, r, 6).translate((x, y, z))

def union(ms):
    return Manifold.batch_boolean(ms, __import__("manifold3d").OpType.Add)

def holes_2d(pts, d):
    return unary_union([Point(p).buffer(d / 2, 24) for p in pts])

def to_trimesh(m):
    mesh = m.to_mesh()
    return trimesh.Trimesh(np.asarray(mesh.vert_properties)[:, :3],
                           np.asarray(mesh.tri_verts), process=False)

REPORT = []

def export(name, m_assembly, m_print, qty=1, material="PLA/PETG", note=""):
    to_trimesh(m_assembly).export(f"stl/assembly/{name}.stl")
    tm = to_trimesh(m_print)
    lo, hi = tm.bounds
    tm.apply_translation([-(lo[0] + hi[0]) / 2, -(lo[1] + hi[1]) / 2, -lo[2]])
    tm.export(f"stl/print/{name}.stl")
    ext = (hi - lo).round(1).tolist()
    fits = ext[0] <= P.BED and ext[1] <= P.BED and ext[2] <= P.BED
    vol = m_print.volume() / 1000.0
    REPORT.append(dict(name=name, qty=qty, material=material, size_mm=ext,
                       fits_a1_mini=bool(fits), volume_cm3=round(vol, 1),
                       est_grams=round(vol * 1.25 * 0.55, 0), note=note,
                       watertight=bool(tm.is_watertight)))
    print(f"{name:18s} {ext}  fits={fits}  {vol:.0f} cm3  watertight={tm.is_watertight}")

Z0 = P.DECK_Z0

# ================================================================ MOTOR POD
c = P.SERVO_CLEAR
sx0, sx1 = P.servo_x_span()
sx0, sx1 = sx0 - c, sx1 + c
zb = P.WHEEL_R - P.SERVO_W / 2 - c           # servo bottom
zt = P.WHEEL_R + P.SERVO_W / 2 + c           # servo top
yi = P.SERVO_FACE_Y - P.SERVO_H              # servo inner face
yo = P.SERVO_FACE_Y                          # servo output face
W = 3.0
Z0 = P.DECK_Z0

def motor_pod():
    parts = []
    # top plate / flange
    fl = LineString([(P.POD_BOLTS_X[0], 66.5), (P.POD_BOLTS_X[1], 66.5)]).buffer(19, 32)
    fl = fl.intersection(box(-200, yi - 3, 200, yo + 2.5))
    fl = fl.union(box(sx0 - W, yi - 3, sx1 + W, yo + 2.5))
    parts.append(extrude(fl, Z0 - zt, zt))
    # end walls along x (grip the servo's short ends)
    parts.append(cube(sx0 - W, sx0, yi + 2, yo + 2.5, zb - 0.01, zt + 0.01))
    parts.append(cube(sx1, sx1 + W, yi + 2, yo + 2.5, zb - 0.01, zt + 0.01))
    # outer lip on the rear part of the output face (keeps horn side free)
    parts.append(cube(sx0 - W, -P.HORN_D / 2 - 2.0, yo, yo + 2.5, zb - 0.01, zt))
    # inner tab: top half of the inner face, servo drops in from below
    parts.append(cube(sx0 - W, sx1 + W, yi - 3, yi, (zb + zt) / 2 + 2, zt))
    # bolt bosses (nut traps underneath) + strap screw bosses
    for bx in P.POD_BOLTS_X:
        for by in P.POD_BOLTS_Y:
            parts.append(cyl(6.0, Z0 - (Z0 - 9), bx, by, Z0 - 9))
    ym = (yi + yo) / 2
    for bx in (sx0 - W - 3.5, sx1 + W + 3.5):
        parts.append(cyl(4.5, Z0 - zb + 0.01, bx, ym, zb - 0.01))
        parts.append(cube(min(bx, (sx0 - W if bx < 0 else sx1 + W)) - 0.01,
                          max(bx, (sx0 - W if bx < 0 else sx1 + W)) + 0.01,
                          ym - 4.5, ym + 4.5, zb - 0.01, Z0))
    # organic gussets from bosses to the sleeve
    for bx in P.POD_BOLTS_X:
        ex = sx0 - W if bx < 0 else sx1 + W
        g = Manifold.hull(union([cyl(3.0, 0.1, bx, P.POD_BOLTS_Y[0], Z0 - 9),
                                 cyl(3.0, 0.1, bx, P.POD_BOLTS_Y[1], Z0 - 9),
                                 cube(ex - 0.5, ex + 0.5, P.POD_BOLTS_Y[0], P.POD_BOLTS_Y[1], zb + 6, zb + 6.1),
                                 cube(min(bx, ex), max(bx, ex), P.POD_BOLTS_Y[0] - 2, P.POD_BOLTS_Y[1] + 2, Z0 - 0.5, Z0)]))
        g = g ^ cube(min(bx, ex) - 6, max(bx, ex) + 6, yi - 3, yo + 2.5, 0, Z0)
        parts.append(g)
    pod = union(parts)
    # servo pocket + horn clearance
    pod -= cube(sx0, sx1, yi, yo + 0.01, zb - 20, zt)
    pod -= cyl(P.HORN_D / 2 + 1.5, 10, 0, 0, 0).rotate((-90, 0, 0)).translate((0, yo - 1, P.WHEEL_R))
    # cable slot in inner tab, wire exit
    pod -= cube(-14, -2, yi - 4, yi + 0.01, (zb + zt) / 2, zt - 3)
    for bx in P.POD_BOLTS_X:
        for by in P.POD_BOLTS_Y:
            pod -= cyl(P.M3 / 2, 20, bx, by, Z0 - 15)
            pod -= hexprism(P.M3_NUT_AF, 3.0, bx, by, Z0 - 9.01)
    for bx in (sx0 - W - 3.5, sx1 + W + 3.5):
        pod -= cyl(P.M3_TAP / 2, 12, bx, ym, zb - 0.1)
    return pod

pod_r = motor_pod()                       # +y side (left wheel)
pod_l = pod_r.mirror((0, 1, 0))
flip = lambda m: m.rotate((180, 0, 0))
export("motor_pod_left", pod_r, flip(pod_r), note="Flange down. Mirror pair.")
export("motor_pod_right", pod_l, flip(pod_l), note="Flange down. Mirror pair.")

# bottom strap that holds the servo up in the pod
strap2d = LineString([(sx0 - W - 3.5, 0), (sx1 + W + 3.5, 0)]).buffer(5, 24)
strap2d = strap2d.difference(holes_2d([(sx0 - W - 3.5, 0), (sx1 + W + 3.5, 0)], P.M3))
strap2d = strap2d.difference(LineString([(sx0 + 6, 0), (sx1 - 6, 0)]).buffer(1.6, 12))
strap = extrude(strap2d, 3.0)
ym = (yi + yo) / 2
export("servo_strap", strap.translate((0, ym, zb - 3.0)), strap, qty=2,
       note="Two needed. M3x16 self-tapping into the pod bosses.")

# ================================================================ WHEEL HUB
def wheel_hub():
    rim_in = P.RIM_R - 4.0
    rim = cyl(P.RIM_R, P.WHEEL_W) - cyl(rim_in, P.WHEEL_W + 2, z=-1)
    for z in (0, P.WHEEL_W - 1.6):
        rim += cyl(P.RIM_R + 1.6, 1.6, z=z) - cyl(rim_in, 3, z=z - 0.5)
    rim -= cyl(P.RIM_R + 0.01, 1.2, z=P.WHEEL_W / 2 - 0.6) - cyl(P.RIM_R - 1.0, 1.2, z=P.WHEEL_W / 2 - 0.6)  # tire key groove
    # swept, tapered spokes
    sp = []
    for k in range(5):
        th0 = 2 * math.pi * k / 5
        circles = []
        for t in np.linspace(0, 1, 14):
            r = 14 + t * (rim_in + 1 - 14)
            th = th0 + 0.75 * t ** 1.3
            w = 5.0 - 2.2 * t
            circles.append(Point(r * math.cos(th), r * math.sin(th)).buffer(w, 20))
        sp.append(unary_union([circles[i].union(circles[i + 1]).convex_hull for i in range(len(circles) - 1)]))
    sp2d = unary_union(sp + [Point(0, 0).buffer(17, 64)])
    spokes = extrude(sp2d, 9.0).refine_to_length(2.0)

    def dish(v):
        v = np.array(v, dtype=float)
        r = np.hypot(v[:, 0], v[:, 1])
        s = np.clip((r - 17) / (rim_in - 17), 0, 1)
        v[:, 2] += 7.5 * (3 * s ** 2 - 2 * s ** 3)
        return v
    spokes = spokes.warp_batch(dish).simplify(0.03)
    boss = cyl(17, 12.0) + cyl(14.5, 2.0, z=12.0, r2=12.0)
    hub = union([rim, spokes, boss])
    hub -= cyl(P.HORN_D / 2 + 0.15, P.HORN_DEPTH + 0.01, z=-0.01)
    hub -= cyl(4.0, 30, z=-1)                           # horn center screw access
    for k in range(4):
        a = math.pi / 4 + k * math.pi / 2
        x, y = P.HORN_BOLT_CIRCLE / 2 * math.cos(a), P.HORN_BOLT_CIRCLE / 2 * math.sin(a)
        hub -= cyl(P.HORN_BOLT_D / 2, 30, x, y, -1)
        hub -= cyl(2.2, 10, x, y, 9.0)                  # head counterbore
    return hub

hub = wheel_hub()
yw = P.SERVO_FACE_Y + P.HORN_PROTRUDE + 0.3
place_l = hub.rotate((-90, 0, 0)).translate((0, yw, P.WHEEL_R))
place_r = hub.rotate((90, 0, 0)).translate((0, -yw, P.WHEEL_R))
export("wheel_hub", union([place_l, place_r]), hub, qty=2,
       note="Two needed. Horn side down. 4 walls, 40% infill.")

# ================================================================ TIRE (TPU)
def tire():
    w = P.WHEEL_W - 3.2 - 0.3
    t = cyl(P.WHEEL_R, w, seg=128) - cyl(P.RIM_R - 0.4, w + 2, z=-1, seg=128)
    t += cyl(P.RIM_R - 0.4, 1.0, z=w / 2 - 0.5, seg=128) - cyl(P.RIM_R - 1.4, 2, z=w / 2 - 1, seg=128)
    grooves = []
    n = 40
    for k in range(n):
        a = 360.0 * k / n
        for sgn in (1, -1):
            g = cube(-1.0, 1.0, -0.5 * w * 0.62, 0.5 * w * 0.62, -2, 2)
            g = g.rotate((0, 0, sgn * 32)).rotate((0, 90, 0))
            g = g.translate((P.WHEEL_R, sgn * 3.2, w / 2 + sgn * w * 0.22)) if False else \
                cube(-2, 2, -1.0, 1.0, -w * 0.3, w * 0.3).rotate((sgn * 32, 0, 0)).translate(
                    (P.WHEEL_R, 0, w / 2 + sgn * w * 0.24))
            grooves.append(g.rotate((0, 0, a + (360 / n / 2 if sgn < 0 else 0))))
    return t - union(grooves)

tr = tire()
export("tire_tpu", union([tr.translate((0, 0, -tr.bounding_box()[2])).rotate((-90, 0, 0)).translate((0, yw + 1.9, P.WHEEL_R)),
                          tr.rotate((90, 0, 0)).translate((0, -yw - 1.9, P.WHEEL_R))]),
       tr, qty=2, material="TPU 95A",
       note="Two needed. TPU 95A, 3 walls, 15% infill, slow (25 mm/s). Stretch over the rim.")

# ================================================================ BALL CASTER
zc = P.BALL_D / 2 + P.CASTER_LIFT
rb = P.BALL_D / 2

def caster_body():
    body = cyl(19.0, Z0 - 4 - zc, z=zc, r2=12.5, seg=64)
    body += cyl(P.CASTER_BOLT_R + 7, 4.0, z=Z0 - 4, seg=64)
    for b in P.caster_bolts(0):
        body += cyl(4.5, 12, b[0], b[1], Z0 - 12)
        body += Manifold.hull(union([cyl(2.5, 0.1, b[0], b[1], Z0 - 12),
                                     cyl(0.8, 0.1, b[0] * 0.55, b[1] * 0.55, zc + 18),
                                     cyl(4.5, 0.1, b[0], b[1], Z0 - 0.1)]))
    body -= Manifold.sphere(rb + 0.3, 96).translate((0, 0, zc))
    body -= cyl(7.0, Z0 - zc - 16, z=zc + 10)              # lighten
    for k in range(3):                                       # windows
        a = 60 + 120 * k
        w = Manifold.hull(union([Manifold.sphere(4.5, 24).translate((0, 0, zc + 24)),
                                 Manifold.sphere(3.0, 24).translate((0, 0, Z0 - 13))]))
        w = w.translate((15, 0, 0)).rotate((0, 0, a))
        body -= w
    for b in P.caster_bolts(0):
        body -= cyl(P.M3 / 2, 20, b[0], b[1], Z0 - 15)
        body -= hexprism(P.M3_NUT_AF, 3.0, b[0], b[1], Z0 - 12.01)
    cap_bolts = [(15.5 * math.cos(math.radians(a)), 15.5 * math.sin(math.radians(a))) for a in (60, 180, 300)]
    for b in cap_bolts:
        body -= cyl(P.M3_TAP / 2, 12, b[0], b[1], zc - 0.1)
    return body, cap_bolts

def caster_cap(cap_bolts):
    zb_ = zc - math.sqrt(rb ** 2 - 11.5 ** 2)
    cap = cyl(19.0, zc - zb_, z=zb_, seg=64)
    cap -= Manifold.sphere(rb + 0.3, 96).translate((0, 0, zc))
    cap -= cyl(11.5, 10, z=zb_ - 1)
    for b in cap_bolts:
        cap -= cyl(P.M3 / 2, 20, b[0], b[1], zb_ - 1)
        cap -= cyl(3.0, 2.5, b[0], b[1], zb_ - 0.01)
    return cap

cb, cbolts = caster_body()
cc = caster_cap(cbolts)
export("caster_body", union([cb.translate((P.CASTER_X, 0, 0)), cb.translate((-P.CASTER_X, 0, 0))]),
       flip(cb), qty=2, note="Two needed. Flange down.")
export("caster_cap", union([cc.translate((P.CASTER_X, 0, 0)), cc.translate((-P.CASTER_X, 0, 0))]),
       cc, qty=2, note="Two needed. Opening down.")

# ================================================================ PROXIES (viewer only)
def proxy(name, m):
    to_trimesh(m).export(f"stl/assembly/_proxy_{name}.stl")

servos = []
for s in (1, -1):
    sv = cube(sx0 + c, sx1 - c, yi, yo, zb + c, zt - c)
    sv += cyl(P.HORN_D / 2, P.HORN_PROTRUDE, 0, 0, 0).rotate((-90, 0, 0)).translate((0, yo, P.WHEEL_R))
    servos.append(sv if s > 0 else sv.mirror((0, 1, 0)))
proxy("servos", union(servos))
proxy("battery", cube(-60, 60, -40, 40, 18, Z0 - 0.5).refine_to_length(200))
balls = union([Manifold.sphere(rb, 48).translate((sx, 0, zc)) for sx in (P.CASTER_X, -P.CASTER_X)])
proxy("balls", balls)
pi = cube(P.PI_X - 42.5, P.PI_X + 42.5, -28, 28, 96 + 6, 96 + 7.6) + \
     cube(P.PI_X - 30, P.PI_X + 10, -22, 22, 96 + 7.6, 96 + 18)
proxy("pi5", pi)

json.dump(REPORT, open("stl/parts.json", "w"), indent=1)

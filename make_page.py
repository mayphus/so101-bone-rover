"""Fill page_template.html with the build reports -> index.html."""
import json
import params as P

import os
INFO = [("print", "Print-friendly", "13% of the design space, 45° overhang rule. Least plastic and supports."),
        ("stiff", "Stiff", "18%, same overhang rule. Thicker bones for heavier payloads and rough ground.")]
variants = []
for key, title, why in INFO:
    if not os.path.exists(f"stl/{key}/summary.json"):
        if key == "light":
            variants.append(dict(key=key, title=title, why=why, rejected=True, frame=[],
                                 reason="Dropped: at 9% material the optimizer never keeps the member to the front caster, which a real robot needs when the arm reaches forward."))
        else:
            variants.append(dict(key=key, title=title, why=why, pending=True, frame=[], reason="Still optimizing."))
        continue
    sm = json.load(open(f"stl/{key}/summary.json"))
    if sm.get("islands", 1) > 1:
        variants.append(dict(key=key, title=title, why=why, rejected=True, grams=sm["grams"], frame=[],
                             reason="Rejected: at this budget the structure breaks into separate islands."))
        continue
    fea = json.load(open(f"stl/{key}/fea.json")) if os.path.exists(f"stl/{key}/fea.json") else None
    variants.append(dict(key=key, title=title, why=why, grams=sm["grams"], compliance=sm["compliance"], fea=fea,
                         overhang_cm2=sm["overhang_cm2"], frame=json.load(open(f"stl/{key}/frame.json"))))
parts = json.load(open("stl/parts.json"))
chassis = variants[0]["frame"]
allp = chassis + parts
grams = sum(p["est_grams"] * p["qty"] for p in allp)
n_prints = sum(p["qty"] for p in allp)
track = 2 * P.WHEEL_Y1
rpm = 110.0                                  # JGB37-520 12 V, 1:90, no load
speed = rpm / 60 * 3.1416 * P.WHEEL_D / 1000
stats = [
    [f"{P.X_NOSE - P.X_TAIL + 10:.0f} × {track:.0f}", "footprint (mm), wheels and caster included"],
    [f"{P.WHEEL_D:.0f} mm", "treaded drive wheels + 2\" caster"],
    [f"~{speed:.2f} m/s", "top speed at 12 V, no load"],
    [f"~{grams/1000:.2f} kg", "printed plastic (Print-friendly), est."],
    [f"{n_prints}", "prints, all fit 180 mm"],
    ["2", "cameras: wrist + mast, no depth camera"],
]
html = open("page_template.html").read()
html = html.replace("__PARTS_JSON__", json.dumps(parts)).replace("__VARIANTS_JSON__", json.dumps(variants)).replace("__STATS_JSON__", json.dumps(stats))
open("index.html", "w").write(html)
print(stats)

---
name: bone-rover
description: Regenerate, check and ship the SO-101 Bone Rover printed frame. Use when the user wants to change a dimension or part (motor, wheel, caster, arm position, cameras, electronics), re-run the topology optimization, rebuild the frame, run FEA, fix a printability or shape issue (thin spots, tails, joints), update the 3D viewer page, rebuild the print package, or slice for the Bambu A1 mini.
---

# Bone Rover pipeline

Everything lives in this repo. All lengths are mm; frame: x forward, y left, z up, origin on the floor under the drive axle.
Use `.venv/bin/python` (create it with `uv venv .venv --python 3.12 && uv pip install --python .venv/bin/python -r requirements.txt`).

## Where things are
- `params.py`: every dimension (motor, wheel, caster, arm `HUB_X`, battery bay, driver `DRV_*`, mast, print cuts `SPLIT_X`/`SPLIT_Y`, joints `COLLAR_*`, `PIECES`). Values marked VERIFY depend on the part the user buys.
- `layout.py`: `keepout()` (where material may not go) and `solid()` (exact mounting features). Shared by optimizer and builder.
- `topo.py`: 3D SIMP. `topo.py <volfrac> <iterations> <out.npz> [am]`. Variants: print = `0.13 50 topo_print.npz am`, stiff = `0.18 50 topo_stiff.npz am`. Env: `GRID_H` (default 6), `INIT=<npz>` warm start (use with ~30 iterations), `K_BOSS`.
- `structure.py <npz> <variant>`: field → isosurface → exact features → cuts → knuckle joints → 4 pieces in `stl/<variant>/{assembly,print}/` + `frame.json`, `summary.json`. `REFINE=1` uses `stl/<variant>/stress_field.npz` from `verify.py`.
- `parts.py` (reference models + camera head), `arm_export.py` (official SO-101 + wrist camera via URDF in `vendor/`).
- `export_web.py` → `web/**.json` meshes; `verify.py <variant>` → FEA, `stl/<variant>/fea.json`, stress maps; `make_page.py` → `index.html`; `package.py <variant>` → `print_package/` + zip.

## Standard rebuild (after changing params/layout)
1. If keep-outs, features or loads changed: re-run `topo.py` for both variants in the background (`OMP_NUM_THREADS=4`, ~45–75 min each, run both in parallel). Geometry-only tweaks that don't move features (joints, smoothing, cut position) don't need it.
2. `parts.py`, `arm_export.py`
3. For each variant: `structure.py` → `export_web.py` → `verify.py` → `REFINE=1 structure.py` → `export_web.py` → `verify.py`
4. `make_page.py`, `package.py print`, update `print_package/README.md` numbers by hand.
5. Preview: `.claude/launch.json` config `rover` (http.server on 8765) → `index.html`.

## Checks that must pass before telling the user it's done
- `summary.json`: `islands == 1`, every piece `bodies == 1`, `fits == true`. Several bodies in a piece = a floating feature (it happened with driver bosses, the mast socket and front caster bosses). Fix the load path (a load case or spring that forces it), not by bolting on a block.
- `verify.py` stops with `FRAME NOT CONNECTED` if more than 40 voxels are loose. A few dropped specks are thin slivers; report them.
- FEA `cg_info` must be 0 and `max_disp_mm` small (a huge value means an unsupported loaded part).
- Stress JSON vertex count must match the mesh JSON (re-run `export_web.py` before `verify.py`).
- Thin spots: measure with a voxel opening (0.4 mm voxels, radius 0.8 / 1.25 mm) rather than eyeballing.

## Design rules learned the hard way
- The optimizer only grows material to loaded or supporting features. Light features (boards, mast) need enough load-case weight or they float. A support that is not loaded can be skipped: the caster plate is modelled as an equal push on every bolt boss (case 12) so all bosses get struts.
- No flat walls or blocky add-ons: the user rejects them. Joints are knuckles grown in the density field (`COLLAR_GROW`, `COLLAR_T`, cosine taper) with M3 screws + heat-set inserts; the cut stays flat. Pick cut planes where few members cross and never through the arm pad (`SPLIT_X` must stay ahead of x ≈ 82).
- The 45° overhang filter grows downward fins; `REFINE=1` removes dead material (<4% of peak stress, or <15% where thinner than 7 mm) and thickens load-carrying parts thinner than 5 mm.
- Export manifold3d meshes as they come (`process=False`); trimesh cleanup can open tiny holes.
- Stress colours in the viewer are relative to each frame's own peak: red is not "weak". Quote absolute MPa and safety factors.
- Tipping: reaching backward past the rear axle with load can tip the base.

## Printing
- A1 mini, PLA Basic or PETG, 0.2 mm layers, 3 walls (2 for quick tests), 15% gyroid, tree supports. Front body prints on its flat cut face; rear upright; pods face plate down. Print one pod first as a fit test.
- OrcaSlicer 2.5 is installed. Its CLI cannot resolve stub presets (bundled profiles are packed); it needs a base project 3MF saved from the GUI with A1 mini + PLA Basic + 0.20mm Standard, then slice that from the CLI with modified settings.
- Never take or use the printer's LAN access code (not even via an env var). The user sends jobs from OrcaSlicer, or runs a send script themselves.

## Page
The 3D viewer page is a private claude.ai Artifact; republish `index.html` with the changed `web/**` JSON files in `files` after every rebuild. Don't publish publicly without asking.

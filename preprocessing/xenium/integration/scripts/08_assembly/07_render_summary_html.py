"""Assemble tabbed summary.html from previews/ PNGs.

No re-rendering — reads PNGs already on disk and embeds them as base64.
Tabs: Embedding | Labels | Clinical | Quality | UOQ | Compartments

Self-contained single HTML file; no external dependencies.
"""
import argparse, base64, os
import yaml

TABS = [
    ("Embedding", ["umap_platform", "umap_patient", "umap_compartment"]),
    ("Labels",    ["umap_l0p5", "umap_position"]),
    ("Clinical",  ["umap_age", "umap_menopause", "umap_brca"]),
    ("Quality",   ["umap_nmp_exclude", "per_sample_qc"]),
    ("UOQ",       ["umap_p3_uoq"]),
    ("Compartments", ["umap_compartments_overview",
                      "umap_epithelial_l0p5", "umap_stromal_l0p5", "umap_immune_l0p5"]),
]


def img_b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()


def panel_html(name, previews_dir):
    png = os.path.join(previews_dir, name + ".png")
    pdf = os.path.join(previews_dir, name + ".pdf")
    if not os.path.exists(png):
        return f'<p style="color:#999">{name} — not yet rendered</p>'
    b64 = img_b64(png)
    pdf_link = (f'<a href="{name}.pdf" style="font-size:11px">PDF</a>'
                if os.path.exists(pdf) else "")
    return (f'<figure style="display:inline-block;margin:8px">'
            f'<img src="data:image/png;base64,{b64}" '
            f'style="max-width:100%;border:1px solid #ddd"><br>'
            f'<figcaption style="font-size:11px;text-align:center">'
            f'{name} {pdf_link}</figcaption></figure>')


def build_html(previews_dir):
    tab_ids = [t[0].replace(" ", "_") for t in TABS]

    buttons = "\n".join(
        f'<button class="tab-btn" onclick="show(\'{tid}\')" id="btn-{tid}">'
        f'{label}</button>'
        for (label, _), tid in zip(TABS, tab_ids)
    )

    sections = "\n".join(
        f'<div id="tab-{tid}" class="tab-section" style="display:none">'
        + "".join(panel_html(p, previews_dir) for p in panels)
        + "</div>"
        for (_, panels), tid in zip(TABS, tab_ids)
    )

    first = tab_ids[0]
    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Spatial HBCA Xenium — Pipeline Preview</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 16px; }}
  .tab-btn {{ padding: 6px 14px; margin: 2px; cursor: pointer;
              border: 1px solid #ccc; background: #f5f5f5; border-radius: 3px; }}
  .tab-btn.active {{ background: #2196F3; color: white; border-color: #1976D2; }}
  .tab-section {{ margin-top: 12px; }}
  figure {{ vertical-align: top; }}
</style>
</head><body>
<h2 style="margin-bottom:8px">Spatial HBCA Xenium — Pipeline Preview (joint_v6)</h2>
<div id="tabs">{buttons}</div>
{sections}
<script>
function show(id) {{
  document.querySelectorAll('.tab-section').forEach(e => e.style.display='none');
  document.querySelectorAll('.tab-btn').forEach(e => e.classList.remove('active'));
  document.getElementById('tab-'+id).style.display='block';
  document.getElementById('btn-'+id).classList.add('active');
}}
show('{first}');
</script>
</body></html>"""


def load_paths(project_root, config_path):
    with open(config_path) as f:
        return {k: os.path.join(project_root, v)
                for k, v in yaml.safe_load(f).items()
                if isinstance(v, str) and not v.startswith("/")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    cfg = load_paths(args.project_root,
                     os.path.join(args.project_root, "config/paths.yaml"))
    previews_dir = cfg["previews_dir"]

    html = build_html(previews_dir)
    out = os.path.join(previews_dir, "summary.html")
    with open(out, "w") as f:
        f.write(html)
    print(f"Wrote {out}", flush=True)


if __name__ == "__main__":
    main()

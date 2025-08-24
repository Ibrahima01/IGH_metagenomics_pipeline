#!/usr/bin/env python3
import os
import subprocess
import shlex
from pathlib import Path
from html import escape

ASSEMBLIES_DIR = Path("Assemblies")
OUT_DIR = Path("Coverage_HTML")
OUT_DIR.mkdir(parents=True, exist_ok=True)
GENOMES_DIR = Path("genomes")

RED = "#8b0000"  # darkred, similar to your screenshot

PAGE_CSS = f"""
body {{ font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace; background: #fafafa; color: {RED}; }}
h1, h2, h3 {{ font-family: system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, Cantarell, Noto Sans, "Helvetica Neue", Arial, "Apple Color Emoji","Segoe UI Emoji"; color: #222; }}
.sample-header {{ margin: 0.5rem 0 1rem; }}
.summary-grid {{
  display: grid; grid-template-columns: 1fr 1fr 1fr; gap: .5rem; margin-bottom: 1rem;
}}
.card {{ background: white; border: 1px solid #eee; border-radius: 10px; padding: .75rem 1rem; }}
pre {{
  background: #fff; border: 1px solid #eee; border-radius: 10px; padding: .75rem 1rem; overflow-x: auto;
  color: {RED};
}}
.details {{ margin-bottom: 1.25rem; }}
.meta {{ color: #444; font-size: .9rem; }}
hr {{ border: none; border-top: 1px solid #eee; margin: 1.25rem 0; }}
a {{ color: #0b6; text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
.badge {{ display:inline-block; background:#eef; color:#225; border:1px solid #dde; border-radius:8px; padding:0 .5rem; margin-left:.5rem; font-size:.8rem; }}
.detail-body {{ display:flex; gap:1rem; align-items:flex-start; }}
.hist {{ flex:1; }}
.blast-table {{ border-collapse:collapse; }}
.blast-table th, .blast-table td {{ border:1px solid #ddd; padding:0.25rem 0.5rem; text-align:right; }}
.blast-table th {{ background:#4caf50; color:white; }}
"""

def run(cmd, check=True, capture=True, text=True):
    """Run shell command and return stdout (strip trailing newline)."""
    res = subprocess.run(shlex.split(cmd), check=check, capture_output=capture, text=text)
    return (res.stdout or "").rstrip("\n")

def ensure_bai(bam: Path):
    bai = bam.with_suffix(bam.suffix + ".bai")
    if not bai.exists():
        run(f"samtools index {bam}")

def first_contig(bam: Path) -> str | None:
    out = run(f"samtools idxstats {bam}")
    for line in out.splitlines():
        fields = line.split("\t")
        if not fields: 
            continue
        if fields[0] != "*":
            return fields[0]
    return None

def mean_coverage_for_contig(bam: Path, contig: str) -> float | None:
    out = run(f"samtools coverage -r {shlex.quote(contig)} {bam}")
    lines = out.splitlines()
    if len(lines) >= 2:
        # samtools coverage table header then one line per region; mean depth is column 7 or 6 depending on version.
        # We want column "meandepth" (usually column 7); coverage% is often column 6. We'll compute mean depth if present, fall back to coverage%.
        header = lines[0].split()
        vals   = lines[1].split()
        def col(name):
            return header.index(name) if name in header else None
        md_i = col("meandepth")
        cov_i = col("coverage")
        if md_i is not None:
            return float(vals[md_i])
        if cov_i is not None:
            return float(vals[cov_i])  # percent, but used only for sorting; still fine
    return None

def hist_for_contig(bam: Path, contig: str) -> str:
    # -A (count all), --hist (ASCII), -w 32 (bin count), -r contig
    return run(f"samtools coverage -A --hist -w 50 -r {shlex.quote(contig)} {bam}")

def blastn_stats(query: Path, subject: Path) -> dict | None:
    fmt = "6 pident length evalue bitscore qcovs scovs"
    try:
        out = run(f"blastn -query {query} -subject {subject} -outfmt '{fmt}'")
    except FileNotFoundError:
        print("⚠️ blastn not found; skipping BLAST stats")
        return None
    except subprocess.CalledProcessError as e:
        print(f"⚠️ blastn failed for {query} vs {subject}: {e}")
        return None
    if not out.strip():
        return None
    vals = out.splitlines()[0].split('\t')
    pident, length, evalue, bitscore, qcovs, scovs = vals
    pident = float(pident)
    length = int(float(length))
    evalue = float(evalue)
    bitscore = float(bitscore)
    qcovs = float(qcovs)
    scovs = float(scovs)
    combined = pident * qcovs / 100
    return {
        "pident": pident,
        "length": length,
        "evalue": evalue,
        "bitscore": bitscore,
        "qcovs": qcovs,
        "scovs": scovs,
        "combined": combined,
    }

def write_sample_page(sample: str, items: list[dict]):
    # items: list of dicts with keys: virus, contig, mean_cov, hist, blast
    # sort by mean_cov desc
    items.sort(key=lambda d: d["mean_cov"], reverse=True)

    html = [f"<!doctype html><html><head><meta charset='utf-8'>",
            f"<title>Coverage histograms — {sample}</title>",
            f"<style>{PAGE_CSS}</style></head><body>",
            f"<h1>Coverage histograms — <code>{escape(sample)}</code></h1>",
            f"<div class='summary-grid'>"
            f"<div class='card'><b>Viruses</b><br>{len(items)}</div>"
            f"<div class='card'><b>Top Mean Coverage</b><br>{(items[0]['mean_cov'] if items else 0):.2f}</div>"
            f"<div class='card'><b>Generated</b><br>samtools coverage --hist</div>"
            f"</div>"]

    for it in items:
        title = f"{it['virus']} <span class='badge'>contig: {escape(it['contig'])}</span> <span class='badge'>mean: {it['mean_cov']:.2f}</span>"
        html.append("<div class='details'>")
        html.append(f"<h3>{title}</h3>")
        html.append("<div class='detail-body'>")
        html.append("<pre class='hist'>")
        html.append(escape(it["hist"]))
        html.append("</pre>")
        blast = it.get("blast")
        if blast:
            html.append("<table class='blast-table'><thead><tr>" +
                        "<th>%IDENTITY</th><th>ALIGN-LENGTH</th><th>E-VALUE</th>" +
                        "<th>BIT-SCORE</th><th>%QUERY</th><th>%SUBJECT</th><th>%COMBINED</th>" +
                        "</tr></thead><tbody>")
            html.append("<tr>" +
                        f"<td>{blast['pident']:.1f}</td>" +
                        f"<td>{blast['length']}</td>" +
                        f"<td>{blast['evalue']:.1e}</td>" +
                        f"<td>{blast['bitscore']:.1f}</td>" +
                        f"<td>{blast['qcovs']:.1f}</td>" +
                        f"<td>{blast['scovs']:.1f}</td>" +
                        f"<td>{blast['combined']:.1f}</td>" +
                        "</tr></tbody></table>")
        html.append("</div>")
        html.append("</div><hr>")
    html.append(f"<p class='meta'>Made with <code>samtools coverage --hist</code> and <code>blastn</code> (when available). Color/style only for readability; numbers are unmodified from tools.</p>")
    html.append("</body></html>")
    out_path = OUT_DIR / f"{sample}.html"
    out_path.write_text("\n".join(html))
    print(f"✅ wrote {out_path}")

def write_index(pages: list[dict]):
    # pages: [{sample, path, viruses, top_cov}]
    pages.sort(key=lambda x: x["top_cov"], reverse=True)
    html = [f"<!doctype html><html><head><meta charset='utf-8'>",
            f"<title>Coverage histograms — index</title>",
            f"<style>{PAGE_CSS}</style></head><body>",
            f"<h1>Coverage histograms — index</h1>",
            f"<p class='meta'>One page per sample. Sorted by highest mean coverage across viruses.</p>",
            "<div class='card'>"]
    html.append("<ul>")
    for p in pages:
        rel = os.path.basename(p["path"])
        html.append(f"<li><a href='{rel}'><code>{escape(p['sample'])}</code></a> — viruses: {p['viruses']} — top mean: {p['top_cov']:.2f}</li>")
    html.append("</ul></div></body></html>")
    (OUT_DIR / "index.html").write_text("\n".join(html))
    print(f"📚 wrote {OUT_DIR/'index.html'}")

def main():
    # Gather per-sample items
    per_sample: dict[str, list[dict]] = {}
    for bam in ASSEMBLIES_DIR.glob("*/*.bam"):
        virus = bam.parent.name
        sample = bam.stem

        try:
            ensure_bai(bam)
            contig = first_contig(bam)
            if not contig:
                continue
            mean_cov = mean_coverage_for_contig(bam, contig)
            if mean_cov is None:
                continue
            hist = hist_for_contig(bam, contig)
            consensus = bam.with_suffix('.fa')
            ref = GENOMES_DIR / f"{virus}.fasta"
            blast = None
            if consensus.exists() and ref.exists():
                blast = blastn_stats(consensus, ref)
        except subprocess.CalledProcessError as e:
            print(f"⚠️ skipping {bam} due to error: {e}")
            continue

        per_sample.setdefault(sample, []).append({
            "virus": virus,
            "contig": contig,
            "mean_cov": float(mean_cov),
            "hist": hist,
            "blast": blast,
        })

    # Write pages
    index_entries = []
    for sample, items in per_sample.items():
        if not items:
            continue
        write_sample_page(sample, items)
        top = max(x["mean_cov"] for x in items)
        index_entries.append({"sample": sample, "path": str(OUT_DIR / f"{sample}.html"), "viruses": len(items), "top_cov": top})

    # Write index
    if index_entries:
        write_index(index_entries)
    else:
        print("No BAMs processed; nothing to write.")

if __name__ == "__main__":
    main()

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import csv


def export_occurrences_csv(hits: List[Dict], out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["motif", "start_1based", "end_1based", "start_0based", "end_0based"]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for h in hits:
            w.writerow({k: h.get(k) for k in fieldnames})

    return out_path


def export_binned_csv(binned: List[Dict], out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["bin_index", "bin_start_1based", "bin_end_1based", "motif", "count"]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in binned:
            w.writerow({k: row.get(k) for k in fieldnames})

    return out_path


def export_summary_csv(stats: Dict[str, float], out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # zapis jednej linijki (nagłówki = klucze)
    keys = list(stats.keys())
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerow(stats)

    return out_path

from typing import Optional, List
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def export_pdf_report(
    out_pdf: str | Path,
    title: str,
    stats: Dict[str, float],
    plot_png: Optional[str | Path] = None,
    extra_lines: Optional[List[str]] = None,
) -> Path:
    out_pdf = Path(out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        title,
        "",
        f"length: {int(stats.get('length', 0))}",
        f"GC_content: {stats.get('GC_content', 0):.4f}",
        f"AT_content: {stats.get('AT_content', 0):.4f}",
        "",
        f"A: {int(stats.get('A', 0))}  C: {int(stats.get('C', 0))}  G: {int(stats.get('G', 0))}  T: {int(stats.get('T', 0))}  N: {int(stats.get('N', 0))}",
    ]
    if extra_lines:
        lines += [""] + extra_lines

    with PdfPages(out_pdf) as pdf:
        
        fig, ax = plt.subplots(figsize=(8.27, 11.69))
        ax.axis("off")
        ax.text(0.05, 0.95, "\n".join(lines), ha="left", va="top", fontsize=12)
        pdf.savefig(fig)
        plt.close(fig)

        if plot_png:
            p = Path(plot_png)
            if p.exists():
                img = plt.imread(str(p))
                fig2, ax2 = plt.subplots(figsize=(11.69, 8.27))
                ax2.axis("off")
                ax2.imshow(img)
                pdf.savefig(fig2)
                plt.close(fig2)

    return out_pdf

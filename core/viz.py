from __future__ import annotations

from pathlib import Path
from typing import List, Dict

import matplotlib
matplotlib.use("Agg")  
import matplotlib.pyplot as plt


def plot_binned_bar(binned: List[Dict], bin_size: int, out_path: str | Path, title: str = "Binned motif counts") -> Path:
    """
    Rysuje wykres słupkowy: liczba trafień motywów w kolejnych binach.
    binned: lista słowników z bin_hits() (bin_index, motif, count).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Zbuduj: bin_index -> {motif -> count}
    by_bin: dict[int, dict[str, int]] = {}
    motifs = set()

    for row in binned:
        b = int(row["bin_index"])
        m = str(row["motif"])
        c = int(row["count"])
        motifs.add(m)
        by_bin.setdefault(b, {})
        by_bin[b][m] = by_bin[b].get(m, 0) + c

    bin_indices = sorted(by_bin.keys())
    motifs = sorted(motifs)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_title(title)
    ax.set_xlabel(f"Bin index (size={bin_size})")
    ax.set_ylabel("Count")

    if not bin_indices:
        ax.text(0.5, 0.5, "No motif hits", ha="center", va="center")
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return out_path

    # Prosty "stacked" bar: dla każdego binu rysujemy kolejne motywy jeden na drugim
    bottoms = [0] * len(bin_indices)
    for m in motifs:
        heights = [by_bin[b].get(m, 0) for b in bin_indices]
        ax.bar(bin_indices, heights, bottom=bottoms, label=m)
        bottoms = [bottoms[i] + heights[i] for i in range(len(bottoms))]

    # Jeżeli binów jest dużo, ukryj etykiety osi X
    if len(bin_indices) > 50:
        ax.set_xticks([])
    else:
        ax.set_xticks(bin_indices)

    ax.legend(title="Motif", ncol=min(3, len(motifs)))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path

def export_interactive_html(hits: List[Dict], out_path: str | Path, title: str = "Motif positions") -> Path:
    """
    Zapisuje interaktywny wykres HTML (plotly): pozycje motywów na osi sekwencji.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        import plotly.express as px
    except Exception:
        raise RuntimeError("Brak plotly. Zainstaluj: pip install plotly")

    if not hits:
        # pusty wykres
        df = {"start_1based": [], "motif": []}
    else:
        df = {
            "start_1based": [h["start_1based"] for h in hits],
            "motif": [h["motif"] for h in hits],
        }

    fig = px.scatter(df, x="start_1based", y="motif", title=title)
    fig.update_layout(xaxis_title="Position (1-based)", yaxis_title="Motif")
    fig.write_html(str(out_path))
    return out_path

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def plot_motif_heatmap_png(
    hits: list[dict],
    motifs: list[str],
    seq_len: int,
    bin_size: int,
    out_path: Path,
    title: str = "Motif heatmap",
) -> Path:
    """
    Heatmapa: wiersze = biny wzdłuż sekwencji, kolumny = motywy, kolor = liczba trafień.
    Zapisuje PNG do out_path.
    """
    motifs = [m.strip().upper() for m in motifs if m.strip()]
    if not motifs:
        raise ValueError("Brak motywów do heatmapy.")

    n_bins = max(1, math.ceil(seq_len / bin_size))
    data = np.zeros((n_bins, len(motifs)), dtype=int)
    motif_to_j = {m: j for j, m in enumerate(motifs)}

    for h in hits:
        m = str(h.get("motif", "")).upper()
        if m not in motif_to_j:
            continue
        start = int(h["start_1based"])
        b = (start - 1) // bin_size
        if 0 <= b < n_bins:
            data[b, motif_to_j[m]] += 1

    fig, ax = plt.subplots(figsize=(max(7, 1.2 * len(motifs)), 5))
    im = ax.imshow(data, aspect="auto", origin="lower")

    ax.set_title(title)
    ax.set_xlabel("Motyw")
    ax.set_ylabel("Bin (wzdłuż sekwencji)")

    ax.set_xticks(range(len(motifs)))
    ax.set_xticklabels(motifs, rotation=45, ha="right")

    # żeby nie było 200 etykiet na osi Y
    step = max(1, n_bins // 12)
    yt = list(range(0, n_bins, step))
    ax.set_yticks(yt)
    ax.set_yticklabels([str(i) for i in yt])

    fig.colorbar(im, ax=ax, label="Liczba wystąpień w binie")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path
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

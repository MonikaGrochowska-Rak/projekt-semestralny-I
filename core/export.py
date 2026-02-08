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
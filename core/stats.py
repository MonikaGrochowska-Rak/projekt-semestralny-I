from __future__ import annotations

from typing import Dict
import math


def basic_stats(sequence: str) -> Dict[str, float]:
    seq = (sequence or "").upper()
    n = len(seq)
    counts = {b: seq.count(b) for b in ["A", "C", "G", "T", "N"]}
    gc = (counts["G"] + counts["C"]) / n if n else 0.0
    at = (counts["A"] + counts["T"]) / n if n else 0.0
    out: Dict[str, float] = {
        "length": n,
        "GC_content": gc,
        "AT_content": at,
        **counts,
    }
    return out


def bin_hits(hits: list[dict], seq_len: int, bin_size: int) -> list[dict]:
    """
    Segmentacja: zliczanie trafień w binach.
    Zwraca listę rekordów:
    bin_index, bin_start_1based, bin_end_1based, motif, count
    """
    if bin_size <= 0:
        raise ValueError("Bin size musi być > 0.")
    if seq_len < 0:
        raise ValueError("seq_len nie może być ujemne.")

    if not hits:
        # zwróć puste — GUI/wykresy poradzą sobie
        return []

    # (bin_index, motif) -> count
    counts: dict[tuple[int, str], int] = {}

    for h in hits:
        start0 = int(h["start_0based"])
        motif = str(h["motif"])
        b = start0 // bin_size
        key = (b, motif)
        counts[key] = counts.get(key, 0) + 1

    out = []
    for (b, motif), cnt in sorted(counts.items(), key=lambda x: (x[0][0], x[0][1])):
        b_start = b * bin_size + 1
        b_end = min((b + 1) * bin_size, seq_len)
        out.append({
            "bin_index": b,
            "bin_start_1based": b_start,
            "bin_end_1based": b_end,
            "motif": motif,
            "count": cnt,
        })

    return out


def num_bins(seq_len: int, bin_size: int) -> int:
    if bin_size <= 0:
        raise ValueError("Bin size musi być > 0.")
    return int(math.ceil(seq_len / bin_size)) if seq_len else 0
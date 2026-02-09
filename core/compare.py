from __future__ import annotations

from typing import Dict, List, Any

from core.motifs import find_multiple_motifs
from core.stats import basic_stats


def motif_counts(sequence: str, motifs: List[str], allow_overlaps: bool = True) -> Dict[str, int]:
    motifs_clean = [m.strip().upper() for m in motifs if m and m.strip()]
    hits = find_multiple_motifs(sequence, motifs_clean, allow_overlaps=allow_overlaps)

    counts: Dict[str, int] = {m: 0 for m in motifs_clean}
    for h in hits:
        m = str(h["motif"])
        counts[m] = counts.get(m, 0) + 1
    return counts


def compare_sequences(seq1: str, seq2: str, motifs: List[str], allow_overlaps: bool = True) -> List[Dict[str, Any]]:
    c1 = motif_counts(seq1, motifs, allow_overlaps=allow_overlaps)
    c2 = motif_counts(seq2, motifs, allow_overlaps=allow_overlaps)

    all_motifs = sorted(set(c1.keys()) | set(c2.keys()))
    rows: List[Dict[str, Any]] = []
    for m in all_motifs:
        a = int(c1.get(m, 0))
        b = int(c2.get(m, 0))
        rows.append({
            "motif": m,
            "count_seq1": a,
            "count_seq2": b,
            "diff_seq2_minus_seq1": b - a,
        })
    return rows


def compare_records(rec1, rec2, motifs: List[str], allow_overlaps: bool = True) -> Dict[str, Any]:
    stats1 = basic_stats(rec1.sequence)
    stats2 = basic_stats(rec2.sequence)
    rows = compare_sequences(rec1.sequence, rec2.sequence, motifs, allow_overlaps=allow_overlaps)

    summary = {
        "id1": rec1.id,
        "id2": rec2.id,
        "len1": int(stats1.get("length", 0)),
        "len2": int(stats2.get("length", 0)),
        "gc1": float(stats1.get("GC_content", 0.0)),
        "gc2": float(stats2.get("GC_content", 0.0)),
    }

    return {
        "summary": summary,
        "rows": rows,
        "stats1": stats1,
        "stats2": stats2,
    }

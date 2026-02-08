from __future__ import annotations

from typing import List, Dict


def find_motif_positions(sequence: str, motif: str, allow_overlaps: bool = True) -> List[Dict]:
    """
    Zwraca listę trafień motywu w sekwencji.
    Pozycje zwracane jako:
    - start_0based / end_0based (0-based)
    - start_1based / end_1based (1-based)
    """
    seq = (sequence or "").upper()
    m = (motif or "").upper().strip()

    if not m:
        return []

    hits = []
    start = 0
    step = 1 if allow_overlaps else len(m)

    while True:
        idx = seq.find(m, start)
        if idx == -1:
            break

        hits.append({
            "motif": m,
            "start_0based": idx,
            "end_0based": idx + len(m) - 1,
            "start_1based": idx + 1,
            "end_1based": idx + len(m),
        })

        start = idx + step

    return hits


def find_multiple_motifs(sequence: str, motifs: List[str], allow_overlaps: bool = True) -> List[Dict]:
    """
    Wyszukuje wiele motywów i zwraca wszystkie trafienia posortowane po pozycji.
    """
    all_hits: List[Dict] = []
    for motif in motifs:
        motif = (motif or "").strip()
        if motif:
            all_hits.extend(find_motif_positions(sequence, motif, allow_overlaps=allow_overlaps))

    all_hits.sort(key=lambda x: (x["start_0based"], x["motif"]))
    return all_hits
from __future__ import annotations

from typing import List, Dict

def reverse_complement(seq: str) -> str:
    """
    Zwraca reverse complement sekwencji DNA (A<->T, C<->G, N->N).
    """
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return (seq or "").translate(table)[::-1]


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


def find_multiple_motifs(
    sequence: str,
    motifs: List[str],
    allow_overlaps: bool = True,
    both_strands: bool = False,
) -> List[Dict]:
    """
    Wyszukuje wiele motywów i zwraca wszystkie trafienia posortowane po pozycji.

    Jeśli both_strands=True:
    - szuka motywu na nici '+' (normalnie)
    - oraz reverse complement motywu na nici '-' (czyli szuka RC w tej samej sekwencji)
    - dodaje pole 'strand' w wynikach: '+' lub '-'
    - dla motywów palindromicznych (motyw == reverse complement) nie dubluje trafień (strand='both')
    """
    all_hits: List[Dict] = []
    seq = (sequence or "").upper()

    for motif in motifs:
        m = (motif or "").strip().upper()
        if not m:
            continue

        # nić +
        plus_hits = find_motif_positions(seq, m, allow_overlaps=allow_overlaps)
        if both_strands:
            for h in plus_hits:
                h["strand"] = "+"
        all_hits.extend(plus_hits)

        if both_strands:
            rc = reverse_complement(m).upper()
            # palindrom: np. CGCG -> RC jest takie samo, nie dublujemy
            if rc == m:
                # oznacz jako 'both' tylko jeśli już były trafienia na '+'
                for h in plus_hits:
                    h["strand"] = "both"
            else:
                minus_hits = find_motif_positions(seq, rc, allow_overlaps=allow_overlaps)
                for h in minus_hits:
                    # motyw wpisany przez użytkownika:
                    h["motif"] = m
                    # opcjonalnie: co realnie znaleźliśmy w sekwencji:
                    h["motif_rc"] = rc
                    h["strand"] = "-"
                all_hits.extend(minus_hits)

    all_hits.sort(key=lambda x: (x["start_0based"], x["motif"], x.get("strand", "+")))
    return all_hits
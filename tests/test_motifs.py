from core.motifs import find_motif_positions, find_multiple_motifs


def test_find_motif_positions_overlaps():
    seq = "ATATATA"
    # "ATA" występuje z overlapami na pozycjach 1,3,5 (1-based): 1-3, 3-5, 5-7
    hits = find_motif_positions(seq, "ATA", allow_overlaps=True)
    starts = [h["start_1based"] for h in hits]
    assert starts == [1, 3, 5]


def test_find_motif_positions_no_overlaps():
    seq = "ATATATA"
    hits = find_motif_positions(seq, "ATA", allow_overlaps=False)
    starts = [h["start_1based"] for h in hits]
    assert starts == [1, 5]


def test_find_multiple_motifs_sorted():
    seq = "ATGxxTATAxxCGCGxxATG".replace("x", "")
    hits = find_multiple_motifs(seq, ["CGCG", "ATG", "TATA"], allow_overlaps=True)

    # powinno być posortowane po pozycji start_0based
    starts0 = [h["start_0based"] for h in hits]
    assert starts0 == sorted(starts0)

    motifs = [h["motif"] for h in hits]
    assert "ATG" in motifs and "TATA" in motifs and "CGCG" in motifs

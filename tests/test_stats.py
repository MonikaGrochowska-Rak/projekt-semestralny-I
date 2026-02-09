from core.stats import basic_stats, bin_hits, num_bins


def test_basic_stats_counts_and_gc():
    seq = "AAGGCCTTNN"
    s = basic_stats(seq)
    assert int(s["length"]) == 10
    assert int(s["A"]) == 2
    assert int(s["G"]) == 2
    assert int(s["C"]) == 2
    assert int(s["T"]) == 2
    assert int(s["N"]) == 2
    # GC = (G+C)/len = 4/10
    assert abs(float(s["GC_content"]) - 0.4) < 1e-9


def test_num_bins():
    assert num_bins(0, 10) == 0
    assert num_bins(1, 10) == 1
    assert num_bins(10, 10) == 1
    assert num_bins(11, 10) == 2


def test_bin_hits_counts_per_bin_and_motif():
    hits = [
        {"motif": "ATG", "start_0based": 0},
        {"motif": "ATG", "start_0based": 1},
        {"motif": "TATA", "start_0based": 9},
        {"motif": "ATG", "start_0based": 10},
    ]
    binned = bin_hits(hits, seq_len=20, bin_size=10)

    # zamieniamy na dict: (bin_index, motif) -> count
    d = {(r["bin_index"], r["motif"]): r["count"] for r in binned}

    assert d[(0, "ATG")] == 2
    assert d[(0, "TATA")] == 1
    assert d[(1, "ATG")] == 1

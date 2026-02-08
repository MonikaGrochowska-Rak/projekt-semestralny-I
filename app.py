from core.io import read_sequence_from_file
from core.motifs import find_multiple_motifs
from core.stats import basic_stats, bin_hits, num_bins
from core.viz import plot_binned_bar

if __name__ == "__main__":
    rec = read_sequence_from_file("data/my_seq1.fasta")
    motifs = ["ATG", "TATA", "CGCG"]
    bin_size = 10

    hits = find_multiple_motifs(rec.sequence, motifs, allow_overlaps=True)
    stats = basic_stats(rec.sequence)
    binned = bin_hits(hits, seq_len=int(stats["length"]), bin_size=bin_size)

    out_png = plot_binned_bar(binned, bin_size=bin_size, out_path="outputs/binned.png", title=f"{rec.id} (len={int(stats['length'])})")

    print("ID:", rec.id)
    print("Length:", int(stats["length"]))
    print("GC:", round(stats["GC_content"], 4))
    print("Total hits:", len(hits))
    print("Num bins:", num_bins(int(stats["length"]), bin_size))
    print("Saved plot:", out_png)

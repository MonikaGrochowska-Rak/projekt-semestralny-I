from core.io import read_sequence_from_file
from core.motifs import find_multiple_motifs
from core.stats import basic_stats, bin_hits, num_bins

if __name__ == "__main__":
    rec = read_sequence_from_file("data/my_seq1.fasta")
    motifs = ["ATG", "TATA", "CGCG"]
    bin_size = 10

    hits = find_multiple_motifs(rec.sequence, motifs, allow_overlaps=True)
    stats = basic_stats(rec.sequence)
    binned = bin_hits(hits, seq_len=stats["length"], bin_size=bin_size)

    print("ID:", rec.id)
    print("Length:", int(stats["length"]))
    print("GC:", round(stats["GC_content"], 4))
    print("Motifs:", motifs)
    print("Total hits:", len(hits))
    print("Num bins:", num_bins(int(stats["length"]), bin_size))
    print("Binned rows:", len(binned))
    print("First 5 binned:")
    for row in binned[:5]:
        print(row)
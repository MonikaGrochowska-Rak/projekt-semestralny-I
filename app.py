from core.io import read_sequence_from_file
from core.motifs import find_multiple_motifs

if __name__ == "__main__":
    rec = read_sequence_from_file("data/my_seq1.fasta")
    motifs = ["ATG", "TATA", "CGCG"]

    hits = find_multiple_motifs(rec.sequence, motifs, allow_overlaps=True)

    print("ID:", rec.id)
    print("Length:", len(rec.sequence))
    print("Motifs:", motifs)
    print("Total hits:", len(hits))

    # pokaż pierwsze 10 trafień
    for h in hits[:10]:
        print(h["motif"], h["start_1based"], h["end_1based"])
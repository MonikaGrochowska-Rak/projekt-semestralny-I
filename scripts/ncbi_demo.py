from core.io import fetch_sequence_from_ncbi
from core.motifs import find_multiple_motifs
from core.stats import basic_stats

if __name__ == "__main__":
    # ZMIEŃ NA SWÓJ EMAIL (nie zapisuj go potem na stałe w kodzie, jeśli nie chcesz)
    EMAIL = "monika.grochowskax@gmail.com"
    ACCESSION = "NC_000913.3"  # E. coli K-12

    rec = fetch_sequence_from_ncbi(ACCESSION, EMAIL)
    stats = basic_stats(rec.sequence)

    motifs = ["ATG", "TATA", "CGCG"]
    hits = find_multiple_motifs(rec.sequence, motifs, allow_overlaps=True)

    print("ID:", rec.id)
    print("Opis:", rec.description[:80])
    print("Length:", int(stats["length"]))
    print("GC:", round(stats["GC_content"], 4))
    print("Total hits:", len(hits))
    print("First 5 hits:", [(h["motif"], h["start_1based"]) for h in hits[:5]])

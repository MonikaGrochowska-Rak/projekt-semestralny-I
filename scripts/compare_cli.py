from __future__ import annotations

import argparse
from pathlib import Path

from core.io import read_sequence_from_file, fetch_sequence_from_ncbi
from core.compare import compare_records
from core.export import export_compare_csv, export_summary_csv


def _parse_motifs(text: str) -> list[str]:
    return [m.strip().upper() for m in (text or "").split(",") if m.strip()]


def main():
    p = argparse.ArgumentParser(description="Compare two DNA sequences by motif counts.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--files", action="store_true", help="Compare two local files")
    src.add_argument("--ncbi", action="store_true", help="Compare two NCBI accessions")

    p.add_argument("--file1", type=str, default="", help="Path to FASTA/TXT (seq1)")
    p.add_argument("--file2", type=str, default="", help="Path to FASTA/TXT (seq2)")

    p.add_argument("--acc1", type=str, default="", help="NCBI accession (seq1)")
    p.add_argument("--acc2", type=str, default="", help="NCBI accession (seq2)")
    p.add_argument("--email", type=str, default="", help="Email required by NCBI")

    p.add_argument("--motifs", type=str, default="ATG,TATA,CGCG", help="Comma-separated motifs")
    args = p.parse_args()

    motifs = _parse_motifs(args.motifs)
    if not motifs:
        raise SystemExit("Podaj motywy, np. --motifs ATG,TATA,CGCG")

    if args.files:
        if not args.file1 or not args.file2:
            raise SystemExit("Dla --files podaj --file1 i --file2.")
        rec1 = read_sequence_from_file(args.file1)
        rec2 = read_sequence_from_file(args.file2)
    else:
        if not args.acc1 or not args.acc2 or not args.email:
            raise SystemExit("Dla --ncbi podaj --acc1, --acc2 oraz --email.")
        rec1 = fetch_sequence_from_ncbi(args.acc1, args.email)
        rec2 = fetch_sequence_from_ncbi(args.acc2, args.email)

    result = compare_records(rec1, rec2, motifs, allow_overlaps=True)
    summary = result["summary"]
    rows = result["rows"]

    out_dir = Path("outputs") / "compare" / f"{summary['id1']}_vs_{summary['id2']}"
    out_dir.mkdir(parents=True, exist_ok=True)

    export_compare_csv(rows, out_dir / "compare.csv")
    export_summary_csv(result["stats1"], out_dir / "stats_seq1.csv")
    export_summary_csv(result["stats2"], out_dir / "stats_seq2.csv")

    print("=== SUMMARY ===")
    print(f"{summary['id1']}  len={summary['len1']}  GC={summary['gc1']:.4f}")
    print(f"{summary['id2']}  len={summary['len2']}  GC={summary['gc2']:.4f}")
    print("\n=== MOTIF COUNTS (seq2 - seq1) ===")
    for r in rows:
        print(f"{r['motif']}: {r['count_seq1']} vs {r['count_seq2']}   diff={r['diff_seq2_minus_seq1']}")

    print(f"\nSaved to: {out_dir}")


if __name__ == "__main__":
    main()

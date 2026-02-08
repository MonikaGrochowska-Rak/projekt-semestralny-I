from __future__ import annotations

from pathlib import Path
from .models import SequenceRecord

_ALLOWED = set("ACGTN")


def _clean_seq(text: str) -> str:
    text = text.upper()
    return "".join(ch for ch in text if ch in _ALLOWED)


def read_sequence_from_file(path: str | Path) -> SequenceRecord:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku: {p}")

    raw = p.read_text(encoding="utf-8", errors="ignore")
    if not raw.strip():
        raise ValueError("Plik jest pusty.")

    # FASTA: ignorujemy linie zaczynające się od '>'
    if raw.lstrip().startswith(">"):
        lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]
        header = lines[0][1:].strip() if lines else ""
        seq_lines = [ln for ln in lines[1:] if not ln.startswith(">")]
        seq = _clean_seq("".join(seq_lines))
        if not seq:
            raise ValueError("Nie udało się odczytać sekwencji z FASTA.")
        rec_id = header.split()[0] if header else p.stem
        desc = header if header else p.name
        return SequenceRecord(id=rec_id, description=desc, sequence=seq)

    # TXT: z całego tekstu wyciągamy A/C/G/T/N
    seq = _clean_seq(raw)
    if not seq:
        raise ValueError("Nie udało się odczytać sekwencji z TXT.")
    return SequenceRecord(id=p.stem, description=p.name, sequence=seq)
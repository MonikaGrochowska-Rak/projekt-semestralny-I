from __future__ import annotations

from pathlib import Path
from .models import SequenceRecord

import os

try:
    from Bio import Entrez, SeqIO
except Exception:
    Entrez = None
    SeqIO = None

import certifi


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

def fetch_sequence_from_ncbi(accession: str, email: str, rettype: str = "gb") -> SequenceRecord:
    """
    Pobiera sekwencję z NCBI (nuccore) po accession.
    rettype: "gb" (GenBank) albo "fasta"
    """
    if Entrez is None or SeqIO is None:
        raise RuntimeError("Brak Biopython. Zainstaluj: pip install biopython")

    accession = accession.strip()
    email = email.strip()

    if not accession:
        raise ValueError("Podaj accession (np. NC_000913.3).")
    if not email:
        raise ValueError("Podaj e-mail (wymagany przez NCBI).")

    # Ustaw certyfikaty (pomaga na Windows przy SSL)
    os.environ["SSL_CERT_FILE"] = certifi.where()
    os.environ["REQUESTS_CA_BUNDLE"] = certifi.where()

    Entrez.email = email
    Entrez.tool = "dna_motif_analyzer_extended"

    # 1) spróbuj GenBank
    try:
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        try:
            record = SeqIO.read(handle, "genbank")
            seq = _clean_seq(str(record.seq))
            if seq:
                return SequenceRecord(id=record.id, description=record.description, sequence=seq)
        finally:
            handle.close()
    except Exception:
        pass

    # 2) fallback FASTA
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
    try:
        record = SeqIO.read(handle, "fasta")
        seq = _clean_seq(str(record.seq))
        if not seq:
            raise ValueError("Pobrano pustą sekwencję z NCBI.")
        return SequenceRecord(id=record.id, description=record.description, sequence=seq)
    finally:
        handle.close()

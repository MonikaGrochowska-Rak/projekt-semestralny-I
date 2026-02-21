# Projekt semestralny – DNA Motif Analyzer (Extended)

Aplikacja w Pythonie do analizy sekwencji DNA pod kątem występowania motywów (np. ATG, TATA, CGCG).

## Funkcje (wariant rozszerzony)
- Wczytanie sekwencji z pliku FASTA/TXT
- Pobieranie sekwencji z NCBI (GenBank)
- Obsługa wielu motywów jednocześnie
- Statystyki sekwencji + segmentacja (binning)
- Wizualizacja rozmieszczenia motywów (wykres binned)
- Interaktywna wizualizacja HTML
- Porównanie 2 sekwencji
- Eksport wyników: CSV + PDF
- Testy: pytest

## Instalacja (Windows)

```bat
py -m venv .venv
.venv\Scripts\activate
py -m pip install -r requirements.txt

## Uruchomienie (GUI)

py app.py

## Testy

py -m pytest

## Przykład NCBI
- Accession: NC_000913.3 (E. coli K-12)
- Motywy: ATG, TATA, CGCG
- Uwaga: NCBI wymaga podania e-mail w zapytaniach (Entrez).

## Przykładowe dane
- `data/my_seq1.fasta`
- `data/my_seq2.fasta`


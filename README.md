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
```bash
py -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
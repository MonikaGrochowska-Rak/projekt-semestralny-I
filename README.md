# Projekt semestralny – DNA Motif Analyzer (Extended)

Aplikacja w Pythonie do analizy sekwencji DNA pod kątem występowania motywów (np. ATG, TATA, CGCG).
Obsługuje wczytanie sekwencji z pliku FASTA/TXT lub pobranie z NCBI (GenBank), statystyki, wizualizacje oraz eksport wyników.

## Funkcje

### Wariant minimalny (MVP)
- [x] Wczytanie sekwencji z pliku FASTA/TXT
- [x] Wyszukiwanie motywu/motywów + pozycje wystąpień
- [x] Statystyki + segmentacja (binning) z użyciem NumPy/Pandas
- [x] Wizualizacja rozmieszczenia motywów (wykres słupkowy/binned)
- [x] GUI: wybór pliku i motywu/motywów
- [x] Eksport wyników do CSV

### Wariant rozszerzony
- [x] Pobieranie sekwencji z NCBI (GenBank)
- [x] Obsługa wielu motywów jednocześnie
- [x] Interaktywna wizualizacja (HTML)
- [x] Porównanie 2 sekwencji
- [x] Eksport: CSV + PDF
- [x] Testy jednostkowe: pytest

## Wymagania
- Windows + Python (zalecane uruchamianie przez `py`)

> Uwaga (Windows): jeśli `python` nie działa w CMD, używaj `py` (Python Launcher).

## Instalacja (Windows)

```bat
py -m venv .venv
.venv\Scripts\activate
py -m pip install -r requirements.txt

## Uruchomienie (GUI)

    python app.py

## Testy

    python -m pytest -q

## Przykład NCBI
- Accession: NC_000913.3 (E. coli K-12)
- Motywy: ATG, TATA, CGCG
- Uwaga: NCBI wymaga podania e-mail w zapytaniach (Entrez).

## Przykładowe dane
- `data/my_seq1.fasta`
- `data/my_seq2.fasta`

## Struktura projektu (skrót)

core/ – logika analizy (motywy, statystyki, porównanie)
tests/ – testy pytest
data/ – przykładowe dane
outputs/ – wygenerowane wyniki
docs/ – dokumentacja/raporty (jeśli dotyczy)


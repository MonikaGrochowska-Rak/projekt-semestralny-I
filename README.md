# Projekt semestralny – DNA Motif Analyzer (Extended)

Aplikacja w Pythonie do analizy sekwencji DNA pod kątem występowania motywów (np. ATG, TATA, CGCG).
Obsługuje wczytanie sekwencji z pliku FASTA/TXT lub pobranie z NCBI (GenBank), statystyki, wizualizacje oraz eksport wyników.

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

## Wymagania
- Windows + Python (zalecane uruchamianie przez `py`)
> Uwaga: jeśli `python` nie działa w CMD, używaj `py` (Python Launcher).

## Instalacja (Windows)
```bat
py -m venv .venv
.venv\Scripts\activate
py -m pip install -r requirements.txt
```

## Uruchomienie (GUI)
```bat
py app.py
```
## Testy
```bat
py -m pytest
```
## Przykład NCBI
```bat
Accession: NC_000913.3 (E. coli K-12)

Motywy: ATG, TATA, CGCG

Uwaga: NCBI wymaga podania e-mail w zapytaniach (Entrez) – ustaw w aplikacji zgodnie z komunikatem w GUI/konfiguracji.
```
## Przykładowe dane
```bat
data/my_seq1.fasta

data/my_seq2.fasta
```
## Struktura projektu (skrót)
```bat
core/ – logika analizy (motywy, statystyki, porównanie)

tests/ – testy pytest

data/ – przykładowe dane

outputs/ – wygenerowane wyniki

docs/ – dokumentacja/raporty (jeśli dotyczy)
```

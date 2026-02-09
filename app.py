from __future__ import annotations

import os
import sys
import threading
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from core.io import read_sequence_from_file, fetch_sequence_from_ncbi
from core.motifs import find_multiple_motifs
from core.stats import basic_stats, bin_hits, num_bins
from core.viz import plot_binned_bar
from core.export import export_occurrences_csv, export_binned_csv, export_summary_csv


def _open_file_default_app(path: Path) -> None:
    try:
        if sys.platform.startswith("win"):
            os.startfile(str(path))  # type: ignore[attr-defined]
        else:
            import webbrowser
            webbrowser.open(path.resolve().as_uri())
    except Exception as e:
        messagebox.showwarning("Nie można otworzyć pliku", str(e))


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DNA Motif Analyzer — Extended (GUI)")
        self.geometry("1050x720")

        # --- state ---
        self._last_out_dir: Path | None = None
        self._last_plot: Path | None = None
        self._last_hits: list[dict] = []
        self._last_binned: list[dict] = []
        self._last_stats: dict = {}
        self._last_rec_id: str = ""

        # --- UI vars ---
        self.source = tk.StringVar(value="file")
        self.file_path = tk.StringVar(value="")
        self.accession = tk.StringVar(value="NC_000913.3")
        self.email = tk.StringVar(value="")
        self.motifs = tk.StringVar(value="ATG, TATA, CGCG")
        self.bin_size = tk.StringVar(value="200")

        self._build_ui()
        self._sync_source()

    def _build_ui(self):
        root = ttk.Frame(self, padding=10)
        root.pack(fill="both", expand=True)

        # ===== Source frame =====
        src = ttk.LabelFrame(root, text="Źródło sekwencji", padding=10)
        src.pack(fill="x")

        ttk.Radiobutton(src, text="Plik (FASTA/TXT)", variable=self.source, value="file",
                        command=self._sync_source).grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(src, text="NCBI (nuccore)", variable=self.source, value="ncbi",
                        command=self._sync_source).grid(row=0, column=1, sticky="w", padx=(12, 0))

        ttk.Label(src, text="Plik:").grid(row=1, column=0, sticky="w", pady=(8, 0))
        self.file_entry = ttk.Entry(src, textvariable=self.file_path, width=70)
        self.file_entry.grid(row=1, column=1, sticky="w", pady=(8, 0))
        ttk.Button(src, text="Wybierz…", command=self._browse).grid(row=1, column=2, padx=8, pady=(8, 0))

        ttk.Label(src, text="Accession:").grid(row=2, column=0, sticky="w", pady=(8, 0))
        self.acc_entry = ttk.Entry(src, textvariable=self.accession, width=25)
        self.acc_entry.grid(row=2, column=1, sticky="w", pady=(8, 0))

        ttk.Label(src, text="E-mail:").grid(row=2, column=2, sticky="e", pady=(8, 0))
        self.email_entry = ttk.Entry(src, textvariable=self.email, width=28)
        self.email_entry.grid(row=2, column=3, sticky="w", padx=(8, 0), pady=(8, 0))

        src.columnconfigure(1, weight=1)

        # ===== Params frame =====
        params = ttk.LabelFrame(root, text="Parametry", padding=10)
        params.pack(fill="x", pady=(10, 0))

        ttk.Label(params, text="Motywy (po przecinku):").grid(row=0, column=0, sticky="w")
        ttk.Entry(params, textvariable=self.motifs, width=60).grid(row=0, column=1, sticky="w", padx=8)

        ttk.Label(params, text="Bin size:").grid(row=0, column=2, sticky="e", padx=(20, 0))
        ttk.Entry(params, textvariable=self.bin_size, width=10).grid(row=0, column=3, sticky="w", padx=8)

        btns = ttk.Frame(params)
        btns.grid(row=1, column=0, columnspan=4, sticky="w", pady=(10, 0))

        self.btn_analyze = ttk.Button(btns, text="Analyze", command=self._analyze_clicked)
        self.btn_analyze.pack(side="left", padx=4)

        self.btn_export = ttk.Button(btns, text="Export CSV", command=self._export_clicked, state="disabled")
        self.btn_export.pack(side="left", padx=4)

        self.btn_open_plot = ttk.Button(btns, text="Open plot", command=self._open_plot_clicked, state="disabled")
        self.btn_open_plot.pack(side="left", padx=4)

        ttk.Button(btns, text="Clear", command=self._clear).pack(side="left", padx=4)

        self.status = tk.StringVar(value="Gotowe.")
        ttk.Label(params, textvariable=self.status).grid(row=2, column=0, columnspan=4, sticky="w", pady=(8, 0))

        # ===== Output frame =====
        out = ttk.Frame(root)
        out.pack(fill="both", expand=True, pady=(10, 0))

        self.summary = tk.Text(out, height=10, wrap="word")
        self.summary.pack(fill="x")
        self.summary.configure(state="disabled")

        cols = ("motif", "start_1based", "end_1based")
        self.tree = ttk.Treeview(out, columns=cols, show="headings", height=14)
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=160, anchor="center")
        self.tree.pack(fill="both", expand=True, pady=(8, 0))

    def _browse(self):
        path = filedialog.askopenfilename(
            title="Wybierz plik FASTA/TXT",
            filetypes=[("FASTA/TXT", "*.fasta *.fa *.fna *.txt"), ("All files", "*.*")]
        )
        if path:
            self.file_path.set(path)

    def _sync_source(self):
        is_file = self.source.get() == "file"
        self.file_entry.configure(state="normal" if is_file else "disabled")
        self.acc_entry.configure(state="disabled" if is_file else "normal")
        self.email_entry.configure(state="disabled" if is_file else "normal")

    def _set_busy(self, busy: bool):
        self.btn_analyze.configure(state="disabled" if busy else "normal")
        self.btn_export.configure(state="disabled" if busy or self._last_out_dir is None else "normal")
        self.btn_open_plot.configure(state="disabled" if busy or self._last_plot is None else "normal")
        self.status.set("Pracuję..." if busy else "Gotowe.")

    def _clear(self):
        self._last_out_dir = None
        self._last_plot = None
        self._last_hits = []
        self._last_binned = []
        self._last_stats = {}
        self._last_rec_id = ""

        for item in self.tree.get_children():
            self.tree.delete(item)

        self.summary.configure(state="normal")
        self.summary.delete("1.0", "end")
        self.summary.configure(state="disabled")

        self.btn_export.configure(state="disabled")
        self.btn_open_plot.configure(state="disabled")
        self.status.set("Wyczyszczono.")

    def _analyze_clicked(self):
        # Walidacja pól
        motifs = [m.strip().upper() for m in self.motifs.get().split(",") if m.strip()]
        if not motifs:
            messagebox.showerror("Błąd", "Podaj co najmniej jeden motyw (np. ATG).")
            return
        try:
            bin_size = int(self.bin_size.get().strip())
            if bin_size <= 0:
                raise ValueError
        except Exception:
            messagebox.showerror("Błąd", "Bin size musi być liczbą całkowitą > 0.")
            return

        # Start wątku (żeby okno się nie wieszało)
        self._set_busy(True)

        def worker():
            try:
                # 1) wczytanie
                if self.source.get() == "file":
                    p = self.file_path.get().strip()
                    if not p:
                        raise ValueError("Wybierz plik FASTA/TXT.")
                    rec = read_sequence_from_file(p)
                else:
                    acc = self.accession.get().strip()
                    em = self.email.get().strip()
                    rec = fetch_sequence_from_ncbi(acc, em)

                # 2) analiza
                hits = find_multiple_motifs(rec.sequence, motifs, allow_overlaps=True)
                stats = basic_stats(rec.sequence)
                binned = bin_hits(hits, seq_len=int(stats["length"]), bin_size=bin_size)

                # 3) zapis domyślny do outputs/gui/<id>_...
                out_dir = Path("outputs") / "gui" / rec.id
                out_dir.mkdir(parents=True, exist_ok=True)

                plot_path = out_dir / "binned.png"
                plot_binned_bar(binned, bin_size=bin_size, out_path=plot_path, title=f"{rec.id} (len={int(stats['length'])})")

                export_occurrences_csv(hits, out_dir / "occurrences.csv")
                export_binned_csv(binned, out_dir / "binned.csv")
                export_summary_csv(stats, out_dir / "summary.csv")

                result = (rec, hits, stats, binned, out_dir, plot_path, motifs, bin_size)
                self.after(0, lambda: self._on_analysis_done(result))
            except Exception as e:
                self.after(0, lambda: self._on_analysis_error(e))

        threading.Thread(target=worker, daemon=True).start()

    def _on_analysis_done(self, result):
        rec, hits, stats, binned, out_dir, plot_path, motifs, bin_size = result

        self._last_out_dir = out_dir
        self._last_plot = plot_path
        self._last_hits = hits
        self._last_stats = stats
        self._last_binned = binned
        self._last_rec_id = rec.id

        # Summary
        counts = {}
        for h in hits:
            counts[h["motif"]] = counts.get(h["motif"], 0) + 1

        lines = [
            f"ID: {rec.id}",
            f"Opis: {rec.description}",
            f"Długość: {int(stats['length'])}",
            f"GC: {stats['GC_content']:.4f}",
            "",
            f"Motywy: {', '.join(motifs)}",
            f"Łącznie trafień: {len(hits)}",
            f"Liczba binów (bin size={bin_size}): {num_bins(int(stats['length']), bin_size)}",
            "",
            "Liczba wystąpień motywów:",
        ]
        if counts:
            for m in sorted(counts.keys()):
                lines.append(f"- {m}: {counts[m]}")
        else:
            lines.append("- brak trafień")

        lines += [
            "",
            f"Zapisano do: {out_dir}",
            f"Wykres: {plot_path}",
        ]

        self.summary.configure(state="normal")
        self.summary.delete("1.0", "end")
        self.summary.insert("end", "\n".join(lines))
        self.summary.configure(state="disabled")

        # Table (nie pokazuj wszystkiego dla ogromnych danych)
        for item in self.tree.get_children():
            self.tree.delete(item)

        max_rows = 2000
        shown = hits[:max_rows]
        for h in shown:
            self.tree.insert("", "end", values=(h["motif"], h["start_1based"], h["end_1based"]))

        if len(hits) > max_rows:
            messagebox.showinfo(
                "Uwaga",
                f"Znaleziono {len(hits)} trafień.\n"
                f"W tabeli pokazuję pierwsze {max_rows}, żeby aplikacja nie zwolniła.\n"
                f"Pełne dane są w CSV: {out_dir / 'occurrences.csv'}"
            )

        self.btn_export.configure(state="normal")
        self.btn_open_plot.configure(state="normal")
        self._set_busy(False)

    def _on_analysis_error(self, e: Exception):
        self._set_busy(False)
        messagebox.showerror("Błąd", str(e))

    def _export_clicked(self):
        if self._last_out_dir is None:
            return

        target = filedialog.askdirectory(title="Wybierz katalog do eksportu CSV")
        if not target:
            return

        target_dir = Path(target) / f"DNA_Motif_Analyzer_{self._last_rec_id}"
        target_dir.mkdir(parents=True, exist_ok=True)

        export_occurrences_csv(self._last_hits, target_dir / "occurrences.csv")
        export_binned_csv(self._last_binned, target_dir / "binned.csv")
        export_summary_csv(self._last_stats, target_dir / "summary.csv")

        messagebox.showinfo("Zapisano", f"CSV zapisane do:\n{target_dir}")

    def _open_plot_clicked(self):
        if self._last_plot is None:
            return
        _open_file_default_app(self._last_plot)


if __name__ == "__main__":
    App().mainloop()

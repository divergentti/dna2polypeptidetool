# Copyright (c) 2025 Jari Hiltunen / GitHub Divergentti
#
# DNA to Polypeptide Message Encoder - GUI Implementation
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QLabel, QVBoxLayout,
    QLineEdit, QPushButton, QTextEdit, QComboBox,
    QMessageBox, QListWidget, QDialog
)
from PySide6.QtGui import QAction, QFont
from PySide6.QtCore import QThread, Signal, QTimer
import sys

# Import core logic from CLI module
from dnaencoder_CLI import CreatePossibilities, IterateFrames, EmbedWords, default_sequence, codon_table_1letter, Config, validate_dna_sequence, validate_words, WordCache

VERSION = "0.0.1 - 04.05.2025"

debug_gui = False


class Worker(QThread):
    """Run long operations in a separate thread."""
    result = Signal(str)
    show_dialog = Signal()

    def __init__(self, func, *args):
        super().__init__()
        self.func = func
        self.args = args

    def run(self):
        """Execute the function and emit the result."""
        try:
            output = self.func(*self.args)
            if output == "SHOW_DIALOG":
                self.show_dialog.emit()
            else:
                self.result.emit(output)
        except Exception as e:
            self.result.emit(f"Error: {str(e)}")
            if debug_gui:
                print("Worker error: %s", e)

class WordListDialog(QDialog):
    """Dialog for displaying and filtering encodable words."""
    def __init__(self, words):
        super().__init__()
        self.setWindowTitle("Encodable Words List")
        self.setMinimumSize(400, 500)

        layout = QVBoxLayout(self)
        self.all_words = words
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Filter words...")
        self.search_input.setToolTip("Enter text to filter the word list.")
        self.word_list = QListWidget()
        self.word_list.addItems(sorted(self.all_words))

        self.search_input.textChanged.connect(self.filter_words)
        layout.addWidget(self.search_input)
        layout.addWidget(self.word_list)

    def filter_words(self, text):
        """Filter the word list based on input text."""
        filtered = [w for w in self.all_words if text.upper() in w.upper()]
        self.word_list.clear()
        self.word_list.addItems(filtered)

class DNAWindow(QMainWindow):
    """Main application window for DNA to polypeptide tool."""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DNA to Polypeptide Encoder " + VERSION)
        self.setMinimumSize(900, 600)

        app_font = QFont("Arial", 11)
        self.setFont(app_font)

        self.possible_words = WordCache().valid_words
        self.amino_letters = set(codon_table_1letter.keys()) - {'*'}

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(12, 12, 12, 12)
        main_layout.setSpacing(14)

        self.operation_combo = QComboBox()
        self.operation_combo.addItems([
            "Select operation",
            "1. Search for a specific word in a DNA sequence",
            "2. Scan DNA for all encodable words",
            "3. Embed 1–3 words into DNA sequence"
        ])
        self.operation_combo.currentIndexChanged.connect(self.update_ui)
        main_layout.addWidget(self.operation_combo)

        self.input_area = QWidget()
        self.input_layout = QVBoxLayout(self.input_area)
        self.input_layout.setSpacing(10)
        main_layout.addWidget(self.input_area)

        self.run_button = QPushButton("Run")
        self.run_button.setEnabled(False)
        self.run_button.setStyleSheet("""
            QPushButton {
                background-color: lightgray;
                border-radius: 6px;
                padding: 6px 12px;
                font-weight: bold;
            }
            QPushButton:enabled {
                background-color: green;
                color: white;
            }
            QPushButton:pressed {
                background-color: darkgreen;
            }
        """)
        self.run_button.clicked.connect(self.execute_action)
        main_layout.addWidget(self.run_button)

        self.clear_button = QPushButton("Clear")
        self.clear_button.clicked.connect(self.clear_inputs)
        main_layout.addWidget(self.clear_button)

        self.progress_label = QLabel("Ready")
        main_layout.addWidget(self.progress_label)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        self.output_text.setLineWrapMode(QTextEdit.NoWrap)
        self.output_text.setMinimumHeight(180)
        main_layout.addWidget(self.output_text)

        self.create_menu()

        self.word_input = QLineEdit()
        self.word_input.setPlaceholderText("Enter word(s)...")
        self.word_input.setFixedHeight(30)
        self.word_input.setToolTip("For embedding, separate multiple words with commas.")

        self.dna_input = QTextEdit()
        self.dna_input.setMinimumHeight(100)
        self.dna_input.setToolTip("Enter a valid DNA sequence using A, G, C, T.")

        self.word_input.textChanged.connect(self.check_run_button_state)
        self.dna_input.textChanged.connect(self.validate_dna_input)
        self.operation_combo.currentIndexChanged.connect(self.check_run_button_state)

    def create_menu(self):
        """Create the About and Help menu."""
        menubar = self.menuBar()
        about_menu = menubar.addMenu("About")

        help_action = QAction("Help", self)
        help_action.triggered.connect(self.show_help)
        about_menu.addAction(help_action)

        basics_action = QAction("DNA Basics", self)
        basics_action.triggered.connect(self.show_basics)
        about_menu.addAction(basics_action)

    def show_help(self):
        """Display help dialog."""
        QMessageBox.information(
            self,
            "Help",
            "DNA to Polypeptide Encoder – Help\n\n"
            "This tool allows you to:\n"
            "1. Search for specific words in a DNA sequence (all six reading frames).\n"
            "2. Scan a DNA sequence for all known encodable English words.\n"
            "3. Embed 1–3 words into a DNA sequence by modifying codons (in any reading frame).\n\n"
            "The default sequence includes the Finnish words:\n"
            "   • KATAINEN (frame 2)\n"
            "   • MINISTERI (frame 3)\n\n"
            "You can find real DNA sequences at:\n"
            "   https://www.ncbi.nlm.nih.gov/\n\n"
            "© 2025 Jari Hiltunen / GitHub Divergentti\n"
            "MIT License"
        )

    def show_basics(self):
        """Display DNA basics dialog."""
        QMessageBox.information(
            self,
            "DNA Basics",
            "DNA Basics\n\n"
            "DNA (deoxyribonucleic acid) is composed of four nucleotides:\n"
            "  • A (adenine)\n"
            "  • T (thymine)\n"
            "  • G (guanine)\n"
            "  • C (cytosine)\n\n"
            "Codons are triplets of nucleotides (e.g., ATG, CCG) that encode amino acids.\n"
            "\nEvery sequence has six possible reading frames:\n"
            "  • Three in the forward strand\n"
            "  • Three in the reverse complement strand\n\n"
            "This tool simulates DNA → amino acid translation using the genetic code.\n\n"
            "Learn more:\n"
            "https://en.wikipedia.org/wiki/Genetic_code"
        )

    def update_ui(self):
        """Update input fields based on selected operation."""
        index = self.operation_combo.currentIndex()
        for i in reversed(range(self.input_layout.count())):
            self.input_layout.itemAt(i).widget().setParent(None)

        if index == 1:
            self.input_layout.addWidget(QLabel("Enter a word to search (or type 'list'):"))
            self.input_layout.addWidget(self.word_input)
            self.input_layout.addWidget(QLabel("DNA Sequence:"))
            self.input_layout.addWidget(self.dna_input)
            self.dna_input.setPlaceholderText(f"Enter DNA sequence (default: {default_sequence[:50]}...)")
        elif index == 2:
            self.input_layout.addWidget(QLabel("Enter a DNA sequence to scan:"))
            self.input_layout.addWidget(self.dna_input)
            self.dna_input.setPlaceholderText(f"Enter DNA sequence to scan (default: {default_sequence[:50]}...)")
        elif index == 3:
            self.input_layout.addWidget(QLabel("Enter 1–3 words to embed, separated by commas:"))
            self.input_layout.addWidget(self.word_input)
            self.input_layout.addWidget(QLabel("DNA sequence to embed into (can be modified):"))
            self.input_layout.addWidget(self.dna_input)
            self.dna_input.setPlaceholderText(f"Enter DNA sequence (default: {default_sequence[:50]}...)")
        self.check_run_button_state()

    def validate_dna_input(self):
        """Validate DNA input and highlight if invalid."""
        dna = self.dna_input.toPlainText().strip()
        try:
            validate_dna_sequence(dna or default_sequence)
            self.dna_input.setStyleSheet("")
        except ValueError:
            self.dna_input.setStyleSheet("border: 1px solid red;")

    def check_run_button_state(self):
        """Enable or disable Run button based on input validity."""
        index = self.operation_combo.currentIndex()
        word_text = self.word_input.text().strip()

        if index == 0:
            self.run_button.setEnabled(False)
            self.run_button.setStyleSheet("background-color: lightgray;")
        elif index in [1, 3] and not word_text:
            self.run_button.setEnabled(False)
            self.run_button.setStyleSheet("background-color: lightgray;")
        else:
            self.run_button.setEnabled(True)
            self.run_button.setStyleSheet("background-color: green; color: white;")

    def clear_inputs(self):
        """Clear all input fields and output."""
        self.word_input.clear()
        self.dna_input.clear()
        self.output_text.clear()
        self.progress_label.setText("Ready")
        self.operation_combo.setCurrentIndex(0)
        self.check_run_button_state()

    def execute_action(self):
        """Start execution in a separate thread."""
        self.run_button.setEnabled(False)
        self.run_button.setStyleSheet("background-color: red; color: white;")
        self.progress_label.setText("Processing...")
        QTimer.singleShot(100, self._start_worker)

    def _start_worker(self):
        """Start the worker thread for the selected operation."""
        index = self.operation_combo.currentIndex()
        words = self.word_input.text().strip()
        dna = self.dna_input.toPlainText().strip() or default_sequence

        if index == 1:
            self.worker = Worker(self._execute_search, words, dna)
        elif index == 2:
            self.worker = Worker(self._execute_scan, dna)
        elif index == 3:
            self.worker = Worker(self._execute_embed, words, dna)
        else:
            self.progress_label.setText("Ready")
            self.check_run_button_state()
            return

        self.worker.result.connect(self._handle_result)
        self.worker.show_dialog.connect(self._show_word_list_dialog)
        self.worker.finished.connect(lambda: self.progress_label.setText("Ready"))
        self.worker.start()

    def _show_word_list_dialog(self):
        """Show the word list dialog."""
        dlg = WordListDialog(self.possible_words)
        dlg.exec()
        self.check_run_button_state()

    def _handle_result(self, result):
        """Display the worker's result and reset UI."""
        self.output_text.setText(result)
        self.check_run_button_state()

    def _execute_search(self, words, dna):
        """Execute word search operation."""
        if words.lower() == "list":
            return "SHOW_DIALOG"
        try:
            finder = IterateFrames(dna)
            results = finder.find_words_in_frames([words.upper()])
            if results:
                return "\n".join([f"🔍 Found '{w}' in {frame}" for w, frame in results])
            return f"'{words}' not found in any reading frame."
        except ValueError as e:
            return f"❌ {e}"

    def _execute_scan(self, dna):
        """Execute scan operation."""
        try:
            min_len = Config.MIN_WORD_LENGTH
            filtered = [w for w in self.possible_words if len(w) >= min_len]
            finder = IterateFrames(dna)
            results = finder.find_words_in_frames(filtered)
            if results:
                return f"✅ Found {len(results)} matches:\n" + "\n".join(f"- {w} ({f})" for w, f in sorted(results))
            return "No known encodable words were found in this DNA sequence."
        except ValueError as e:
            return f"❌ {e}"

    def _execute_embed(self, words, dna):
        """Execute embed operation."""
        word_list = [w.strip().upper() for w in words.split(",") if w.strip()]
        if not (1 <= len(word_list) <= 3):
            return "❌ You must enter 1 to 3 words."
        if not validate_words(word_list, self.amino_letters):
            invalid = [c for w in word_list for c in w if c not in self.amino_letters]
            return (
                f"❌ Some words contain characters not encodable via amino acid 1-letter codes.\n"
                f"Invalid letters: {', '.join(sorted(set(invalid)))}\n\n"
                f"Valid letters: {', '.join(sorted(self.amino_letters))}"
            )

        try:
            embedder = EmbedWords(dna)
            results, timed_out = embedder.try_embed_multiple_words(word_list)
            if not results:
                message = "❌ No suitable embedding found."
                if timed_out:
                    message += f" (timed out after {Config.TIMEOUT} seconds)"
                return message

            lines = [f"✅ Generated {len(results)} sequence(s) with all words embedded:"]
            for idx, r in enumerate(results[:5], 1):
                embedded = ", ".join(f"{w} ({d} frame {f})" for w, d, f in r["embedded"])
                lines.append(f"\n🔢 Option {idx} — embedded: {embedded}")
                lines.append(f"🧬 Modified sequence start:\n{r['seq'][:90]}...")
            return "\n".join(lines)
        except ValueError as e:
            return f"❌ {e}"

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DNAWindow()
    window.show()
    sys.exit(app.exec())

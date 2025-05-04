# Copyright (c) 2025 Jari Hiltunen / GitHub Divergentti
#
# DNA to Polypeptide Message Encoder - CLI
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
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY332333 ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Purpose: Encode messages into polypeptides using amino acid 1-letter codes.
Use: NLTK words matched to valid amino acid sequences.
"""

__all__ = ["Userinterface", "CreatePossibilities", "IterateFrames", "EmbedWords", "default_sequence"]

import nltk
import os
import appdirs
import sys
import json
import shutil
from nltk.corpus import words
from nltk import download
import random
import time
import logging
from tenacity import retry, stop_after_attempt, wait_fixed
from utils import reverse_complement

VERSION = "0.0.1 - 04.05.2025"

# Setup logging
logging.basicConfig(
    filename=os.path.join(appdirs.user_data_dir("DNA To Polypeptide Encoder", "Divergentti"), "app.log"),
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# Debug flags
debug_files = False
debug_generation = False

# --- Configuration ---
class Config:
    """Centralized configuration for the application."""
    MAX_CANDIDATES = 10000
    TIMEOUT = 10
    MIN_WORD_LENGTH = 3
    MIN_DNA_LENGTH = 6
    DATA_DIR = appdirs.user_data_dir("DNA To Polypeptide Encoder", "Divergentti")
    POSSIBLE_WORDS_FILE = os.path.join(DATA_DIR, "possiblewords.json")

# --- Directory setup ---
def get_data_dir():
    """Get or create the application data directory."""
    os.makedirs(Config.DATA_DIR, exist_ok=True)
    return Config.DATA_DIR

POSSIBLE_WORDS_FILE = Config.POSSIBLE_WORDS_FILE

if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    BUNDLE_DIR = sys._MEIPASS
else:
    BUNDLE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_POSSIBLE_WORDS_FILE = os.path.join(BUNDLE_DIR, "possiblewords.json")

# --- Word Cache ---
class WordCache:
    """Singleton for caching valid words."""
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.valid_words = []
            cls._instance.load_words()
        return cls._instance

    def load_words(self):
        """Load valid words from JSON or generate if missing."""
        try:
            with open(POSSIBLE_WORDS_FILE, 'r') as f:
                self.valid_words = json.load(f)
                logging.info("Loaded words from %s", POSSIBLE_WORDS_FILE)
        except (FileNotFoundError, json.JSONDecodeError):
            try:
                with open(DEFAULT_POSSIBLE_WORDS_FILE, 'r') as f:
                    self.valid_words = json.load(f)
                shutil.copyfile(DEFAULT_POSSIBLE_WORDS_FILE, POSSIBLE_WORDS_FILE)
                logging.info("Copied default words to %s", POSSIBLE_WORDS_FILE)
            except Exception as e:
                logging.error("Error loading possible words: %s", e)
                creator = CreatePossibilities()
                creator.generate()
                self.valid_words = creator.valid_words

NLTK_DATA_DIR = os.path.join(Config.DATA_DIR, 'nltk_data')
os.makedirs(NLTK_DATA_DIR, exist_ok=True)
nltk.data.path = [NLTK_DATA_DIR]
os.environ['NLTK_DATA'] = NLTK_DATA_DIR

# --- Constants ---
default_sequence = (
    "TACTTCAAGGCGGAAAAATGATCAACATTAGCACAGAAAGAATTTAATAAAAGCGACGGCGATTAACGAAAACTAATTTAATTTAATTTTTGGGAAAAAA TTTT"
)

default_sequence_contains = "KATAINEN,2,MINISTERI,3"

codon_table_1letter = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    '*': ['TAA', 'TAG', 'TGA']
}

# --- Input Validation ---
def validate_dna_sequence(dna):
    """Validate DNA sequence for correct characters and length."""
    dna_clean = dna.upper().replace(" ", "")
    if len(dna_clean) < Config.MIN_DNA_LENGTH:
        raise ValueError("DNA sequence must be at least 6 nucleotides long.")
    if any(c not in {'A', 'G', 'C', 'T'} for c in dna_clean):
        raise ValueError("Invalid sequence: only A, G, C, T characters are allowed.")
    return dna_clean

def validate_words(words, amino_letters):
    """Validate words for encodability."""
    return all(all(c in amino_letters for c in w) for w in words)

# --- Classes ---
class CreatePossibilities:
    """Generate valid words that can be encoded as amino acid sequences."""
    def __init__(self):
        self.valid_words = []

    @retry(stop=stop_after_attempt(3), wait=wait_fixed(2))
    def download_nltk_words(self):
        """Download NLTK words with retry mechanism."""
        download('words', download_dir=NLTK_DATA_DIR)

    def generate(self):
        """Generate valid words from NLTK corpus."""
        try:
            self.download_nltk_words()
            word_list = set(words.words())
        except Exception as e:
            logging.error("Failed to load NLTK words: %s", e)
            return

        amino_letters = set(codon_table_1letter.keys()) - {'*'}
        for word in word_list:
            upper = word.upper()
            if all(letter in amino_letters for letter in upper):
                self.valid_words.append(upper)

        if debug_generation:
            logging.info("Generated %d valid amino-acid words", len(self.valid_words))

        with open(POSSIBLE_WORDS_FILE, 'w') as f:
            json.dump(self.valid_words, f)
            logging.info("Saved valid words to %s", POSSIBLE_WORDS_FILE)

class IterateFrames:
    """Handle DNA sequence frame iteration and translation to amino acids."""
    def __init__(self, dna_sequence):
        """Initialize with a validated DNA sequence."""
        self.forward = validate_dna_sequence(dna_sequence)
        self.reverse = reverse_complement(self.forward)

    def get_frames(self, seq):
        """Get all three reading frames for a sequence."""
        return [seq[i:] for i in range(3)]

    def translate_frame(self, frame_seq):
        """Translate a DNA frame to an amino acid sequence."""
        codons = [frame_seq[i:i + 3] for i in range(0, len(frame_seq) - 2, 3)]
        aa_seq = ''
        for codon in codons:
            matched = False
            for aa, codon_list in codon_table_1letter.items():
                if codon in codon_list:
                    aa_seq += aa
                    matched = True
                    break
            if not matched:
                aa_seq += '?'
        return aa_seq

    def find_words_in_frames(self, possible_words):
        """Find encodable words in all frames."""
        results = []
        for i, frame_seq in enumerate(self.get_frames(self.forward)):
            aa_seq = self.translate_frame(frame_seq)
            for word in possible_words:
                if word in aa_seq:
                    results.append((word, f"forward frame {i + 1}"))

        for i, frame_seq in enumerate(self.get_frames(self.reverse)):
            aa_seq = self.translate_frame(frame_seq)
            for word in possible_words:
                if word in aa_seq:
                    results.append((word, f"reverse frame {i + 1}"))

        return results

class EmbedWords:
    """Embed words into a DNA sequence by modifying codons."""
    def __init__(self, dna_sequence):
        """Initialize with a validated DNA sequence."""
        self.original_dna = validate_dna_sequence(dna_sequence)
        self.amino_to_codons = codon_table_1letter.copy()

    def try_embed_multiple_words(self, words, timeout=Config.TIMEOUT, max_candidates=Config.MAX_CANDIDATES):
        """Attempt to embed multiple words into the DNA sequence."""
        candidates = []
        start_time = time.time()

        def embed_one(seq, word):
            """Embed a single word into a sequence."""
            results = []
            for direction, strand in [("forward", seq), ("reverse", reverse_complement(seq))]:
                for frame in range(3):
                    aa_seq = ''
                    codons = [strand[i:i + 3] for i in range(frame, len(strand) - 2, 3)]
                    for codon in codons:
                        for aa, codon_list in self.amino_to_codons.items():
                            if codon in codon_list:
                                aa_seq += aa
                                break
                        else:
                            aa_seq += '?'

                    for i in range(len(aa_seq) - len(word) + 1):
                        if time.time() - start_time > timeout:
                            logging.warning("Embedding timed out for word: %s", word)
                            return results
                        new_seq = list(strand)
                        for j, aa in enumerate(word):
                            # Prefer codons that match existing sequence
                            codon = random.choice(self.amino_to_codons[aa])
                            codon_start = frame + (i + j) * 3
                            if codon_start + 3 <= len(new_seq):
                                new_seq[codon_start:codon_start + 3] = list(codon)
                        new_seq_str = ''.join(new_seq)
                        final_seq = new_seq_str if direction == "forward" else reverse_complement(new_seq_str)
                        results.append({
                            'new_seq': final_seq,
                            'direction': direction,
                            'frame': frame + 1
                        })
            return results

        queue = [{'seq': self.original_dna, 'embedded': []}]
        for word in words:
            next_queue = []
            for entry in queue:
                if time.time() - start_time > timeout:
                    break
                res = embed_one(entry['seq'], word)
                for r in res:
                    next_queue.append({
                        'seq': r['new_seq'],
                        'embedded': entry['embedded'] + [(word, r['direction'], r['frame'])]
                    })
            next_queue = next_queue[:max_candidates]
            queue = next_queue[:max_candidates]
            if not queue:
                logging.info("No candidates left after embedding %s", word)
                break

        # Filter results that truly contain all words
        filtered = []
        for entry in queue:
            if time.time() - start_time > timeout:
                break
            seq = entry['seq']
            frame_checker = IterateFrames(seq)
            found = frame_checker.find_words_in_frames(words)
            found_words = {w for w, _ in found}
            if all(w in found_words for w in words):
                filtered.append(entry)

        return filtered

class Userinterface:
    """Command-line interface for the DNA encoder."""
    def __init__(self):
        self.possible_words = WordCache().valid_words
        self.amino_letters = set(codon_table_1letter.keys()) - {'*'}

    def is_encodable(self, word):
        """Check if a word can be encoded with amino acid codes."""
        return all(c in self.amino_letters for c in word)

    def show_menu(self):
        """Display the CLI menu."""
        print("\n=== DNA to Polypeptide Message Encoder ===")
        print("Choose an action:")
        print("1. Search for a specific word in a DNA sequence")
        print("2. Scan a DNA sequence for all known encodable words >2 characters")
        print("3. Try to embed 1–3 words into a DNA sequence (modifying it)")
        print("4. Exit")

    def search_specific_word(self):
        """Search for a specific word in a DNA sequence."""
        print("Enter a word to search (or 'list' to see examples):")
        word = input("> ").strip().upper()

        if word == 'LIST':
            print("Examples of possible words:")
            print(", ".join(self.possible_words[:100]))
            return

        if word in self.possible_words:
            print("✅ The word is in the English dictionary and encodable using amino acids.")
        elif self.is_encodable(word):
            print("ℹ️ The word is NOT in the English dictionary, but it IS encodable with amino acid codes.")
        else:
            print("❌ The word CANNOT be encoded using amino acid 1-letter codes.")
            return

        dna = input("Enter a DNA sequence to search (or leave blank for default):\n> ").strip()
        if not dna:
            dna = default_sequence

        try:
            finder = IterateFrames(dna)
            results = finder.find_words_in_frames([word])
            if results:
                for w, frame in results:
                    print(f"🔍 Found '{w}' in frame {frame}")
            else:
                print(f"'{word}' was not found in any reading frame.")
        except ValueError as e:
            print(f"❌ {e}")

    def find_all_possible_words_in_sequence(self):
        """Scan a DNA sequence for all encodable words."""
        dna = input("Paste the DNA sequence to scan:\n> ").strip()
        if not dna:
            print("No sequence provided.")
            return

        try:
            dna_clean = validate_dna_sequence(dna)
            min_length = Config.MIN_WORD_LENGTH
            filtered_words = [w for w in self.possible_words if len(w) >= min_length]

            print(f"🔎 Scanning for all known English words ≥{min_length} letters that can be encoded...")
            finder = IterateFrames(dna_clean)
            results = finder.find_words_in_frames(filtered_words)

            if results:
                print(f"✅ Found {len(results)} matches:")
                for word, frame in sorted(results):
                    print(f"- {word} (frame {frame})")
            else:
                print("No known encodable words were found in this DNA sequence.")
        except ValueError as e:
            print(f"❌ {e}")

    def embed_words_into_sequence(self):
        """Embed 1-3 words into a DNA sequence."""
        print("Enter 1–3 words to hide, separated by commas:")
        word_input = input("> ").strip()
        words = [w.strip().upper() for w in word_input.split(",") if w.strip()]

        if not (1 <= len(words) <= 3):
            print("❌ You must enter 1 to 3 words.")
            return

        for word in words:
            if not self.is_encodable(word):
                print(f"❌ '{word}' cannot be encoded using amino acid 1-letter codes.")
                return

        dna = input("Enter a DNA sequence to embed into (or leave blank for default):\n> ").strip()
        if not dna:
            dna = default_sequence

        try:
            embedder = EmbedWords(dna)
            results = embedder.try_embed_multiple_words(words)
            if results:
                print(f"✅ Generated {len(results)} sequence(s) with all words embedded:")
                for idx, result in enumerate(results[:5], 1):
                    print(f"\n🔢 Option {idx} — embedded: {result['embedded']}")
                    print(f"🧬 Modified sequence start:\n{result['seq'][:90]}...")
            else:
                print("❌ No suitable embedding found.")
        except ValueError as e:
            print(f"❌ {e}")

    def run(self):
        """Run the CLI interface in a loop until exit."""
        while True:
            self.show_menu()
            choice = input("> ").strip()
            if choice == "1":
                self.search_specific_word()
            elif choice == "2":
                self.find_all_possible_words_in_sequence()
            elif choice == "3":
                self.embed_words_into_sequence()
            elif choice == "4":
                print("Exiting...")
                break
            else:
                print("Invalid choice. Please choose 1, 2, 3, or 4.")

if __name__ == "__main__":
    ui = Userinterface()
    ui.run()
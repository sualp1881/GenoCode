import ast
import argparse
import sys


class GenoEncoder(ast.NodeVisitor):
    """Encodes Python Abstract Syntax Trees (AST) into biological DNA sequences."""

    def __init__(self):
        self.dna_sequence = []

        # Mapping AST Nodes to Genetic Bases
        self.GENETIC_MAP = {
            # A (Adenine): Control Flow & Conditional Logic
            "If": "A",
            "For": "A",
            "While": "A",
            "Try": "A",
            "ExceptHandler": "A",
            "With": "A",
            "Assert": "A",
            # T (Thymine): Structural Definitions & Imports
            "FunctionDef": "T",
            "AsyncFunctionDef": "T",
            "ClassDef": "T",
            "Import": "T",
            "ImportFrom": "T",
            # C (Cytosine): Operations, Assignments & Comparisons
            "Assign": "C",
            "AugAssign": "C",
            "AnnAssign": "C",
            "BinOp": "C",
            "UnaryOp": "C",
            "Compare": "C",
            "BoolOp": "C",
            # G (Guanine): Expressions, Callables & Returns
            "Call": "G",
            "Return": "G",
            "Yield": "G",
        }

    def generic_visit(self, node):
        node_name = type(node).__name__
        if node_name in self.GENETIC_MAP:
            self.dna_sequence.append(self.GENETIC_MAP[node_name])
        super().generic_visit(node)

    def get_dna(self) -> str:
        return "".join(self.dna_sequence)


class GenoAligner:
    """Performs Needleman-Wunsch Dynamic Programming alignment on AST DNA sequences."""

    def __init__(self, match_score=1, mismatch_score=-1, gap_penalty=-2):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty

    def align(self, seq1: str, seq2: str) -> dict:
        n, m = len(seq1), len(seq2)
        score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

        for i in range(n + 1):
            score_matrix[i][0] = i * self.gap_penalty
        for j in range(m + 1):
            score_matrix[0][j] = j * self.gap_penalty

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                char1, char2 = seq1[i - 1], seq2[j - 1]
                match = score_matrix[i - 1][j - 1] + (
                    self.match_score if char1 == char2 else self.mismatch_score
                )
                delete = score_matrix[i - 1][j] + self.gap_penalty
                insert = score_matrix[i][j - 1] + self.gap_penalty

                score_matrix[i][j] = max(match, delete, insert)

        # Backtracking phase
        align1, align2 = [], []
        i, j = n, m

        while i > 0 and j > 0:
            score_curr = score_matrix[i][j]
            score_diag = score_matrix[i - 1][j - 1]
            score_up = score_matrix[i - 1][j]

            char1, char2 = seq1[i - 1], seq2[j - 1]
            step_match = (
                self.match_score if char1 == char2 else self.mismatch_score
            )

            if score_curr == score_diag + step_match:
                align1.append(char1)
                align2.append(char2)
                i -= 1
                j -= 1
            elif score_curr == score_up + self.gap_penalty:
                align1.append(char1)
                align2.append("-")
                i -= 1
            else:
                align1.append("-")
                align2.append(char2)
                j -= 1

        while i > 0:
            align1.append(seq1[i - 1])
            align2.append("-")
            i -= 1

        while j > 0:
            align1.append("-")
            align2.append(seq2[j - 1])
            j -= 1

        aligned_seq1 = "".join(reversed(align1))
        aligned_seq2 = "".join(reversed(align2))

        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
        total_len = len(aligned_seq1)
        similarity = (matches / total_len * 100) if total_len > 0 else 0.0

        return {
            "score": score_matrix[n][m],
            "aligned_seq1": aligned_seq1,
            "aligned_seq2": aligned_seq2,
            "similarity_percent": round(similarity, 2),
            "mutations": self._parse_mutations(aligned_seq1, aligned_seq2),
        }

    def _parse_mutations(self, seq1: str, seq2: str) -> list:
        mutations = []
        for idx, (c1, c2) in enumerate(zip(seq1, seq2)):
            if c1 == c2:
                continue
            if c1 == "-":
                mutations.append(
                    (idx, f"Insertion: Gene '{c2}' added to codebase.")
                )
            elif c2 == "-":
                mutations.append(
                    (idx, f"Deletion: Gene '{c1}' removed from codebase.")
                )
            else:
                mutations.append(
                    (idx, f"Mutation: Gene '{c1}' mutated into '{c2}'.")
                )
        return mutations


def extract_dna(file_path: str) -> str:
    """Reads a Python file and extracts its AST DNA sequence."""
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            code = f.read()
        tree = ast.parse(code)
        encoder = GenoEncoder()
        encoder.visit(tree)
        return encoder.get_dna()
    except Exception as e:
        print(f"[Error] Failed to process {file_path}: {e}")
        sys.exit(1)


def display_report(file1: str, file2: str, result: dict):
    """Prints a structured ASCII report of the genetic alignment."""
    print("\n" + "=" * 60)
    print("               GENOCODE ALIGNMENT REPORT               ")
    print("=" * 60)
    print(f" Source File A : {file1}")
    print(f" Source File B : {file2}")
    print("-" * 60)
    print(f" Genetic Similarity : {result['similarity_percent']}%")
    print(f" Alignment Score    : {result['score']}")
    print("-" * 60)
    print(f" Sequence A : {result['aligned_seq1']}")
    print(f" Sequence B : {result['aligned_seq2']}")
    print("-" * 60)
    print(" Detected Genetic Mutations:")
    if result["mutations"]:
        for idx, desc in result["mutations"]:
            print(f"  [Index {idx:02d}] {desc}")
    else:
        print("  None. The codebases have identical DNA structures.")
    print("=" * 60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="GenoCode: AST-to-DNA Sequence Mapper and Alignment Tool"
    )
    parser.add_argument("file1", help="Path to the original Python file")
    parser.add_argument("file2", help="Path to the modified Python file")

    args = parser.parse_args()

    dna1 = extract_dna(args.file1)
    dna2 = extract_dna(args.file2)

    aligner = GenoAligner()
    result = aligner.align(dna1, dna2)

    display_report(args.file1, args.file2, result)


if __name__ == "__main__":
    main()
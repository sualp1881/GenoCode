# GenoCode 🧬

**GenoCode** is an experimental CLI tool that bridges **Bioinformatics** and **Software Architecture**. It maps Python Abstract Syntax Trees (AST) into biological DNA sequences (A, T, C, G) and aligns them using the **Needleman-Wunsch Dynamic Programming algorithm** to analyze code evolution and architectural similarity.

---

## 💡 How It Works

1. **AST to DNA Encoding**: Ignores variable names, comments, and formatting. Translates core logic into genetic bases:
   - **A (Adenine)**: Control Flow (`If`, `For`, `While`, `Try`)
   - **T (Thymine)**: Structures & Definitions (`FunctionDef`, `ClassDef`, `Import`)
   - **C (Cytosine)**: Operations & Assignments (`Assign`, `BinOp`, `Compare`)
   - **G (Guanine)**: Expressions & Calls (`Call`, `Return`)

2. **Sequence Alignment**: Runs the **Needleman-Wunsch DP Algorithm** ($O(N \times M)$) to calculate genetic similarity percentages and pinpoint exact structural mutations (**Insertions**, **Deletions**, **Mutations**).

---

## 🛠️ Installation & Usage

Clone the repository and run via Python 3:

```bash
git clone https://github.com/sualp1881/GenoCode.git
cd GenoCode

# Compare two Python files
python genocode.py sample_a.py sample_b.py

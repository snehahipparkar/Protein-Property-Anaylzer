 🧬 Protein Property Analyzer (Interactive Version)

**Author:** Sneha Hipparkar  
A simple and interactive **Bioinformatics mini-project** written in Python that calculates the **physicochemical properties** of a protein sequence.  
Users can either provide a **FASTA file** or directly enter a **protein accession number** (e.g., `P0DTC2`, `XP_019944888.2`) to fetch data from **NCBI**.  
All results are automatically saved and opened from a **results** folder on your Desktop.
#  Features
- Accepts user input — either **accession number** or **FASTA file path**  
- Fetches protein sequence automatically from **NCBI**  
- Calculates:
  - Amino acid composition  
  - Molecular weight  
  - Isoelectric point (pI)  
  - Aromaticity  
  - Instability Index  
  - GRAVY (hydrophobicity) score  
- Saves results to a **Desktop/results** folder  
- Automatically opens CSV, PNG, and TXT reports after analysis  
 # Technologies Used
- **Python**
- **Biopython** → Sequence retrieval and analysis  
- **Matplotlib** → Visualization  
- **Pandas** → Data handling and CSV creation  
- **Webbrowser** → Auto-open output files  
#  How to Run
### 1. Clone this repository
git clone https://github.com/<snehahipparkar>/Protein-Property-Analyzer.git
### 2. Install dependencies
pip install biopython matplotlib pandas
### 3. Run the program
python protein_analyzer_interactive.py
### 4. Enter when prompted:
A protein accession number (e.g., P0DTC2) Or a local FASTA file path (e.g., protein.fasta)
Results will be saved and opened from:
C:\Users\<YourName>\Desktop\results\


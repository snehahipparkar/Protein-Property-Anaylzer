# protein_analyzer_interactive.py
# Author: Sneha Hipparkar
# Interactive version — user inputs accession number or FASTA file path

import os
from Bio import SeqIO, Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import pandas as pd
import webbrowser

# ====================================================
# CONFIGURATION
# ====================================================

#  Always save results to Desktop
os.chdir(os.path.join(os.path.expanduser("~"), "Desktop"))

Entrez.email = "snehahipparkar21@gmail.com"  # <-- replace with your actual email

# Create a results folder if it doesn’t exist
os.makedirs("results", exist_ok=True)

# ====================================================
# STEP 1: ASK USER FOR INPUT
# ====================================================
print("\n Welcome to the Protein Property Analyzer!")
print("You can enter either:")
print(" - A protein accession number (e.g., P0DTC2, BAB16373)")
print(" - OR a local FASTA file path (e.g., protein.fasta)\n")

input_data = input("Enter accession number or FASTA filename: ").strip()

# ====================================================
# STEP 2: FETCH PROTEIN SEQUENCE
# ====================================================
print(f"\n🔍 Fetching protein data for: {input_data}")

try:
    if os.path.isfile(input_data):
        # Read from FASTA file
        seq_record = next(SeqIO.parse(input_data, "fasta"))
        protein_seq = str(seq_record.seq)
        protein_id = seq_record.id
        print(f" Loaded sequence from FASTA: {protein_id}")
    else:
        # Fetch from NCBI using accession number
        handle = Entrez.efetch(db="protein", id=input_data, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()
        protein_seq = str(seq_record.seq)
        protein_id = seq_record.id
        print(f" Downloaded sequence for accession: {protein_id}")
except Exception as e:
    print(f" Error fetching sequence: {e}")
    exit()

# ====================================================
# STEP 3: ANALYZE SEQUENCE
# ====================================================
print("\n Analyzing protein sequence...")

analysis = ProteinAnalysis(protein_seq)
aa_comp = analysis.count_amino_acids()
mw = analysis.molecular_weight()
pi = analysis.isoelectric_point()
arom = analysis.aromaticity()
instability = analysis.instability_index()
gravy = analysis.gravy()

# ====================================================
# STEP 4: DISPLAY RESULTS
# ====================================================
print("\n Results:")
print(f"Protein ID: {protein_id}")
print(f"Molecular Weight: {mw:.2f} Da")
print(f"Isoelectric Point (pI): {pi:.2f}")
print(f"Aromaticity: {arom:.3f}")
print(f"Instability Index: {instability:.2f}")
print(f"GRAVY Score: {gravy:.3f}\n")

# ====================================================
# STEP 5: SAVE RESULTS (CSV + Plot)
# ====================================================
df = pd.DataFrame.from_dict(aa_comp, orient="index", columns=["Count"])
df["Percentage"] = df["Count"] / df["Count"].sum() * 100

csv_path = os.path.abspath(f"results/aa_composition_{protein_id}.csv")
plot_path = os.path.abspath(f"results/aa_composition_{protein_id}.png")

df.to_csv(csv_path)
print(f" Saved amino acid composition to: {csv_path}")

# Plot amino acid frequency
plt.figure(figsize=(10,6))
df_sorted = df.sort_values(by="Count", ascending=False)
plt.bar(df_sorted.index, df_sorted["Count"], color="teal", edgecolor="black")
plt.title(f"Amino Acid Composition: {protein_id}")
plt.xlabel("Amino Acids")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig(plot_path)
plt.close()

print(f" Saved bar chart to: {plot_path}")

# ====================================================
# STEP 6: AUTO-OPEN RESULTS
# ====================================================
print("\n Opening results for you...")
webbrowser.open(csv_path)   # opens CSV
webbrowser.open(plot_path)  # opens image

# ====================================================
# STEP 7: GENERATE SUMMARY REPORT (TXT)
# ====================================================
summary_path = os.path.abspath(f"results/summary_{protein_id}.txt")

with open(summary_path, "w") as f:
    f.write("==============================================\n")
    f.write("        PROTEIN PROPERTY ANALYSIS REPORT      \n")
    f.write("==============================================\n\n")
    f.write(f"Protein ID: {protein_id}\n")
    f.write(f"Sequence Length: {len(protein_seq)} amino acids\n\n")
    f.write("----- Computed Properties -----\n")
    f.write(f"Molecular Weight: {mw:.2f} Da\n")
    f.write(f"Isoelectric Point (pI): {pi:.2f}\n")
    f.write(f"Aromaticity: {arom:.3f}\n")
    f.write(f"Instability Index: {instability:.2f}\n")
    f.write(f"GRAVY Score: {gravy:.3f}\n\n")
    f.write("----- Top 5 Most Abundant Amino Acids -----\n")
    for aa, count in sorted(aa_comp.items(), key=lambda x: x[1], reverse=True)[:5]:
        f.write(f"{aa}: {count}\n")
    f.write("\nAll data files and plots are saved in the 'results' folder on Desktop.\n")

print(f" Saved summary report to: {summary_path}")
webbrowser.open(summary_path)  # open summary automatically

print("\n Analysis complete! All results are available in your Desktop 'results' folder.")

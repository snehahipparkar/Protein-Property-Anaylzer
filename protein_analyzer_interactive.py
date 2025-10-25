# ====================================================
# protein_analyzer_complete.py
# Author: Sneha Hipparkar
# Bioinformatics Project — Single & Batch Protein Analysis with UniProt Metadata & Correlation
# ====================================================

import os
import requests
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import webbrowser

# ====================================================
# CONFIGURATION
# ====================================================
Entrez.email = "snehahipparkar21@gmail.com"  # Replace with your email
os.chdir(os.path.join(os.path.expanduser("~"), "Desktop"))  # Save results to Desktop
os.makedirs("results", exist_ok=True)

# ====================================================
# FUNCTION: FETCH METADATA FROM UNIPROT
# ====================================================
def fetch_uniprot_metadata(acc):
    """Fetch metadata (organism, name, sequence) from UniProt if available."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            organism = data["organism"]["scientificName"]
            protein_name = data["proteinDescription"]["recommendedName"]["fullName"]["value"]
            sequence = data["sequence"]["value"]
            return protein_name, organism, sequence
    except Exception:
        pass
    return None, None, None


# ====================================================
# FUNCTION: SINGLE PROTEIN ANALYSIS
# ====================================================
def analyze_single_protein():
    print("\n Welcome to the Protein Property Analyzer!")
    print("Input options:")
    print(" - Protein accession number (e.g., NP_001029195.1)")
    print(" - UniProt ID (e.g., P0DTC2)")
    print(" - Local FASTA file (e.g., myprotein.fasta)\n")

    input_data = input("Enter accession, UniProt ID, or FASTA filename: ").strip()
    protein_seq, protein_id, protein_name, organism = None, input_data, "N/A", "N/A"

    # Try fetching from UniProt first
    print(f"\n Fetching protein data for: {input_data}")
    protein_name, organism, sequence = fetch_uniprot_metadata(input_data)
    if sequence:
        protein_seq = sequence
        print(f" Retrieved from UniProt: {protein_name} ({organism})")
    else:
        try:
            if os.path.isfile(input_data):
                seq_record = next(SeqIO.parse(input_data, "fasta"))
                protein_seq = str(seq_record.seq)
                protein_id = seq_record.id
                print(f"Loaded from FASTA file: {protein_id}")
            else:
                handle = Entrez.efetch(db="protein", id=input_data, rettype="fasta", retmode="text")
                seq_record = SeqIO.read(handle, "fasta")
                handle.close()
                protein_seq = str(seq_record.seq)
                protein_id = seq_record.id
                print(f" Retrieved from NCBI: {protein_id}")
        except Exception as e:
            print(f" Error fetching sequence: {e}")
            return

    # Perform analysis
    analysis = ProteinAnalysis(protein_seq)
    aa_comp = analysis.count_amino_acids()
    mw = analysis.molecular_weight()
    pi = analysis.isoelectric_point()
    arom = analysis.aromaticity()
    instability = analysis.instability_index()
    gravy = analysis.gravy()

    print("\n Results:")
    print(f"Protein ID: {protein_id}")
    print(f"Organism: {organism}")
    print(f"Protein Name: {protein_name}")
    print(f"Length: {len(protein_seq)} aa")
    print(f"Molecular Weight: {mw:.2f} Da")
    print(f"Isoelectric Point (pI): {pi:.2f}")
    print(f"Aromaticity: {arom:.3f}")
    print(f"Instability Index: {instability:.2f}")
    print(f"GRAVY Score: {gravy:.3f}\n")

    # Save amino acid composition
    df = pd.DataFrame.from_dict(aa_comp, orient="index", columns=["Count"])
    df["Percentage"] = df["Count"] / df["Count"].sum() * 100

    csv_path = os.path.abspath(f"results/aa_composition_{protein_id}.csv")
    plot_path = os.path.abspath(f"results/aa_composition_{protein_id}.png")
    df.to_csv(csv_path)

    # Plot amino acid composition
    plt.figure(figsize=(10,6))
    df_sorted = df.sort_values(by="Count", ascending=False)
    plt.bar(df_sorted.index, df_sorted["Count"], color="teal", edgecolor="black")
    plt.title(f"Amino Acid Composition: {protein_id}")
    plt.xlabel("Amino Acids")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()

    # Summary report
    summary_path = os.path.abspath(f"results/summary_{protein_id}.txt")
    with open(summary_path, "w") as f:
        f.write("==============================================\n")
        f.write("        PROTEIN PROPERTY ANALYSIS REPORT      \n")
        f.write("==============================================\n\n")
        f.write(f"Protein ID: {protein_id}\n")
        f.write(f"Protein Name: {protein_name}\n")
        f.write(f"Organism: {organism}\n")
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

    print(f"\n Results saved to Desktop/results/\n")
    webbrowser.open(csv_path)
    webbrowser.open(plot_path)
    webbrowser.open(summary_path)


# ====================================================
# FUNCTION: BATCH ANALYSIS
# ====================================================
def analyze_batch():
    print("\n Batch Mode Activated — Analyze multiple proteins together.")
    file_path = input("Enter filename containing accession/UniProt IDs (e.g., protein_list.txt): ").strip()

    if not os.path.exists(file_path):
        print(f" File '{file_path}' not found.")
        return

    with open(file_path) as f:
        accessions = [line.strip() for line in f if line.strip()]

    print(f" Loaded {len(accessions)} proteins.")
    all_data = []

    for acc in accessions:
        print(f"\n Processing: {acc}")
        protein_name, organism, sequence = fetch_uniprot_metadata(acc)

        if not sequence:
            try:
                handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                sequence = str(record.seq)
                organism, protein_name = "N/A", "N/A"
            except Exception as e:
                print(f" Could not retrieve {acc}: {e}")
                continue

        try:
            analysis = ProteinAnalysis(sequence)
            mw = analysis.molecular_weight()
            pi = analysis.isoelectric_point()
            arom = analysis.aromaticity()
            instability = analysis.instability_index()
            gravy = analysis.gravy()
            length = len(sequence)

            all_data.append({
                "Accession": acc,
                "Protein Name": protein_name,
                "Organism": organism,
                "Length (aa)": length,
                "Molecular Weight (Da)": round(mw, 2),
                "Isoelectric Point (pI)": round(pi, 2),
                "Aromaticity": round(arom, 3),
                "Instability Index": round(instability, 2),
                "GRAVY": round(gravy, 3)
            })
        except Exception as e:
            print(f" Error analyzing {acc}: {e}")
            continue

    df = pd.DataFrame(all_data)
    csv_path = os.path.abspath("results/batch_protein_analysis.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n Results saved to: {csv_path}")

    # Generate correlation plots
    plt.figure(figsize=(7,5))
    plt.scatter(df["Length (aa)"], df["GRAVY"], color="teal", alpha=0.7)
    plt.title("Protein Length vs GRAVY Score")
    plt.xlabel("Protein Length (aa)")
    plt.ylabel("GRAVY Score")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("results/length_vs_gravy.png")
    plt.close()

    plt.figure(figsize=(7,5))
    plt.scatter(df["Length (aa)"], df["Isoelectric Point (pI)"], color="orange", alpha=0.7)
    plt.title("Protein Length vs Isoelectric Point (pI)")
    plt.xlabel("Protein Length (aa)")
    plt.ylabel("Isoelectric Point (pI)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("results/length_vs_pI.png")
    plt.close()

    summary_path = os.path.abspath("results/batch_summary.txt")
    with open(summary_path, "w") as f:
        f.write("==============================================\n")
        f.write("     BATCH PROTEIN PROPERTY ANALYSIS REPORT    \n")
        f.write("==============================================\n\n")
        f.write(f"Total Proteins Analyzed: {len(df)}\n\n")
        f.write("Average Property Values:\n")
        f.write(df[["Molecular Weight (Da)", "Isoelectric Point (pI)", "GRAVY", "Instability Index"]].mean().to_string())
        f.write("\n\nPlots:\n - length_vs_gravy.png\n - length_vs_pI.png\n")

    print(" Generated correlation plots and summary report.")
    webbrowser.open(csv_path)
    webbrowser.open("results/length_vs_gravy.png")
    webbrowser.open("results/length_vs_pI.png")
    webbrowser.open(summary_path)


# ====================================================
# MAIN MENU
# ====================================================
print("""
====================================================
  Protein Property Analyzer — Complete Edition
====================================================
Choose a mode:
1️⃣  Single Protein Analysis (NCBI / UniProt / FASTA)
2️⃣  Batch Protein Analysis (with UniProt & plots)
====================================================
""")

choice = input("Enter your choice (1 or 2): ").strip()
if choice == "1":
    analyze_single_protein()
elif choice == "2":
    analyze_batch()
else:
    print(" Invalid choice. Please enter 1 or 2.")

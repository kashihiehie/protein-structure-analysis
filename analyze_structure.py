#!/usr/bin/env python3
"""
Protein Structure Analysis Tool
Analyzes AlphaFold predicted protein structures (.pdb or .cif)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    from Bio.PDB import PDBParser, MMCIFParser
    print("✅ BioPython installed")
except ImportError:
    print("❌ BioPython not installed. Run: python3 -m pip install biopython")
    exit(1)

def parse_structure(file_path):
    """Automatically parse either PDB or mmCIF (.cif) structure files"""
    file_path = Path(file_path)
    if file_path.suffix.lower() == '.cif':
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
        
    structure = parser.get_structure('protein', str(file_path))
    return structure

def calculate_plddt_score(structure):
    """Extract pLDDT confidence scores per residue (using CA atoms)"""
    scores = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id('CA'):
                    b_factor = residue['CA'].get_bfactor()
                    scores.append(b_factor)
    return np.array(scores)

def calculate_radius_of_gyration(structure):
    """Calculate protein radius of gyration using C-alpha backbone atoms"""
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id('CA'):
                    coords.append(residue['CA'].coord)
    
    if not coords:
        return 0.0

    coords = np.array(coords)
    centroid = coords.mean(axis=0)
    distances_sq = ((coords - centroid) ** 2).sum(axis=1)
    rg = np.sqrt(distances_sq.mean())
    return rg

def main():
    data_dir = Path("data")
    
    # Search for .cif or .pdb files inside data and any subfolders
    structure_files = list(data_dir.rglob("*.cif")) + list(data_dir.rglob("*.pdb"))
    
    print("=" * 60)
    print("PROTEIN STRUCTURE ANALYSIS")
    print("=" * 60)
    
    if not structure_files:
        print(f"❌ Error: No structure files (.cif or .pdb) found in '{data_dir.resolve()}'")
        return
        
    structure_file = structure_files[0]
    print(f"\n📂 Loading structure: {structure_file.relative_to(data_dir)}")
    
    # Parse structure
    structure = parse_structure(structure_file)
    
    # Get per-residue pLDDT scores
    plddt_scores = calculate_plddt_score(structure)
    
    if len(plddt_scores) == 0:
        print("❌ Error: No C-alpha (CA) atoms found in the structure file.")
        return

    print(f"\n📊 Structure Statistics:")
    print(f"   Total residues: {len(plddt_scores)}")
    print(f"   Mean pLDDT: {plddt_scores.mean():.2f}")
    print(f"   Min pLDDT:  {plddt_scores.min():.2f}")
    print(f"   Max pLDDT:  {plddt_scores.max():.2f}")
    
    # Radius of gyration
    rg = calculate_radius_of_gyration(structure)
    print(f"   Radius of Gyration: {rg:.2f} Å")
    
    # Create output directory
    Path("figures").mkdir(exist_ok=True)
    
    # Plot pLDDT scores per residue
    print("\n📈 Creating visualizations...")
    plt.figure(figsize=(12, 6))
    plt.plot(range(1, len(plddt_scores) + 1), plddt_scores, linewidth=2, color='steelblue')
    plt.axhline(70, color='orange', linestyle='--', alpha=0.7, label='Low Confidence (70)')
    plt.axhline(90, color='green', linestyle='--', alpha=0.7, label='High Confidence (90)')
    plt.xlabel('Residue Position')
    plt.ylabel('pLDDT Confidence Score')
    plt.title('AlphaFold Confidence Scores Along Protein Sequence')
    plt.ylim([0, 100])
    plt.legend(loc='lower right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('figures/plddt_scores.png', dpi=300)
    plt.close()
    print("   ✅ Saved: figures/plddt_scores.png")
    
    # Histogram distribution
    plt.figure(figsize=(10, 6))
    plt.hist(plddt_scores, bins=20, color='steelblue', edgecolor='black')
    plt.xlabel('pLDDT Score')
    plt.ylabel('Frequency (Residue Count)')
    plt.title('Distribution of AlphaFold Confidence Scores')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig('figures/plddt_histogram.png', dpi=300)
    plt.close()
    print("   ✅ Saved: figures/plddt_histogram.png")
    
    print("\n✅ Analysis complete!")

if __name__ == "__main__":
    main()
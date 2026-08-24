cat > README.md << 'EOF'
# Protein Structure Analysis

Using AI to predict and analyze protein 3D structures with AlphaFold.

## What is This Project?

This project uses AlphaFold2 (via ColabFold) to predict how proteins fold into their 3D shapes, then analyzes the predicted structures to understand their properties and biological significance.

## What is AlphaFold?

AlphaFold is an artificial intelligence system created by Google DeepMind that can predict protein 3D structures from amino acid sequences with remarkable accuracy. This breakthrough solved a 50-year-old problem in structural biology and won the 2024 Nobel Prize in Chemistry.

## How It Works

1. **Input**: Amino acid sequence (just letters like ATGC...)
2. **AlphaFold AI**: Analyzes the sequence and predicts 3D structure
3. **Analysis**: Calculate confidence scores and structural properties
4. **Visualization**: Generate charts showing the results

## Key Concepts

### pLDDT Score
pLDDT (predicted Local Distance Difference Test) measures AlphaFold's confidence (0-100):
- **90-100**: Very high confidence (reliable structure)
- **70-90**: Confident (generally reliable)
- **50-70**: Low confidence (may be unreliable)
- **<50**: Very unreliable (should not be trusted)

### Radius of Gyration
Measures how compact the protein is. Lower values = tighter, more compact structure.

## Files

- `data/protein_structure.pdb` - Predicted 3D structure
- `figures/plddt_scores.png` - Confidence scores along sequence
- `figures/plddt_histogram.png` - Score distribution
- `analyze_structure.py` - Analysis script

## Usage

### Step 1: Get Structure from AlphaFold
Visit: https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb
- Paste your protein sequence
- Run AlphaFold
- Download the protein_structure.pdb file
- Save to `data/` folder

### Step 2: Analyze Structure
```bash
python analyze_structure.py
```

This generates charts and statistics in the `figures/` folder.

## Biological Significance

Understanding protein structures enables:
- **Drug Discovery**: Design drugs that bind to specific proteins
- **Disease Research**: Understand mutations and their effects
- **Synthetic Biology**: Design new proteins with custom functions
- **Enzyme Engineering**: Improve enzyme efficiency

## Author
Kashika Vaish | https://github.com/kashihiehie

## References
- AlphaFold Paper: https://www.nature.com/articles/s41586-020-2828-1
- ColabFold: https://github.com/sokrypton/ColabFold
- PDB Format: https://www.rcsb.org/docs/programmatic-access/file-download-services
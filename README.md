# FoldTrust

**Ensemble-first RNA structure reliability reports for disease-relevant windows.**

FoldTrust is a personal portfolio project that demonstrates RNA structure analysis expertise by turning ViennaRNA's minimum free energy (MFE) structure and base-pair probabilities into a **structure report card**. It identifies which stems you can trust, which look crisp in the MFE but are floppy in the thermodynamic ensemble, and recommends next steps.

[![CI](https://github.com/kanekalla/foldtrust/actions/workflows/ci.yml/badge.svg)](https://github.com/kanekalla/foldtrust/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## The Problem

RNA structure prediction tools typically show only the **minimum free energy (MFE) structure**—a single cartoon of the most stable fold. But that cartoon can be misleading:

- Some stems that appear stable in the MFE have **low pair probabilities** in the Boltzmann ensemble
- Alternative structures may compete
- Local regions can be floppy even when the global MFE looks crisp

For **antisense oligonucleotide (ASO) design**, **splice-switching therapeutics**, and **antiviral targeting**, understanding ensemble reliability is essential.

## The Solution

FoldTrust runs ViennaRNA's partition function to compute base-pair probabilities across the full ensemble, then:

1. Parses MFE stems
2. Calculates mean pair probability for each stem
3. Flags stems as **FIRM** (P ≥ 0.85), **SOFT** (0.5 ≤ P < 0.85), or **FLOPPY** (P < 0.5)
4. Generates an HTML + Markdown report with heatmap and verdict

**Key principle:** _Ask for base-pair probabilities (or the ensemble), not only the MFE._

## Disease Case Gallery

FoldTrust ships with five curated disease-relevant RNA cases demonstrating real-world therapeutic and diagnostic relevance:

### 1. SARS-CoV-2 Frameshift Element (sars2-fse)
**Disease:** COVID-19  
**Gene:** SARS-CoV-2 ORF1ab frameshifting pseudoknot  
**Teaching point:** Competing structures; MFE alone insufficient for antiviral design  
**References:** [DOI:10.1126/science.abc3546](https://doi.org/10.1126/science.abc3546)

### 2. SMN2 ISS-N1 (smn2-iss-n1)
**Disease:** Spinal Muscular Atrophy (SMA)  
**Gene:** SMN2 intron 7 (Spinraza/nusinersen target)  
**Teaching point:** ASO site selection requires local pairing probability assessment  
**References:** [DOI:10.1056/NEJMoa1702752](https://doi.org/10.1056/NEJMoa1702752)

### 3. CFTR 5' UTR (cftr-5utr)
**Disease:** Cystic Fibrosis  
**Gene:** CFTR 5' untranslated region  
**Teaching point:** Lung disease; UTR structure affects translation/splicing hypotheses  
**References:** [DOI:10.1152/physrev.00025.2017](https://doi.org/10.1152/physrev.00025.2017)

### 4. MAPT Exon 10 (mapt-e10)
**Disease:** Frontotemporal dementia, tauopathies  
**Gene:** MAPT (tau) exon 10 splice site  
**Teaching point:** Disease splicing + structure-aware oligo design  
**References:** [DOI:10.1093/hmg/10.10.1029](https://doi.org/10.1093/hmg/10.10.1029)

### 5. HCV IRES Domain II (hcv-ires-dii)
**Disease:** Hepatitis C  
**Gene:** HCV IRES ribosome entry site  
**Teaching point:** Structured viral RNA; ensemble flexibility near functional loops  
**References:** [DOI:10.1006/jmbi.1999.2918](https://doi.org/10.1006/jmbi.1999.2918)

## Installation

### Prerequisites

**macOS:**
```bash
brew install viennarna
```

**Ubuntu/Debian:**
```bash
sudo apt-get install vienna-rna
```

### Install FoldTrust

```bash
git clone https://github.com/kanekalla/foldtrust.git
cd foldtrust
pip install -e .
```

Or install with dev dependencies for testing:
```bash
pip install -e ".[dev]"
```

## Usage

### Quick Demo

Run all five disease cases:

```bash
foldtrust demo
```

Reports will be generated in `examples/out/<case-name>/`.

### Single Sequence Report

```bash
foldtrust report data/cases/sars2-fse/sequence.fa -o output/sars2-fse
```

Open `output/sars2-fse/report.html` in your browser.

### Batch Processing

Process all cases in a directory:

```bash
foldtrust batch data/cases -o output/batch
```

### Custom Sequence

Create a FASTA file `my_rna.fa`:

```
>my_rna_sequence
GGGAAACCCUUUUGGGGAAAACCCCUUUUU
```

Then run:

```bash
foldtrust report my_rna.fa -o output/my_rna
```

## Output

Each report includes:

- **Sequence metadata:** Length, GC%, MFE energy
- **Stem table:** Position, length, mean pair probability, flag (FIRM/SOFT/FLOPPY)
- **Heatmap:** Upper-triangle base-pair probability matrix
- **Verdict:** TRUST / REDESIGN / NEED PROBING
- **Methods:** Energy model documentation and limitations

### Example: Identifying Floppy Stems

The **SMN2 ISS-N1** case demonstrates a region where the MFE structure contains stems that have low ensemble support (floppy stems), making them potentially better targets for antisense oligonucleotides—accessibility is higher when the structure is less stable.

## Scientific Honesty: Limitations

Thermodynamic models operate under simplified assumptions:
- **37°C in 1M NaCl** (not cellular conditions)
- No co-transcriptional folding effects
- No RNA-binding proteins or modifications
- No consideration of longer-range tertiary contacts

**Experimental validation** via SHAPE, DMS probing, or functional assays remains essential. FoldTrust is a hypothesis generation tool, not a substitute for wet-lab data.

## Testing

Run the test suite:

```bash
pytest
```

Tests include:
- Unit tests for stem parsing and classification
- Smoke test with a simple hairpin
- Case structure validation

**Note:** Tests that invoke ViennaRNA are skipped if `RNAfold` is not installed.

## Development

### Linting

```bash
ruff check src/ tests/
black --check src/ tests/
```

### Auto-format

```bash
black src/ tests/
```

## Requirements

- **Python:** 3.11+
- **ViennaRNA:** 2.x (command-line `RNAfold` with partition function support)

No GPU or deep learning training required.

## Repository Structure

```
foldtrust/
├── src/foldtrust/         # Python package
│   ├── cli.py             # Typer CLI
│   ├── core.py            # Main processing logic
│   ├── vienna.py          # ViennaRNA integration
│   ├── report.py          # HTML/Markdown generation
│   ├── visualization.py   # Heatmap plotting
│   └── utils.py           # Utilities
├── data/cases/            # Five curated disease cases
│   ├── sars2-fse/
│   ├── smn2-iss-n1/
│   ├── cftr-5utr/
│   ├── mapt-e10/
│   └── hcv-ires-dii/
├── tests/                 # Pytest suite
├── examples/out/          # Demo output (generated)
├── pyproject.toml         # Package metadata
└── README.md
```

## Roadmap / Future Work

- Web app deployment for interactive reports
- Docker container with ViennaRNA pre-installed
- Integration with SHAPE/DMS reactivity data for constrained folding
- Full-transcriptome scanning mode
- ASO design optimizer with accessibility scoring
- Support for ML-based structure predictors (e.g., RNAfold + AlphaFold-RNA ensemble)

## Author

**Kishore Anekalla**  
Personal portfolio project  
GitHub: [@kanekalla](https://github.com/kanekalla)

## License

MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use FoldTrust in academic work, please cite the underlying methods:

- **ViennaRNA Package 2.0:** Lorenz et al., _Algorithms Mol Biol_ 6:26 (2011). [DOI:10.1186/1748-7188-6-26](https://doi.org/10.1186/1748-7188-6-26)
- **Turner energy parameters:** Mathews et al., _PNAS_ 101(19):7287 (2004). [DOI:10.1073/pnas.0401799101](https://doi.org/10.1073/pnas.0401799101)

And reference each disease case by its published DOI (see case `meta.yaml` files).

---

**FoldTrust:** Because the MFE cartoon is only the beginning of the story.

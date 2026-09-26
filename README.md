# FoldTrust

**Ensemble-first RNA structure reliability reports for disease-relevant windows.**

FoldTrust turns ViennaRNA's minimum free energy (MFE) structure and base-pair probabilities into a **structure report card**. It identifies which stems you can trust, which look crisp in the MFE but are floppy in the thermodynamic ensemble, and recommends next steps.

Scientific notes (question, theory, results, interpretation): [NOTES.md](NOTES.md)

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

FoldTrust ships with five curated disease-relevant RNA cases with verified NCBI provenance (see `docs/benchmark/layer0_cases.md`):

### 1. SARS-CoV-2 Frameshift Element (sars2-fse)
**Coordinates:** NC_045512.2:13462-13542(+), 81 nt  
**Disease:** COVID-19  
**Biology:** Programmed −1 ribosomal frameshift element (slippery site + 3-stem pseudoknot) controls ORF1a/ORF1ab polyprotein ratio  
**Teaching point:** Competing structures; MFE alone insufficient for antiviral design  
**References:** [Kelly et al. 2020](https://doi.org/10.1074/jbc.AC120.013449), [Bhatt et al. 2021](https://doi.org/10.1126/science.abf3546)

### 2. SMN2 ISS-N1 (smn2-iss-n1)
**Coordinates:** NG_008728.1:31999-32152(+), 154 nt  
**Disease:** Spinal Muscular Atrophy (SMA)  
**Biology:** SMN2 exon 7 (classic numbering) + intron 7 ISS-N1 (nusinersen/Spinraza antisense target at intron 7 +10..+27)  
**Teaching point:** ASO site selection requires local pairing probability assessment, not just MFE cartoons  
**References:** [Singh et al. 2006](https://doi.org/10.1128/MCB.26.4.1333-1346.2006), [Finkel et al. 2017](https://doi.org/10.1056/NEJMoa1702752)

### 3. CFTR 5' UTR and Start Codon Region (cftr-5utr)
**Coordinates:** NM_000492.4:1-200(+), 200 nt  
**Disease:** Cystic Fibrosis  
**Biology:** Complete 5' UTR (1-70) + first 130 nt of CDS; structure around start codon affects translation efficiency  
**Teaching point:** UTR structure influences ribosome loading; start codon accessibility matters for expression  
**References:** [Riordan et al. 1989](https://doi.org/10.1126/science.2475911), [Zielenski & Tsui 1991](https://doi.org/10.1016/0888-7543(91)90503-7)

### 4. MAPT Exon 10 5'ss Stem-Loop (mapt-e10)
**Coordinates:** NG_007398.2:120818-121000(+), 183 nt  
**Disease:** Frontotemporal dementia with parkinsonism (FTDP-17), tauopathies  
**Biology:** Exon 10 (classic numbering, encodes tau R2 repeat) + 5' splice site stem-loop; mutations that destabilize the stem increase exon 10 inclusion (pathogenic 4R/3R ratio shift)  
**Teaching point:** cis-regulatory RNA structures control alternative splicing; ensemble analysis reveals mutation effects  
**References:** [Grover et al. 1999](https://doi.org/10.1074/jbc.274.21.15134), [Varani et al. 1999](https://doi.org/10.1073/pnas.96.14.8229), [Hutton et al. 1998](https://doi.org/10.1038/31508)

### 5. HCV IRES Domain II (hcv-ires-dii)
**Coordinates:** AF009606.1:44-118(+), 75 nt  
**Disease:** Hepatitis C  
**Biology:** HCV genotype 1a internal ribosome entry site domain II (5' NTR nt 44-118); structured element for cap-independent translation  
**Teaching point:** 2 FIRM / 3 SOFT, verdict NEED PROBING (case_metrics.csv); every predicted stem is lost when genomic flanks are added (Layer 5), so this window needs probing data before use  
**References:** [Honda et al. 1999](https://doi.org/10.1128/JVI.73.2.1165-1174.1999), [Lukavsky et al. 2003](https://doi.org/10.1038/nsb1004)

**Provenance note:** Original SMN2, CFTR, and MAPT sequences were invented (no NCBI source). HCV was the wrong region (domain III at 148-415, not domain II at 44-118). All were rebuilt from verified NCBI coordinates. See `docs/benchmark/layer0_cases.md` for audit details.

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

### Benchmark Results

**Six benchmark layers** (see [`BENCHMARK.md`](BENCHMARK.md)) test FIRM/SOFT/FLOPPY tiers using unit tests, reference structures, calibration analysis, SHAPE agreement, and robustness checks.

```bash
# Run all layers (~3 min, 16GB machine)
foldtrust benchmark all --output benchmarks/outputs
```

**Key findings:**
- **Layer 2:** MEA F1 = 0.563 on 600 ArchiveII/Rfam/bpRNA structures
- **Layer 3:** FIRM tier PPV = 0.674; AUROC = 0.889; ECE = 0.066
- **Layer 4:** SHAPE correlation ρ 0.29–0.55 for 4 of 5 datasets (Pyle 0.16 n.s.); SARS-CoV-2 FSE
- **Layer 5:** FIRM retention 0.95 at 25/30 °C and 1.00 at 42 °C; 0.60–0.80 with alternate parameter sets

Higher tiers have higher reference agreement. Drug-discovery applications (hypothesis-level): ASO design (SMN2, MAPT), small-molecule targeting (FSE, HCV IRES). See [`docs/benchmark/impact.md`](docs/benchmark/impact.md).

**Limitations:** Thermodynamic model only; FLOPPY tier uninformative (PPV 0.15); parameter/context sensitive.

Full methods, data sources, and reproduce instructions: [`BENCHMARK.md`](BENCHMARK.md)
- References with DOIs for all methods and datasets

## Output

Each report includes:

- **Sequence metadata:** Length, GC%, MFE energy
- **Stem table:** Position, length, mean pair probability, flag (FIRM/SOFT/FLOPPY)
- **Heatmap:** Upper-triangle base-pair probability matrix
- **Verdict:** TRUST / REDESIGN / NEED PROBING
- **Methods:** Energy model documentation and limitations

### Example: Identifying Floppy Stems

The SMN2 ISS-N1 case illustrates a region where most MFE stems have low ensemble support (6 SOFT / 1 FLOPPY of 9; case_metrics.csv), which could make it more accessible to antisense oligonucleotides (hypothesis).

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

- **Benchmark:** Six-layer benchmark completed (see `BENCHMARK.md` and `NOTES.md` § 9 for full results): unit tests, reference structure accuracy (600 structures from ArchiveII/Rfam/bpRNA), calibration analysis, SHAPE agreement (SARS-CoV-2 FSE, 5 datasets), and robustness checks (temperature, parameter sets, window context).
- Web app deployment for interactive reports
- Docker container with ViennaRNA pre-installed
- Integration with SHAPE/DMS reactivity data for constrained folding
- Full-transcriptome scanning mode
- ASO design optimizer with accessibility scoring
- Support for ML-based structure predictors (e.g., RNAfold + AlphaFold-RNA ensemble)

## Author

**Kishore Anekalla**  
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

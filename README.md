# 🧬 nf-tumor-clonality

<div align="center">

![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-brightgreen?style=for-the-badge&logo=nextflow)
![Docker](https://img.shields.io/badge/Docker-Multi--Conda--Env-2496ED?style=for-the-badge&logo=docker)
![Python](https://img.shields.io/badge/Python-Data%20Engineering-3776AB?style=for-the-badge&logo=python)
![R](https://img.shields.io/badge/R-ClonEvol%20Reporting-276DC3?style=for-the-badge&logo=r)
#![License](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg?style=for-the-badge)

**A scalable Nextflow DSL2 pipeline for tumor subclonality analysis.**  
Automated data engineering via Python · Multi-Conda-environment Dockerization · Reproducible phylogenetic reporting.

</div>

---

## 📖 Overview

This repository hosts a **production-ready Nextflow DSL2 pipeline** that automates the full workflow of tumor subclonality analysis — from raw somatic variant calls and copy number data through to phylogenetic tree reconstruction and final reporting.

The pipeline integrates legacy computational tools (PhyloWGS, ClonEvol) with modern workflow orchestration and containerization, directly addressing the **reproducibility crisis** common in cancer phylogenetics. A single Docker image (`divinwilson/nf_clonality:v1`) hosts multiple isolated Conda environments, with each pipeline stage injecting the correct environment's `PATH` at runtime — a practical engineering solution to deep dependency conflicts between Python 2, Python 3, and R toolchains.

---

## ✨ Key Features

| Feature | Description |
|---|---|
| 🔁 **Nextflow DSL2** | Modular workflow with clean module imports and channel-based data flow |
| 🐍 **Python Data Engineering** | Custom scripts for SSM prioritisation, input formatting, and result parsing |
| 🐳 **Single-Image, Multi-Env Docker** | One image hosts `pybedtools`, `PhyloWGS`, and `ClonEvol` Conda envs — activated per-process |
| 📊 **Reproducible Reporting** | Automated clonality trees and mutation assignment tables via R (ClonEvol) |
| 🔖 **Gene-of-Interest Tracking** | Priority SSM file generation flags user-defined genes throughout the workflow |
| ♻️ **Resume Capability** | Nextflow caching — interrupted runs resume without re-computation |

---

## 🏗️ Repository Structure

```
nf-tumor-clonality/
├── main.nf                                           # Pipeline entry point (DSL2)
├── nextflow.config                                   # Params, Docker & process config
│
├── modules/
│   ├── preparation.nf                                # GENERATE_PRIORITY_SSM, PREPARE_PHYLOWGS_INPUTS
│   ├── phylowgs.nf                                   # RUN_MULTIEVOLVE
│   ├── formatting.nf                                 # WRITE_RESULTS, FORMAT_RESULTS
│   └── clonevol.nf                                   # MAKE_CLONEVOL_INFILES, CLONEVOL
│
├── bin/
│   ├── python_generate_priority_SSM_file.py          # Priority SSM generation & gene annotation
│   ├── python_generate_out_treefiles.py              # PhyloWGS output parsing & formatting
│   ├── python_generate_ClonEvol_infile.py            # ClonEvol input file generation
│   └── Rscript_clonevol.R                            # Clonal tree visualisation (ClonEvol)
│
└── utils/
    └── gencode.v45.chr_patch_hapl_scaff.basic.annotation.bed   # Gene annotation reference
```

---

## 🔬 Pipeline Workflow

The pipeline is composed of **7 sequential processes** organised across 4 modules:

```
  VCF files + CNV files + Gene-of-Interest list
                    │
                    ▼
    ┌───────────────────────────────┐
    │   GENERATE_PRIORITY_SSM       │  [preparation.nf | pybedtools env]
    │   Python                      │  Intersects VCFs with gene annotations;
    │                               │  generates a priority SSM TSV with gene info
    └───────────────┬───────────────┘
                    │
                    ▼
    ┌───────────────────────────────┐
    │   PREPARE_PHYLOWGS_INPUTS     │  [preparation.nf | PhyloWGS env]
    │   Python                      │  Reformats SSM + CNV data into
    │                               │  PhyloWGS-compatible input files
    └───────────────┬───────────────┘
                    │
                    ▼
    ┌───────────────────────────────┐
    │   RUN_MULTIEVOLVE             │  [phylowgs.nf | PhyloWGS env]
    │   PhyloWGS                    │  MCMC-based subclone inference;
    │                               │  samples posterior tree distribution
    └───────────────┬───────────────┘
                    │
                    ▼
    ┌───────────────────────────────┐
    │   WRITE_RESULTS               │  [formatting.nf | PhyloWGS env]
    │   PhyloWGS                    │  Serialises tree summaries and
    │                               │  mutation assignments to JSON / TSV
    └───────────────┬───────────────┘
                    │
                    ▼
    ┌───────────────────────────────┐
    │   FORMAT_RESULTS              │  [formatting.nf | pybedtools env]
    │   Python                      │  Merges tree outputs with priority SSM
    │                               │  and gene info into analysis-ready TSVs
    └───────────────┬───────────────┘
                    │
                    ▼
    ┌───────────────────────────────┐
    │   MAKE_CLONEVOL_INFILES       │  [clonevol.nf | pybedtools env]
    │   Python                      │  Transforms formatted results into
    │                               │  ClonEvol-compatible input files
    └───────────────┬───────────────┘
                    │
                    ▼
    ┌───────────────────────────────┐
    │   CLONEVOL                    │  [clonevol.nf | ClonEvol env]
    │   R                           │  Infers clonal ordering; generates
    │                               │  publication-ready phylogenetic tree PDFs
    └───────────────────────────────┘
```

---

## 📦 Docker & Environment Design

A key engineering decision is the use of a **single Docker image hosting three Conda environments**, with the correct environment activated per-process via `containerOptions` in `nextflow.config`:

```
divinwilson/nf_clonality:v1
│
├── /opt/conda/envs/pybedtools/    ← GENERATE_PRIORITY_SSM
│                                     FORMAT_RESULTS
│                                     MAKE_CLONEVOL_INFILES
│
├── /opt/conda/envs/PhyloWGS/      ← PREPARE_PHYLOWGS_INPUTS
│                                     RUN_MULTIEVOLVE
│                                     WRITE_RESULTS
│
└── /opt/conda/envs/ClonEvol/      ← CLONEVOL
```

This design **eliminates Python 2 / Python 3 / R dependency conflicts** without requiring separate image builds and registry pushes per tool. The correct `PATH` is injected at the process level inside `nextflow.config` using `containerOptions = '--env PATH=...'`.

---

## 🚀 Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [Docker](https://www.docker.com/)
- Java 11+

### Clone the Repository

```bash
git clone https://github.com/<your-username>/nf-tumor-clonality.git
cd nf-tumor-clonality
```

### Run the Pipeline

```bash
nextflow run main.nf \
  --sample_id   "SAMPLE_001" \
  --projectDir  "/path/to/your/data"
```

The pipeline resolves `cnv_files`, `vcf_files`, and `goi_file` automatically from `projectDir` using the naming conventions in `nextflow.config`. Output is written to `results/<sample_id>/` by default.

---

## ⚙️ Configuration

All parameters are defined in `nextflow.config` and can be overridden at runtime via `--param value`:

```groovy
params {
    sample_id     = "sample1"
    projectDir    = "/path/to/data"
    outdir        = "results/${params.sample_id}"

    // Input file glob patterns (auto-resolved from projectDir)
    cnv_files     = "${params.projectDir}/CN_DATA/${params.sample_id}*__cnvs.txt"
    vcf_files     = "${params.projectDir}/VCF_DATA/${params.sample_id}*_reformated.vcf"
    vcf_types     = "pcawg_consensus"
    goi_file      = "${params.projectDir}/gene_of_interest.txt"

    // PhyloWGS MCMC subsampling
    subsample_num = 5000
}
```

| Parameter | Description | Default |
|---|---|---|
| `--sample_id` | Sample identifier used in all output filenames | `sample1` |
| `--projectDir` | Root data directory containing `CN_DATA/` and `VCF_DATA/` | Required |
| `--outdir` | Output directory | `results/<sample_id>` |
| `--vcf_types` | VCF parser type passed to PhyloWGS input script | `pcawg_consensus` |
| `--subsample_num` | Mutations to subsample for MCMC | `5000` |

---

## 📊 Output Files

All results are written to `results/<sample_id>/`:

```
results/
└── sample1/
    ├── sample1_priority_ssm_file.tsv           # Priority SSM (pre-annotation)
    ├── sample1_priority_ssm_with_geneinfo.tsv  # Priority SSM with gene annotations
    ├── sample1_ssm_data.txt                    # PhyloWGS SSM input
    ├── sample1_cnv_data.txt                    # PhyloWGS CNV input
    ├── sample1_params.json                     # PhyloWGS run parameters
    ├── sample1.muts.json(.gz)                  # PhyloWGS mutation posteriors
    ├── sample1.treesumm.json(.gz)              # PhyloWGS tree summaries
    ├── sample1.mutass.zip                      # Mutation assignments (compressed)
    ├── trees.zip                               # Full posterior tree archive
    ├── sample1__treesum.tsv                    # Parsed tree summary table
    ├── sample1__mutass_with_geneinfo.tsv       # Mutation assignments + gene info
    ├── sample1_tree.tsv                        # Final tree structure table
    ├── sample1_variant.tsv                     # Variant-level results
    ├── sample1_clonality_tree.pdf              # ClonEvol phylogenetic tree
    ├── sample1_clonality_tree_clean.pdf        # ClonEvol tree (publication-ready)
    └── mutass/                                 # Per-tree mutation assignment files
```

---

## 🛠️ Engineering Highlights

This project was built to demonstrate applied bioinformatics engineering skills in a production context:

**Workflow Orchestration** — Nextflow DSL2 with modular process imports and channel transformation operators (`map`, `combine`, `join`, `collect`) to route data cleanly across seven processes. Intermediate channels are explicitly constructed to join outputs from non-adjacent processes before passing to downstream modules.

**Pragmatic Containerisation** — Rather than maintaining and pushing multiple Docker images, a single image hosts three Conda environments. The correct environment is activated per process group by injecting `PATH` via `containerOptions` in `nextflow.config` — a clean solution to real-world legacy tool dependency conflicts.

**Data Engineering as Glue** — Four custom scripts in `bin/` bridge tools that were never designed to interoperate, handling gene annotation intersections (`pybedtools`), format conversion, result parsing, and structured output generation across the Python 2 → Python 3 → R boundary.

**Reproducibility by Design** — All parameters are externalised to config, Docker is enforced for every process, and Nextflow's content-addressed caching makes every run resumable, traceable, and portable across execution environments.

**Separation of Concerns** — Business logic lives in `bin/`, process definitions in `modules/`, and orchestration in `main.nf`, keeping each layer independently readable and maintainable.

---

## 📚 References

- [PhyloWGS](https://github.com/morrislab/phylowgs) — Deshwar et al., *Genome Biology* (2015)
- [ClonEvol](https://github.com/hdng/clonevol) — Nguyen et al., *Annals of Oncology* (2017)
- [Nextflow](https://www.nextflow.io/) — Di Tommaso et al., *Nature Biotechnology* (2017)

---

## 📄 License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

---

<div align="center">
Built to demonstrate production-grade bioinformatics engineering · <a href="https://github.com/<your-username>">@your-username</a>
</div>

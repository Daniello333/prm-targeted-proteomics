# Targeted proteomics PRM planning — MSc research methods

[![R helper checks](https://github.com/Daniello333/prm-targeted-proteomics/actions/workflows/test.yml/badge.svg)](https://github.com/Daniello333/prm-targeted-proteomics/actions/workflows/test.yml)

> **Publication status:** Runnable R helper-method excerpt from a five-stage MSc research workflow. This repository is **not yet the entire peptide-selection / PRM scheduling pipeline**.

This project develops a computational workflow to propose peptide precursor lists for **parallel reaction monitoring (PRM)**, based on DIA-NN evidence, peptide sequence constraints, precursor properties and chromatographic retention-time information.

## Research workflow (broader MSc project)

1. In-silico digestion of the target-protein FASTA panel and peptide uniqueness checks
2. Filtering and ranking by identification evidence and fragment-ion scores from DIA-NN results
3. Selection of Top 2 / Top 3 precursor candidates per protein
4. Retention-window crowding assessment and alternative-peptide selection
5. iRT-based 60-minute to 30-minute gradient retention-time transfer and list extension

## Code currently published

| File | Purpose |
| --- | --- |
| [`selected_prm_methods.R`](selected_prm_methods.R) | Trypsin/Lys-C sequence cleavage (as implemented in the thesis), peptide neutral masses with fixed cysteine carbamidomethylation, precursor m/z calculations and iRT transfer |
| [`test_selected_prm_methods.R`](test_selected_prm_methods.R) | Base-R checks for these standalone helpers |

### Quick start

Requires R (base installation).

```bash
Rscript test_selected_prm_methods.R
```

## Reported study results

| Candidate coverage | Initial candidate lists | iRT-extended lists |
| --- | ---: | ---: |
| Proteins with selected precursors | 596 | 609 |
| Top 2 peptides | 1,087 | 1,129 |
| Top 3 peptides | 1,493 | 1,563 |

The exported iRT-transfer validation reports **MAE 0.333 min** and **RMSE 0.424 min**. These results are from the complete thesis analyses, **not recomputed by the standalone excerpt**.

## Reproducibility and data note

- Full reproduction requires the five complete research stages, original target-panel FASTA, DIA-NN outputs and other experiment-specific metadata.
- Large experimental inputs and institution-owned data should not be published without permission.
- The generated candidate lists still require instrument-method review and experimental validation; they are **not** validated PRM assays.
- The larger project is described in the [MSc research case study](https://github.com/Daniello333/Daniello333/blob/main/projects/prm-proteomics.md).

**Author:** Daniel Drabik · Jagiellonian University · [Portfolio](https://github.com/Daniello333)

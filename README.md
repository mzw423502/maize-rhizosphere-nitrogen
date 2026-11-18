# Nitrogen responses in the maize rhizosphere across field trials

Research compendium for *A shared community response to nitrogen in the maize rhizosphere spans field trials despite a local reversal and few consistent family-level indicators*. [Version 1.0.0 and downloadable simulation datasets](https://github.com/mzw423502/maize-rhizosphere-nitrogen/releases/tag/v1.0.0).

| Material | Location and contents |
|---|---|
| Field observations | [Empirical data](15_FINAL_AUDIT/PHASE5_D3_RANK_TRANSFER_20261010_V1/inputs/) contain integer family counts, labels and analysis metadata for 3,306 A, 2,123 B and 93 C profiles. Match the first column of the count table to the public run identifiers. B genotype-class files are subsets of B. |
| Annual model inputs | [Model data](15_FINAL_AUDIT/) include original R count matrices, annual selection and metadata, experiment B metadata and 214 family-level A/B synthesis records. |
| Paper tables | [Tables S1–S15](tables/Supplementary_Tables_S1-S15.zip). S4–S7 cover family evidence, target scores, source training and annual annotations. S8–S15 cover simulation estimators, decisions, intervals, calibration and decision utility. |
| Figures | [Manuscript-order figure inputs and original vector PDFs](figures/) for Figures 1–5 and S1–S3. |
| Simulation data | [26 scenario archives](https://github.com/mzw423502/maize-rhizosphere-nitrogen/releases/tag/v1.0.0), c00 through c25. Each holds 200 datasets, their native estimator outputs and their derived analysis records. Each dataset has three count tables (A1, A2 and B). |
| Source reads | Original public NCBI deposits for [experiment A](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA772177/), [experiment B](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1092665/), and [experiment C](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA669400/). Use the run identifiers above to identify the analyzed subset. |

## Use the compendium

1. Download the source archive for v1.0.0 or clone the repository and check out the tag. It contains tables, figure inputs, empirical observations and code.
2. To work with simulation records, extract the scenario ZIPs of interest **at the checkout root**. Match the same filename stem across data, estimates and analysis, for example c00_r0000. Each dataset records its seed, cell, replicate, effect configuration and integer counts. Native outputs and derived results include the original input hash. The data/formal manifest in the benchmark directory lists all inputs and their hashes.
3. Run **python verify.py** (Python 3.11 or later) from the checkout to check the repository files and empirical sample counts. With all scenario ZIPs downloaded into one directory, run **python verify.py PATH_TO_DOWNLOADS** to check the archived files.
4. To redraw Figures 1–5 and S1–S3, install matplotlib, numpy and pillow, then run **python plotting/render_figures.py PATH_TO_OUTPUT**. Rendering uses the packaged numbers. Fonts and panel spacing can vary across graphics environments. The original vector PDFs provide the reference artwork.
5. [Reproduction instructions](REPRODUCTION.md) describe the field and simulation code and its dependencies. The release covers the three field experiments and the synthetic benchmark.

## Citation

Cite the [versioned release](https://github.com/mzw423502/maize-rhizosphere-nitrogen/releases/tag/v1.0.0) and the original data producers for their BioProject deposits. CITATION.cff gives the authors in manuscript order. Zimeng Niu and Zhongwu Ma contributed equally.

To download just scenarios 0 and 25, use **python download_simulations.py downloads 0 25**. Omit the scenario numbers to fetch all 26. The downloader checks the SHA-256 of every file.

# Scientific code and reproduction

## Field study

Four treatment-allocation units in each of A and B and eight units in C determine the field comparisons. The B genotype-class transfers reuse B observations. Repeated years, settings and subsamples do not add independently allocated treatments.

The methods.py file in 16_MANUSCRIPT/ALTERNATIVE_ANGLE_20261011/code implements rank and covariance calculations. Its companion compute.py reads the packaged counts, run metadata, frozen source-training eligibility and saved family components. For this distribution its input checks use evidence/PUBLIC_INPUT_HASHES.json. To execute the four components, run the following command from the checkout root with a fresh output location. The analysis requires Python, NumPy and SciPy.

    python 16_MANUSCRIPT/ALTERNATIVE_ANGLE_20261011/code/compute.py formal PATH_TO_NEW_OUTPUT

Tables S1–S7 report the corresponding quantities. The implementations directory preserves the R model routines for ANCOM/ALDEx, the D2 calculator, transfer functions and the two-experiment synthesis. Some original fitting workers depend on their project controller and package runtime. The distributed saved-count and saved-component analyses above operate without that controller. The field data in this distribution start at family counts; the original reads are available from NCBI.

## End-to-end benchmark

The recorded environment used Python 3.12.10, NumPy 2.5.1, SciPy 1.18.0 and R 4.4.3 with modeest 2.5.0, foreach 1.5.2, matrixStats 1.5.0 and jsonlite 2.0.0. LinDA was sourced from the pinned MicrobiomeStat 1.4 revision. The fetch_pinned_linda.py helper downloads that revision and checks the recorded upstream SHA-256 hashes. The upstream project's terms apply to its files.

To regenerate the benchmark, copy code/ and upstream/ from METHODOLOGY_E2E_ITERATION_20261011 into a **new empty benchmark directory** and run these commands there after installing dependencies. To inspect the existing results, use the downloaded archives directly. Pilot and formal are distinct executions.

    python code/fetch_pinned_linda.py
    python code/generate.py formal
    Rscript code/fit_linda.R formal
    python code/analyze.py formal
    python code/audit.py formal

data/formal contains the simulated counts, estimates/formal the paired inputs and native LinDA output, and analysis/formal the derived annual and synthesis records. The sample-independence comparator, paired-unit CLR and LinDA targets are distinguished in the source. Within each scenario the first 100 datasets are used for simulation calibration and the second 100 for evaluation of decision utility.

## Saved-component decision utility

Extract all 26 scenario ZIPs at the repository root and keep the benchmark and utility directories as siblings. Inside METHODOLOGY_DECISION_VALUE_FINAL_20261011 run **python code/utility.py formal**, followed by **python code/audit_utility.py formal**. Use a fresh checkout for each calculation. The benchmark's provenance folder records hashes of the original saved components. Table S12 reports the decision-utility summaries, S13 the calibration estimates and S14 the dataset and taxon decisions.

This release draws on existing outputs. Integrity checks and figure rendering for public distribution are separate from the original scientific computations. The original native fitting process produced all 5,200 expected output records; its operating-system exit code was unavailable. The saved numerical outputs were independently audited in the original study.

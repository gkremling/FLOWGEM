# FLOWGEM: FLOW-based GEeneration for Missing data

This repository contains the official implementation of the paper:

**[Generative Modeling under Non-Monotonic MAR Missingness via Approximate Wasserstein Gradient Flows](https://arxiv.org/abs/2604.04567)**
Gitte Kremling, Jeffrey Näf, Johannes Lederer

*FLOWGEM is a principled and practical method for generating complete datasets from data with values Missing at Random (MAR). It uses a Wasserstein Gradient Flow with an approximated velocity field to minimize the expected KL divergence.*


[comment]: <> (Installation, requirements.txt)

## Usage
The main experiments are run via two notebooks:
- `main_realdata.ipynb` --- runs and evaluates different imputation/generation methods on real datasets
- `main_simulations.ipynb` --- runs and evaluates methods for the simulation study described in the paper

The R script `main_realdata.R` prepares the real datasets by amputing them and inserting missing values (typically does not need to be re-run).

## Repository Structure
`Bayesian/` contains the implementation of the Bayesian approach from [Chérief-Abdellatif and Näf (2026)](https://arxiv.org/abs/2603.23449)
`MIRI-Imputation/` contains the implementations of MIRI, MissDiff and NewImp from [Yu et al. (2025)](https://arxiv.org/abs/2505.11749)
`data/` contains the real datasets
`results/` contains the saved results (included for real datasets; excluded for simulations due to file size)

## Citation

If you find this work useful, please cite:
```bibtex
@misc{kremling2026generativemodelingnonmonotonicmar,
      title={Generative Modeling under Non-Monotonic MAR Missingness via Approximate Wasserstein Gradient Flows}, 
      author={Gitte Kremling and Jeffrey Näf and Johannes Lederer},
      year={2026},
      url={https://arxiv.org/abs/2604.04567}, 
}
```

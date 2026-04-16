# FLOWGEM: FLOW-based GEeneration for Missing data

This repository contains the official implementation of the paper:

**[Generative Modeling under Non-Monotonic MAR Missingness via Approximate Wasserstein Gradient Flows](https://arxiv.org/abs/2604.04567)**
Gitte Kremling, Jeffrey Näf, Johannes Lederer

*FLOWGEM is a principled and practical method for generating complete datasets from data with values Missing at Random (MAR). It uses a Wasserstein Gradient Flow with an approximated velocity field to minimize the expected KL divergence.*


[comment]: <> (Installation, requirements.txt)

## Usage

The method is implemented in `flowgem.py` and can be applied to a dataset with missing values as follows:
```python
from flowgem import sample_flowgem
Xnew = sample_flowgem(X0, X_obs, M)
```
**Arguments:**
- `X0` --- initial particles
- `X_obs` --- observed dataset with missing values
- `M` --- missingness mask (0 = missing, 1 = observed)

Additional arguments allow tuning the hyperparameters, though we found the defaults to work well across a wide range of datasets.

## Empirical results in the paper

The experimental results in the paper can be reproduced via two notebooks:
- `main_realdata.ipynb` --- runs and evaluates different imputation/generation methods on real datasets
- `main_simulations.ipynb` --- runs and evaluates methods for the simulation study described in the paper

The R script `main_realdata.R` prepares the real datasets by amputing them and inserting missing values, while `load_uci_datasets.ipynb` can be used to download the UCI datasets. Both, typically do not need to be re-run.

## Repository structure
- `Bayesian/` contains the implementation of the Bayesian approach from [Chérief-Abdellatif and Näf (2026)](https://arxiv.org/abs/2603.23449)
- `MIRI-Imputation/` contains the implementations of MIRI, MissDiff and NewImp from [Yu et al. (2025)](https://arxiv.org/abs/2505.11749)
- `datasets/` contains the real datasets (the complete ones as well as the fully-observed and amputed parts after splitting)
- `results/` contains the saved results (simulation results are included, while real data results are excluded due to file size, but can be reproduced as described above)

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

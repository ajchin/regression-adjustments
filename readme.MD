## Code for interference regression adjustments

This repository contains code for producing the simulations and figures found in 

* A. Chin, [Regression adjustments for estimating the global treatment effect in experiments with interference](https://www.degruyter.com/view/j/jci.ahead-of-print/jci-2018-0026/jci-2018-0026.xml), Journal of Causal Inference, May 2019.

### Data

The directories `cai-data` and `data` serve as placeholder folders for the data.

* `.mat` files from the Facebook100 (FB100) dataset should be placed in the `data` folder.  The Facebook100 (FB100) dataset is publicly available from the Internet Archive at https://archive.org/details/oxford-2005-facebook-matrix and other public repositories.
* `cai-data` should contain data from the paper ``Social networks and the decision to insure'' (2015).  The data is available for download from the [publication website](https://www.aeaweb.org/articles?id=10.1257/app.20130442).

## Reproducing results

* The `functions` folder contains useful helpers and utilities described as follows:

- `covariate_functions.R` contains a few sample functions for generating possible interference features.  These are the $X_i(\textbf{W_{-i})$ functions described in the paper.

* The `figures` directory contains all figures from the paper.
* The `results` directory is a placeholder location for the simulation output stored as `csv` files.
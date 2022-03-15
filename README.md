README Document for Reproducing Results in
==========================================
"Identification and estimation of multinomial choice models with latent special covariates"
=============================================
Nail Kashaev
nkashaev@uwo.ca

Data Sources
============

The paper uses one data sets from

-   Benoit, D. F., Van Aelst, S., and Van den Poel, D. (2016). Outlier-robust bayesian multinomial choice
modeling. Journal of Applied Econometrics, 31(7):1445–1466. <https://doi.org/10.1002/jae.2482>

Software
========

A version 1.6.4 of the `Julia` programming language was used in coding the analysis files. For details about how to install `Julia` on different platforms and make `Julia` programs executable from the command line see <https://julialang.org/downloads/platform/>. After installation of `Julia 1.6.4.` run `using Pkg` and `Pkg.instantiate()` in the `Julia` terminal after setting the replication folder as the main one.

Hardware
========

The code was run on Mac mini (M1, 2020) with 16 Gb of RAM

Content
=======

-   `Application`  -- the folder contains the analysis files to replicate the results in Section 6.

-   `Simulations`  -- the folder contains the analysis files to replicate the results in Section 5 and Appendix B.

-   `Manifest.toml` and `Project.toml`  -- toml files with all necessary `Julia` packages.


Below I describe the content of every folder.

`Applicaton`
============

`data`
-----------

-    `csv` files contain data on margarine purchases.

-   `creating_table2.jl`  -- the code generates Table 2 from page 17. 

`parametric_estimation`
-----------

-    `application_parametric_main.jl` -- the code generates the results from Section 6 (Parameric Estimation), pages 18-19.

-    `parametric_functions.jl` -- the functions used in `application_parametric_main.jl`.

`results`
-----------

This folder contains all estimation results from Section 6 and Table 2.

-    `estimates_and_se_param.csv` -- parametric estimates and standard errors.

-    `estimates_and_se.csv` -- semiparametric estimates and standard errors.

-    `table2.csv` -- Table 2 from oage 17.

`semiparametric_estimation`
-----------

-    `application_semiparametric_main.jl` -- the code generates the results from Section 6 (Semiparameric Estimation), pages 19-20.

-    `semiparametric_functions.jl` -- the functions used in `application_semiparametric_main.jl`.



`Simulations`
============

`logit estimator`
-----------
This folder constains simulation results used in construction of Table 4 on page 35
-    `Results` -- csv files with simulation results

-   `common_functions.jl`  -- the functions used in simulations.

-   `logit_main.jl`  -- the main code. 

`probit estimator`
-----------
This folder constains simulation results used in construction of Table 5 on page 35
-    `Results` -- csv files with simulation results

-   `common_functions.jl`  -- the functions used in simulations.

-   `probit_main.jl`  -- the main code. 


`semiparametric estimator`
-----------
This folder constains simulation results used in construction of Tables 1 and 3 on pages 16 and 35.
-    `Results` -- csv files with simulation results

-   `common_functions.jl`  -- the functions used in simulations.

-   `np_main.jl`  -- the main code. 

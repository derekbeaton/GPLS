Generalized Partial Least Squares
================

This repository is for our work on generalized partial least squares.
The preprint can be found here: [Beaton, D., ADNI, Saporta, G., Abdi, H.
(2019). A generalization of partial least squares regression and
correspondence analysis for categorical and mixed data: An application
with the ADNI data.
bioRxiv, 598888.](https://www.biorxiv.org/content/10.1101/598888v2).

The preprint was generated through RMarkdown with the `rticles` package,
and is fully reproducible *for us* at least\! You can reproduce it so
long as you have access to ADNI data. See the [Manuscript](/Manuscript/)
directory for the manuscript’s content and supplements (RMarkdown, PDF).

Other materials can be found in the [Code](/Code/) directory. A detailed
RMD file provides how to reproduce the data and is found
[here](/Code/prep_data_code_notes.rmd) with an HTML
[file](/Code/prep_data_code_notes.html).

The GPLS package lives in the [Package](/Package/) directory. It is
documented and has a dependency on the GSVD package (found
[here](https://github.com/derekbeaton/GSVD)). See the
[Package](/Package/) directory for installation instructions.

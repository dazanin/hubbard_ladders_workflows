WORKFLOW SCRIPTS FOR "Pair Correlations in Doped Hubbard Ladders"
=================================================================


1. CONTENTS
-----------

This package contains the evaluation scripts for reproducing the figures and
analysis in the paper M. Dolfi, B. Bauer, S. Keller, M. Troyer, "Pair
Correlations in Doped Hubbard Ladders" (to appear).


The package is organized as follow:
```
-- LICENSE_1_0.txt                         # License file
-- README.md                               # This readme
-- data_evaluated/                         # Folder containing evaluated data in
                                           # text format
-- data_raw/                               # Folder containing raw data in hdf5
                                           # format
-- density_correlations.ipynb              # (**) IPython notebook for
                                           # density-density correlations
-- download_evaluated_data.sh              # Script downloading a preset of
                                           # evaluated data (see below)
-- download_raw_data.sh                    # Script downloading the raw data
                                           # (see below)
-- energy.ipynb                            # (**) IPython notebook for the
                                           # energy extrapolation
-- local_density.ipynb                     # (**) IPython notebook for local
                                           # density plots, as well as density amplidues analsys
-- pair_correlations.ipynb                 # (**) IPython notebook for pair
                                           # correlations
-- scripts/                                # Scripts for data evaluation and
                                           # plotting
-- scripts/amplitudes.py                   # (**) Performs the density-amplitudes
                                           # analysis
-- scripts/corr_helpers.py                 # helper functions
-- scripts/density.py                      # (**) Extracts the density and fit
                                           # the density
-- scripts/density_correlations.py         # (**) Extracts the density-density
                                           # correlations
-- scripts/energy.py                       # (**) Extract the energy
-- scripts/extrapolate.py                  # (**) Extrapolation functions for the
                                           # energy
-- scripts/extrapolate_local.py            # (**) Extrapolation functions for
                                           # vector observables (density and correlation functions)
-- scripts/load.py                         # (**) Load functions which load the
                                           # evaluated data or perform the data analysis if not yet cached
-- scripts/load_raw_data.py                # (*) Load functions which load the measurements
                                           # from the raw data
-- scripts/pairfield_correlations.py       # (**) Extracts the pairing
                                           # correlations
-- scripts/pyalps_dset/                    # subset of the pyalps library to perform data analysis
                                           # correlations
-- scripts/utils.py                        # Utilities

(*) requires pyalps in the Python PATH
(**) requires pyalps only if new obserservables have to be extracted from the raw data
```


2. DEPENDENCIES
---------------

For reproducing the figures in the paper from the evaluated data:
* Python 2.7
* IPython
* Numpy
* SciPy
* Matplotlib

Additionally for repeating the analysis from the raw data and/or extend the evaluation:
* Pyalps (http://alps.comp-phys.org)


3. OBTAIN THE DATA
------------------

We have two sets of data:
 - **Extracted data** is stored in `data_extracted/`. If the directory is empty,
  you can download pre-extracted observables from the data DOI.
  
  doi: http://dx.doi.org/10.7910/DVN/I5ANSU
  The script `download_extracted_data.sh` will perform the operation for you.
  
  Extracted data contain the all observables used in the paper in text format.
  Correlation functions have already been exported in the form C(M, r) with M
  being the bond dimension and r=|i-j| the distance between two rungs. As
  described in the paper we average over several pair i,j. Correlation
  functions for r=|i-j| with a fixed i=18 or i=20 are also exported to
  reproduce the comparison in the paper.
  
  More observables or other correlation function analysis can be extracted
  from the raw data (see next part).
  
  Data is stored in text format which can be read with your preferred tool,
  the initial lines starting with # are inteded as comments.

 - **Raw data** is stored in `data_raw/`. If the directory is empty, you need to
  download the raw data from the the data DOI.

  doi: http://dx.doi.org/10.7910/DVN/I5ANSU
  The script `download_raw_data.sh` will perform the operation for you.
  
  Raw data is stored in HDF5 according to the ALPS Schema
  (http://alps.comp-phys.org). The loader script in `scripts/load_raw_data.py`
  is an example on how to read the raw data.

In this package we provide data for these parameters:
* Bond dimension M = 800, 1200, 1600, 2000, 2800, 3200, 3600, 4000, 4800
* System size L = 32, 48, 64, 80, 96, 128, 160, 192
* Averange filling n = 0.875, 0.9375


4. START THE WORKFLOWS
----------------------

Evaluation workflows to reproduce the figures in the paper as available as
IPython notebooks and simple Python scripts. In what follows we describe only
how to run the notebooks; for the Python scripts the procedure is analogous.

1. Launch the notebook from the root of this package:
  ```bash
  ipython notebook
  ```
  This should open a browser window, where you can see the content of the
package.

2. Open the notebooks (files ending in .ipynb) by clicking on them.


5. LICENSE
----------

Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)

We ask that you acknowledge the use of our data and/or analysis workflows
by citing the following paper:
M. Dolfi, B. Bauer, S. Keller, M. Troyer, "Pair Correlations in Doped Hubbard
Ladders" (to appear).


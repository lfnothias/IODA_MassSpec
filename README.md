# IODA Mass Spec

This repository offers a collection of notebooks and scripts to run Iterative Optimised Data Acquisition (IODA), including MS2Planner - [(See article, Zuo, Cao, Nothias, Mohimani, *Bioinformatics*, Volume 37, Issue Supplement_1, July 2021, Pages i231–i236](https://academic.oup.com/bioinformatics/article/37/Supplement_1/i231/6319686)) and  [see GitHub repository](https://github.com/mohimanilab/MS2Planner).

### A. Generation custom exclusion list for IODA on Orbitrap MS

1. Feature detection and alignemnent from mzML/RAW files with pyOpenMS.
2. Create an exclusion list for the Exactive/Fusion serie MS or Exploris/Trihybrid serie MS.
3. Run the IODA experiments on the mass spectrometer.



### B. Generate a targeted list for targeted IODA experiments on Orbitrap MS

1. Feature detection and alignemnent from mzML/RAW files with pyOpenMS.
2. Optimise target ion distribution accross multiple iterative experiments.
3. Create an targeted list for the Exactive/Fusion serie MS or Exploris/Trihybrid serie MS.
4. Run the IODA experiments on the mass spectrometer.


## Running on your own data

Start the notebook on the cloud with Binder (GESIS server) -> [![Binder with GESIS](https://mybinder.org/badge_logo.svg)](https://notebooks.gesis.org/binder/v2/gh/lfnothias/IODA_MassSpec/2303_pyopenms?urlpath=lab/tree/IODA_notebooks_welcome.ipynb)

Binder with GESIS server [FROZEN] -> [![Binder with GESIS](https://mybinder.org/badge_logo.svg)](https://notebooks.gesis.org/binder/v2/gh/lfnothias/IODA_MassSpec/53276f8e456d0ab788da3cc9b33059934622e507)

**Interactive interface with Binder**

Binder standard server [LATEST] -> [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lfnothias/IODA_MassSpec/master?urlpath=lab/tree/home/jovyan/IODA_notebooks_welcome.ipynb)

Binder standard server [FROZEN] -> [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lfnothias/IODA_MassSpec/53276f8e456d0ab788da3cc9b33059934622e507)

**View-only interface** (non-interactive): [`IODA_notebooks_welcome.ipynb`](https://nbviewer.jupyter.org/github/lfnothias/IODA_MassSpec/blob/master/IODA_notebooks_welcome.ipynb)

You will need to do two things:

1. Open `IODA_TOPPAS_mztab_generation.ipynb` to run the feature finding to create the mztab file
2. Open `IODA_exclusion_from_mztab.ipynb` to create the exclusion file to download


## Testing

![Unit Test](https://github.com/lfnothias/IODA_MassSpec/workflows/Unit%20Test/badge.svg)

Unit tests are run using github actions. To run them manually:

```make test-unit```

The actual tests are in the ```test``` directory.

## Running the Docker in Linux

Run the following:

`make build-standalone`

`make run_standalone_notebook`

Access the notebook at `localhost:9000`

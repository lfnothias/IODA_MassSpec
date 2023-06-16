# IODA Mass Spec

This repository offers a collection of notebooks and scripts to run Iterative Optimised Data Acquisition (IODA), including MS2Planner - [(See article, Zuo, Cao, Nothias, Mohimani, *Bioinformatics*, Volume 37, Issue Supplement_1, July 2021, Pages i231–i236](https://academic.oup.com/bioinformatics/article/37/Supplement_1/i231/6319686)) and  [see GitHub repository](https://github.com/mohimanilab/MS2Planner).

### A. Generation custom exclusion list for IODA on Orbitrap MS

1. Feature detection and alignemnent from mzML/RAW files with pyOpenMS.
2. Create an exclusion list for the Exactive/Fusion serie MS or Exploris/Trihybrid serie MS.
3. Run the IODA experiments on the mass spectrometer.



### B. Generate a targeted list for targeted IODA experiments on Orbitrap MS

1. Feature detection and alignemnent from mzML/RAW files with pyOpenMS.
2. Optimise target ion distribution accross multiple iterative experiments.
3. Create an targeted list for the Exactive/Fusion serie MS or Exploris/Tribrid serie MS.
4. Run the IODA experiments on the mass spectrometer.


# Running IODA on the cloud with Binder

Note that this a temporary web-instance. Nothing will be saved on the server. Download all required results.

**Interactive interface with Binder-GESIS [First choice]** 

Start the notebook on the cloud with Binder on the GESIS server (more RAM for the Curve mode) -> [![Binder with GESIS](https://mybinder.org/badge_logo.svg)](https://notebooks.gesis.org/binder/v2/gh/lfnothias/IODA_MassSpec/master?urlpath=lab/tree/IODA_notebooks_welcome.ipynb)

**Interactive interface with Binder [back up and no curve mode]**

Binder standard server [LATEST] -> [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lfnothias/IODA_MassSpec/master?urlpath=lab/tree/IODA_notebooks_welcome.ipynb)

**View-only interface** (non-interactive): [`IODA_notebooks_welcome.ipynb`](https://nbviewer.jupyter.org/github/lfnothias/IODA_MassSpec/blob/master/IODA_notebooks_welcome.ipynb)

# Running IODA locally

First, if not already installed, install [git](https://git-scm.com/downloads) and [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Clone or download this repository. From a terminal, do:

```
git clone --recursive https://github.com/lfnothias/IODA_MassSpec.git
```

Navigate to the `IODA_MassSpec` folder:

```
cd IODA_MassSpec
```
Create the conda environment using the conda env create command with your environment.yml file:

```
conda env create -f environment_IODA_MassSpec.yml
```
Activate the IODA_MassSpec conda environment 
```
conda activate IODA_MassSpec
```

Now run the following to install pyOpenMS 3.0:

For Linux:

```
pip install pyopenms_wheels/pyopenms-3.0.0.dev20230316-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
```	


For Mac, see the correct wheel build and instructions at [https://pypi.cs.uni-tuebingen.de/](https://pypi.cs.uni-tuebingen.de/). For example for ARM:

```	
pip install pyopenms_wheels/pyopenms-3.0.0.dev20230615-cp310-cp310-macosx_13_0_arm64.whl
```	

For Windows, see the correct wheel build and instructions at [https://pypi.cs.uni-tuebingen.de/](https://pypi.cs.uni-tuebingen.de/)

Alternatively, go to the OpenMS action webpage [like here](https://github.com/OpenMS/OpenMS/actions/), download the pyOpenMS 3 wheel file corresponding to your system (note that you need to be logged in GitHub). Save it in the IODA_MassSpec folder, and install it similarly as above.


Start the jupyter lab session and notebook(you can use other ipython interface if needed):

```
jupyter lab IODA_notebooks_welcome.ipynb
```


## Update the submodule

```
git submodule update --remote
```

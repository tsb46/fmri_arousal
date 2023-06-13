# Global BOLD and Peripheral Physiology Repo
This repository contains the code necessary to replicate and extend the analyses of Bolt et al. (2023) - 'A Unified Physiological Process Links Global Patterns of Functional MRI, Respiratory Activity, and Autonomic Signaling' (https://www.biorxiv.org/content/10.1101/2023.01.19.524818v1). 

The primary analyses of the study are contained within the command-line Python utilities (see below), and the Jupyter notebook (analysis.ipynb) in the base directory of the repo. Eight separate fMRI datasets were used for the study and the default pipeline pulls and preprocesses all eight. If one would like to focus their analysis on a smaller subset of the datasets, the preprocessing and analysis scripts allow you to restrict your analyses to those datasets.

# Installation and Requirements

## Python  
The code in this repo was run with Python 3.9, but should be compatible with other versions (but see note below on MATLAB). 

In the base directory of the repo, pip install all necessary packages:
```
pip install -r requirements.txt
```

## MATLAB and EEGLAB

### MATLAB
For one dataset, NATVIEW, MATLAB is required for the EEGLAB preprocessing. **If the NATVIEW dataset is not used, there is no need for MATLAB**. Our EEG preprocessing was conducted in MATLAB R2023a. The matlab EEGLAB preprocesing script is integrated into the default Python preprocessing pipeline (preprocess.py) via the [MATLAB Engine API for Python](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html). This puts constraints on the Python version that must be used (see [MATLAB Engine Python Version Requirements](https://www.mathworks.com/support/requirements/python-compatibility.html)), and its implementation can be a headache. Alternatively, one can run the EEGLAB preprocessing for the NATVIEW dataset in a standalone MATLAB script ('utils/preprocess_natview.m') and pass this information to the Python preprocessing pipeline via a command-line argument (see below).

### EEGLAB







# Global BOLD and Peripheral Physiology Repo
This repository contains the code necessary to replicate and extend the analyses of Bolt et al. (2023) - 'A Unified Physiological Process Links Global Patterns of Functional MRI, Respiratory Activity, and Autonomic Signaling' (https://www.biorxiv.org/content/10.1101/2023.01.19.524818v1). 

The primary analyses of the study are contained within the command-line Python utilities (see below), and the Jupyter notebook (analysis.ipynb) in the base directory of the repo. Eight separate fMRI datasets were used for the study. If one would like to focus their analysis on a smaller subset of the datasets, the preprocessing and analysis scripts allow you to restrict your analyses to those datasets. 


# Installation and Requirements

## Python  
The code in this repo was run with Python 3.9, but should be compatible with other V3 versions (but see note below on MATLAB). 

In the base directory of the repo, pip install all necessary packages:
```
pip install -r requirements.txt
```

## FSL
The preprocessing script uses Nipype to call FSL utilities. Ensure that FSL is available from your command-line. Our study used FSL 5.0.9.

## MATLAB and EEGLAB

**Note: MATLAB is only needed for preprocessing of the NATVIEW dataset.**

### MATLAB
For the NATVIEW dataset, MATLAB is required for the EEGLAB preprocessing. Our EEG preprocessing was conducted in MATLAB R2023a. The matlab EEGLAB preprocesing script is integrated into the default Python preprocessing pipeline (preprocess.py) via the [MATLAB Engine API for Python](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html). To use the EEGLAB preprocessing script on the NATVIEW dataset, one additional pip install of the MATLAB engine API is necessary:

```
python -m pip install matlabengine
```

The MATLAB engine API puts constraints on the Python version that can be used (see [MATLAB Engine Python Version Requirements](https://www.mathworks.com/support/requirements/python-compatibility.html)), and its implementation can be a headache. Alternatively, one can run the EEGLAB preprocessing for the NATVIEW dataset in a standalone MATLAB script ('utils/preprocess_natview.m') and pass this information to the Python preprocessing pipeline via a command-line argument (see below). Instructions are provided in the MATLAB script.

### EEGLAB
[EEGLAB](https://github.com/sccn/eeglab) and the [fMRIB plugin](https://github.com/sccn/fMRIb) are required for NATVIEW EEG preprocessing. If the default Python pipeline (preprocess.py) is used with the MATLAB Engine API, the EEGLAB directory must be placed in the [external](external/) directory so that EEGLAB functions are readable from the script. Detailed instructions can be found in the README.md in the [external](external/) directory.

If you are using the standalone MATLAB script (see above), you place the EEGLAB directory in any location, provided that it is added to the MATLAB path before running the script.

# Get Data
Most datasets for this study are publically available (excluding the 'chang', 'chang_bh', 'chang_cue' datasets, but hosting solutions are being explored). A data pull utility is provided for half of the datasets (see below). At present, the others will need to manually requested/downloaded and placed into their appropriate location in the [data](data/) directory. 

Below is a table of the datasets used in the study. For each dataset, their repo label (the label used in the code), their manuscript label, their tasks, a documentation/manuscript link, the format they are downloaded in (raw vs. preprocessed; i.e. where the preprocessing starts from), and whether they are pulled in the data pull utility script:


| Dataset (Repo Label) | Paper Label  | Tasks                      | Link                                                               | pull_data.py       | format       |
| -------------------- | ------------ | -------------------------- | ------------------------------------------------------------------ | ------------------ | ------------ |
|        chang         | ME-REST      | Resting-state              | [link](https://elifesciences.org/articles/62376)                   |                    | preprocessed |
|        chang_bh      | ME-TASK      | Auditory Cue - Deep Breath |                                                                    |                    | preprocessed |
|        chang_cue     | ME-TASK-CUE  | Auditory Cue - No Breath   |                                                                    |                    | preprocessed |
|        hcp           | HCP-REST     | Resting-state              | [link](https://www.humanconnectome.org/study/hcp-young-adult)      | :white_check_mark: | preprocessed |
|        natview       | NATVIEW-REST | Resting-state              | [link](https://www.biorxiv.org/content/10.1101/2022.11.23.517540v1)| :white_check_mark: | raw          |
|        nki           | NKI-TASK     | Breath Hold                | [link](http://fcon_1000.projects.nitrc.org/indi/enhanced/)         | :white_check_mark: | raw          |
|        spreng        | ME-REST-SUPP | Resting-state              | [link](https://openneuro.org/datasets/ds003592/versions/1.0.13)    |                    | preprocessed |
|        yale          | YALE-REST    | Resting-state              | [link](https://openneuro.org/datasets/ds003673/versions/2.0.1)     | :white_check_mark: | raw          |

To pull the hcp, natview, nki and yale datasets, move into the [data](data/) directory and run the following command-line Python script:

```
cd data
python pull_data.py -d all
```

To pull a specific dataset, simply provide the repo label of the dataset (e.g. natview) to the '-d' argument:

```
python pull_data.py -d repo_label
```

# Preprocess
A single command-line Python script is provided to preprocess all datasets:

```
python preprocess.py -d all
```

The above script will assume that all datasets have been downloaded. If you want to preprocess a specific dataset, provide the repo label of the dataset to the '-d' argument:

```
python preprocess.py -d repo_label
```

# Analysis
The majority of the analyses reported in the manuscript are performed in the [analysis.ipynb](analysis.ipynb) Jupyter notebook. In order to run the analyses in the notebook, outputs from PCA and CPCA are needed. To get all the necessary outputs:

```
python run_all.py
```
**Note:** this script assumes all datasets have been downloaded and preprocessed. If you are interested in a subset of datasets, simply comment out/remove the dataset labels in the run_all.py script:
```
# specify datasets
datasets = [
	'chang', 
	'chang_bh', 
	'chang_cue', 
	'natview',
	'nki', 
	'hcp', 
	'spreng', 
	'yale'
]
```

**Note:** this script also runs several supplementary analyses, including physiological cross-correlations, peak-averaging of BOLD time courses around physiological signal peaks, and block averaging of task-fMRI datasets (nki, chang_bh, and chang_cue). These are not strictly necessary for the analysis.ipynb. If you are not interested in these outputs, you can comment out these lines from the run_all.py script.



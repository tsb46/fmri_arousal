# Global BOLD and Peripheral Physiology Repo
This repository contains the code necessary to replicate and extend the analyses of Bolt et al. (2024) - 'Widespread neural and autonomic system synchrony across the brain-body axis' (https://www.biorxiv.org/content/10.1101/2023.01.19.524818v2). 

The primary analyses of the study are contained within the command-line Python utilities (see below), and the Jupyter notebook (analysis.ipynb) in the base directory of the repo. Ten separate fMRI datasets were used for the study. If one would like to focus their analysis on a smaller subset of the datasets, the preprocessing and analysis scripts allow you to restrict your analyses to those datasets. 


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

## FMRI Datasets

Most FMRI datasets for this study are publically available, excluding the 'chang' datasets (hosting solutions are being explored) and 'toronto' datasets (available through data use agreement). A data pull utility (pull_data.py) is provided for publically-available datasets.

Below is a table of the datasets used in the study. For each dataset, their repo label (the label used in the code), their manuscript label, their tasks, a documentation/manuscript link, the format they are downloaded in (raw vs. preprocessed; i.e. where the preprocessing starts from), and whether they are pulled in the data pull utility script:


| Dataset (Repo Label) | Paper Label  | Tasks                      | Link                                                               | pull_data.py       | format       |
| -------------------- | ------------ | -------------------------- | ------------------------------------------------------------------ | ------------------ | ------------ |
|        chang         | ME-REST      | Resting-state              | [link](https://elifesciences.org/articles/62376)                   |                    | raw |
|        chang_bh      | ME-TASK      | Auditory Cue - Deep Breath |                                                                    |                    | raw |
|        chang_cue     | ME-TASK-CUE  | Auditory Cue - Button Response   |                                                                    |              | raw |
|        hcp           | HCP-REST     | Resting-state              | [link](https://www.humanconnectome.org/study/hcp-young-adult)      | :white_check_mark: | preprocessed |
|        natview       | NATVIEW-REST | Resting-state              | [link](https://www.biorxiv.org/content/10.1101/2022.11.23.517540v1)| :white_check_mark: | raw          |
|        nki           | NKI-TASK     | Breath Hold                | [link](http://fcon_1000.projects.nitrc.org/indi/enhanced/)         | :white_check_mark: | raw          |
|        nki_rest      | NKI-REST     | Resting-state              | [link](http://fcon_1000.projects.nitrc.org/indi/enhanced/)         | :white_check_mark: | raw          |
|        spreng        | ME-REST-SUPP | Resting-state              | [link](https://openneuro.org/datasets/ds003592/versions/1.0.13)    | :white_check_mark: | raw          |
|        yale          | YALE-REST    | Resting-state              | [link](https://openneuro.org/datasets/ds003673/versions/2.0.1)     | :white_check_mark: | raw          |
|        toronto       | Clamped CO2/Free Breathing   | Resting-state   | [link](https://www.sciencedirect.com/science/article/pii/S1053811920303608)     |  | raw          |

To pull the hcp, natview, spreng, nki and yale datasets, move into the [data](data/) directory and run the following command-line Python script:

```
cd data
python pull_data.py -d all
```

To pull a specific dataset, simply provide the repo label of the dataset (e.g. natview) to the '-d' argument:

```
python pull_data.py -d repo_label
```

## Atlases

In **Figure 3**, the spatiotemporal pattern of the global fMRI signal is compared to several publically-avalaible brain atlases, including probability atlases of large cerebral veins (Huck et al., 2019), arteries (Mouches & Forket, 2019), and cerebrospinal fluid (CSF; Lorio et al., 2016), as well as noradrenaline transporter (NAT) binding potential (Ding et al., 2010; Hansen et al., 2022). Links to download atlases:

* [Cerebral vein atlas ('PartialVolume')](https://figshare.com/articles/dataset/VENAT_Probability_map_nii_gz/7205960)
* [Cerebral artery atlas ('Vessel Occurence Probablity Atlas')](https://springernature.figshare.com/collections/A_statistical_atlas_of_cerebral_arteries_generated_using_multi-center_MRA_datasets_from_healthy_subjects/4215089)
* [Cerebrospinal Fluid Probability (CSF) Density (index=2, indexing starts at zero)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4819722/)
* [Noradrenaline transporter (NAT) binding potential](https://github.com/netneurolab/hansen_receptors)

Once downloaded, place these files into the 'results' directory for use in the [analysis.ipynb](analysis.ipynb) Jupyter notebook.

# Preprocess

## Standard Preprocessing

A single command-line Python script is provided to preprocess all datasets:

```
python preprocess.py -d all
```

The above script will assume that all datasets have been downloaded. If you want to preprocess a specific dataset, provide the repo label of the dataset to the '-d' argument:

```
python preprocess.py -d repo_label
```

If you decided to preprocess the NATVIEW EEG dataset with the stand-alone MATLAB script (utils/natview_preprocess.py), make sure you pass the '-no_eeglab' flag to the script:

```
python preprocess.py -d all -no_eeglab
```

## Multi-Echo T2* and SO Preprocessing
In **Figure 3** of the manuscript, signal decay (T2*) and initial signal intensity (S0) spatiotemporal dynamics of the global fMRI signal are estimated from multiecho fMRI datasets (ME-REST and ME-REST-SUPP). Multiecho estimation and preprocessing of T2* and S0 signals is performed by the command-line script 'multiecho_pipe.py' in the 'multiecho' directory. To run the multiecho preprocessing pipeline (from the 'multiecho' directory) on a multiecho fMRI dataset ('dataset'):

```
python multiecho_pipe.py -d 'dataset'
```


# Analysis
The majority of the analyses reported in the manuscript are performed in the [analysis.ipynb](analysis.ipynb) Jupyter notebook. In order to run the analyses in the notebook, outputs from PCA, CPCA are needed. To get all the necessary outputs:

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
	'nki_rest',
	'hcp', 
	'spreng', 
	'yale',
	'toronto'
]
```

# References

Ding, Y.-S., Singhal, T., Planeta-Wilson, B., Gallezot, J.-D., Nabulsi, N., Labaree, D., Ropchan, J., Henry, S., Williams, W., Carson, R. E., Neumeister, A., & Malison, R. T. (2010). PET imaging of the effects of age and cocaine on the norepinephrine transporter in the human brain using (S,S)-[11C]O-methylreboxetine and HRRT. Synapse, 64(1), 30–38. https://doi.org/10.1002/syn.20696

Hansen, J. Y., Shafiei, G., Markello, R. D., Smart, K., Cox, S. M. L., Nørgaard, M., Beliveau, V., Wu, Y., Gallezot, J.-D., Aumont, É., Servaes, S., Scala, S. G., DuBois, J. M., Wainstein, G., Bezgin, G., Funck, T., Schmitz, T. W., Spreng, R. N., Galovic, M., … Misic, B. (2022). Mapping neurotransmitter systems to the structural and functional organization of the human neocortex. Nature Neuroscience, 25(11), 1569–1581. https://doi.org/10.1038/s41593-022-01186-3

Huck, J., Wanner, Y., Fan, A. P., Jäger, A.-T., Grahl, S., Schneider, U., Villringer, A., Steele, C. J., Tardif, C. L., Bazin, P.-L., & Gauthier, C. J. (2019). High resolution atlas of the venous brain vasculature from 7 T quantitative susceptibility maps. Brain Structure & Function, 224(7), 2467–2485. https://doi.org/10.1007/s00429-019-01919-4

Lorio, S., Fresard, S., Adaszewski, S., Kherif, F., Chowdhury, R., Frackowiak, R. S., Ashburner, J., Helms, G., Weiskopf, N., Lutti, A., & Draganski, B. (2016). New tissue priors for improved automated classification of subcortical brain structures on MRI. Neuroimage, 130, 157–166. https://doi.org/10.1016/j.neuroimage.2016.01.062

Mouches, P., & Forkert, N. D. (2019). A statistical atlas of cerebral arteries generated using multi-center MRA datasets from healthy subjects. Scientific Data, 6, 29. https://doi.org/10.1038/s41597-019-0034-5




# Download EEGLAB and FMRIB Plug-in
In order to preprocess the simultaneous EEG data from the NATVIEW dataset, we need to use a MATLAB script that calls functions from EEGLAB and one of its plugins (FMRIB). Of course, you will need access to MATLAB, which requires a license. We looked for open-source alternatives in Python, but based on a limited search, we found no production-ready gradient-artifact correction tools (though this could change at any moment).

We can download EEGLAB through Git (or manually). Ensure that you are in the 'external' directory and use this command to clone the EEGLAB repo:

```
git clone --recurse-submodules https://github.com/sccn/eeglab.git
```

To download the FMRI plugin, and place in the correct location expected by EEG lab, run the following command (again, ensure you are in the 'external' directory):

```
git -C eeglab/plugins/ clone https://github.com/sccn/fMRIb.git 
```

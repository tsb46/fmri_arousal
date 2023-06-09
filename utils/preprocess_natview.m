%% Description
% There are two ways to preprocess the NATVIEW eeg data (.set files).
% NATVIEW. The most straightforward way is use the 'preprocess.py' script
% in the main directory. This script uses the Matlab engine API for Python:
%
% https://www.mathworks.com/help/matlab/matlab-engine-for-python.html
%
% The script calls the 'eeglab_preprocess.py' script from Python and no
% extra steps are needed. However, installation and set-up of the API can
% be a headache, depending on the version of MATLAB and Python installed
% locally. 
% 
% Alternatively one can run this script to preprocess the EEG data
% strictly in MATLAB, and one can finish the rest of the preprocessing of
% the NATVIEW dataset with the 'preprocess.py' script without calling the
% Matlab engine API. Simply call the 'preprocess.py' script with the
% following command in the terminal after you have run this MATLAB script:
%
% python preprocess.py -d natview -no_eeglab
%
% Important Note: run this script from this directory (utils)
%% Add path of base directory and eeglab directory
addpath(genpath('../external/eeglab'))
%% load natview subject file and define directories
subjects = readtable('../data/dataset_natview/subject_list_natview.csv');
input_dir = '../data/dataset_natview/eeg/raw';
output_dir = '../data/dataset_natview/eeg/proc1_eeg';
eeg_format = 'sub-%s_ses-0%d_task-rest_eeg.set';

%% Loop through subjects and preprocess
for row = 1:height(subjects)
    if subjects{row,1} < 10
        subj = strcat('0',num2str(subjects{row,1}));
    else
        subj = num2str(subjects{row,1});
    end
    % get eeg file path
    eeg_fp = sprintf(eeg_format,subj, subjects{row,2});
    eeg_in = strcat(input_dir,'/',eeg_fp);
    % run eeglab preprocessing
    eeglab_preprocess(eeg_in, output_dir);
end
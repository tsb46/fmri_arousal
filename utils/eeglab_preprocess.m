function EEG = eeglab_preprocess(fileNameSET,outputDir)
%% Description
% This is a modified preprocessing script derived from the NATVIEW team's
% eeg preprocessing pipeline:
% 
% https://github.com/NathanKlineInstitute/NATVIEW_EEGFMRI/blob/main/eeg/natview_eeg_preprocess_pipeline.m
% 
% We have modified the script to suit the objectives of our project.
% Specifically, to speed up preprocessing we downsample from 5000 to 500 Hz
% before preprocessing. In addition, we do not apply ASR artifact
% rejection to ensure no data segments are excluded.

%--------------------------------------------------------------------------
% INPUT:
%       fileNameSET - Filename of SET file
%
%         outputDir - Output directory for preprocessed EEG file(s)

%% Set environment
addpath(genpath('external/eeglab')) % path to EEGLAB

%% Input Handling
[~,filename,~] = fileparts(fileNameSET);

%% STEP 0: Load data into EEGLAB
% ensure chars are converted to strings
EEG = pop_loadset(fileNameSET); % Load SET file into MATLAB
EEG.data = double(EEG.data);

%% Non-EEG Channel specification
ECGChan = find(strcmp({EEG.chanlocs.labels},'ECG'));
EOGLChan = find(strcmp({EEG.chanlocs.labels},'EOGL'));
EOGUChan = find(strcmp({EEG.chanlocs.labels},'EOGU'));
electrodeExclude = [ECGChan,EOGLChan,EOGUChan];

%% STEP 1: Resample before pre-processing to 500Hz
resample_freq = 500;
EEG = pop_resample(EEG,resample_freq);

%% STEP 2: Gradient Artifact Removal
EEG = pop_fmrib_fastr(EEG,[],[],[],'R128',1,0,[],[],[],[],electrodeExclude,'auto'); % Remove gradient artifact

%% STEP 3a: QRS Detection
% This step detects QRS complexes in the ECG channel. If the function fails
% to find QRS complexes, QRS detection is performed on every EEG channel;
% the channel chosen for QRS detection equals the mode of QRS counts
try
    EEG = pop_fmrib_qrsdetect(EEG,ECGChan,'QRS','no'); % FMRIB Toolbox QRS Detection
catch
    nChannels = EEG.nbchan;
    channelEEG = 1:nChannels;

    QRSCount = zeros(nChannels,1);
    channelError = zeros(nChannels,1);

    for nn = 1:nChannels
        try
            EEG_QRS = pop_fmrib_qrsdetect(EEG,channelEEG(nn),'QRS','no');
            eventLatency = extract_eventLatency(EEG_QRS,'QRS');

            if(length(eventLatency) > (EEG.xmax - 50) || nn ~= 32)
                QRSCount(nn) = length(eventLatency);
            end
        catch
            channelError(nn) = 1;
        end
    end

    channelEEG(QRSCount == 0) = [];
    QRSCount(QRSCount == 0) = [];

    [QRSCount_mode, QRSCount_modeNum] = mode(QRSCount);
    [QRSCount_sort, QRSCount_sort_idx] = sort(QRSCount);

    % Select mode of QRS count if 3 or more, else select median QRS count
    if(QRSCount_modeNum >= 3)
        QRSCount_mode_idx = find(QRSCount == QRSCount_mode);
    else
        if(length(QRSCount) == 1)
            QRSCount_mode_idx = 1;
        else
            QRSCount_mode_idx = QRSCount_sort_idx(find(diff(QRSCount_sort > median(QRSCount)))); %#ok<FNDSB>
        end
    end

    QRS_channel = channelEEG(QRSCount_mode_idx);

    EEG = pop_fmrib_qrsdetect(EEG,QRS_channel(1),'QRS','no'); % FMRIB Toolbox QRS Detection
    % disp(find(channelError==1));
end

%% STEP 3b: Pulse Artifact Removal
PAType = 'median'; % Template for pulse artifact (Default: median)

EEG = pop_fmrib_pas(EEG,'QRS',PAType); % Pulse Artifact removal


%% STEP 4: Downsample EEG data to 250Hz
resample_freq = 250;
EEG = pop_resample(EEG,resample_freq);

%% STEP 5: Remove non-EEG channels (i.e., EOG and ECG)
% first pop ECG and save out for further preprocessing
ECG = pop_select(EEG, 'channel', ECGChan);
ecg_out_split = split(filename, '_');
ecg_out_prefix = join(ecg_out_split(1:3), '_');
ecg_out = strcat(ecg_out_prefix{1}, '_ecg.set');
pop_saveset(ECG,'filename',ecg_out,'filepath',outputDir, 'savemode', 'onefile');

% remove non-EEG channels
EEG  = pop_select(EEG,'nochannel',electrodeExclude);

%% STEP 6: Bandpass filter data
freq_lo = 0.3;
freq_hi = 50;
EEG = pop_eegfiltnew(EEG,'locutoff',freq_lo,'hicutoff',freq_hi);

%% STEP 7: Remove bad channels
EEG = pop_clean_rawdata(EEG,'FlatlineCriterion',5,...
                            'ChannelCriterion',0.8,...
                            'LineNoiseCriterion',4,...
                            'Highpass',[0.75 1.25],...
                            'BurstCriterion','off',...
                            'WindowCriterion','off',...
                            'BurstRejection','off',...
                            'Distance','Euclidian',...
                            'WindowCriterionTolerances','off');

%% STEP 8: Rereference data using average reference
EEG = pop_reref(EEG,[]);
                    
%% STEP 9: Compute ICA, flat IC using ICLabel, and remove ICs highly correlated with muscle and eye artifacts
EEG = pop_runica(EEG,'icatype','runica','concatcond','on','options',{'pca',-1});
EEG = pop_iclabel(EEG,'default');
EEG = pop_icflag(EEG,[NaN NaN; 0.8 1; 0.8 1; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);
EEG = pop_subcomp(EEG,[]);

%% Write Out
pop_saveset(EEG,'filename',filename,'filepath',outputDir, 'savemode', 'onefile');

end
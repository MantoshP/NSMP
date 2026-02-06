%% NSMP Study: Comprehensive EEG processing pipeline 
% Processing steps:
% 1. Initiate 
%       - set up directories
%       - get & update INFO output structure
% 2. Import & preprocessing of continuous EEG data
%       - assignment of electrode coordinates
%        - crop around triggers, DC removal, downsampling
% 3. Reload artifact free continuous data post-frequency filtering, electrode interpolation & rereferencing on GUI
% 
% 4. Import task data
%      - attach Unique trial ID
% 
% 5. Classify epochs based on task and form updated datasets( 'Baseline','DelaySelect','NoTMSSelect','CatchSelect','DelayNonSelect','NoTMSNonSelect','CatchNonSelect')
%
% %%%% Pipeline splits into pre-TMS analysis and Post-TMS   
% 6. Time-frequency analysis analysis of beta power change
% 
% 7. TEP analysis

%% 1 - Initialize
clear; close all; clc;
fprintf('Section 1: Initialize\n')
[FOLDER, INFO, subject_idx] = NSMP_EEG_Init();
fprintf('Section finished.\n\n')

%% 2 - Import data
fprintf('Section 2: Import & preprocess continuous data\n')

% Import continuous EEG datasets, pre-process & save for letswave
[DATASET, INFO] = NSMP_EEG_Import_data(FOLDER, INFO, subject_idx);
% [DATASET, INFO] = NSMP_EEG_Import_data_MultipleFolders(FOLDER, INFO, subject_idx);
% [DATASET, INFO] = NSMP_EEG_Import_data_lessthan4triggers(FOLDER, INFO, subject_idx);

% Continue or save
continue_or_save(DATASET)

% === MANUAL LETSWAVE STAGE ===
fprintf('Open Letswave, inspect data, interpolate channels, rereference.\n');
fprintf('Save datasets with suffix reref.\n');
pause;   % or letswave + wait4files

% === 🔑 RELOAD STEP (THIS IS WHAT YOU ARE MISSING) ===
DATASET = NSMP_EEG_ReloadProcessed(FOLDER, DATASET, 'reref');

% Import task information
[DATASET, INFO] = NSMP_EEG_Import_task(FOLDER, INFO, DATASET, subject_idx);

% Classify epochs based on task

[DATASET] = NSMP_EEG_ClassifyEpochs2(FOLDER, INFO, DATASET, subject_idx);



% Example usage right after:
idx  = DATASET.epochs.DelaySelect;           % [block, epochIdx]
rows = idx(idx(:,1) == 1, :);
X    = DATASET.raw(1).data_seg(:, :, rows(:,2));  % [channels x time x nEpochs]
t    = DATASET.raw(1).t;



%% 3 - VEP check
fprintf('Section 3: Analyze visual-evoked response\n')
[DATASET, INFO, DATA] = NSMP_EEG_VEP(FOLDER, INFO, DATASET, subject_idx);

% Continue or save
continue_or_save(DATASET)

%% FUNCTIONS
function continue_or_save(DATASET)
    answer = questdlg('Will you continue processing or do you want to save DATASET for later?', 'Next step',...
    'Continue', 'Save & stop', 'Continue');
    switch answer
        case 'Continue'
            fprintf('Section finished.\n\n')
        case 'Save & stop'
            filename = sprintf('DATASET_%s.mat', date);
            save(filename, 'DATASET')
            fprintf('DATASET saved to current directory.\nSave to close.\n\n')
    end
end
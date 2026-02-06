function [DATASET, INFO] = NSMP_EEG_ClassifyEpochs2(FOLDER, INFO, DATASET, subject_idx)
% NSMP_EEG_ClassifyEpochs2
%
% Safe Letswave-compatible implementation with:
% - VEP-style preprocessing
% - Task-aligned epoching
% - Trial UID traceability WITHOUT touching header.events pre-segmentation

fprintf('\n--- Epoching, Classifying, Merging, Preprocessing (Task) ---\n');

% === SAFETY CHECK: ensure correct input stage ===
assert( ...
    isfield(DATASET, 'raw') && ...
    isfield(DATASET.raw, 'header') && ...
    isfield(DATASET.raw(1).header, 'name') && ...
    contains(lower(DATASET.raw(1).header.name), 'reref'), ...
    'ClassifyEpochs2 expects rereferenced data (header.name must contain ''reref'').' ...
);

%% PARAMETERS
params = setup_parameters();
params = validate_and_complete_params(params);

%% INIT
DATASET = initialize_epoch_structure(DATASET, params);
DATASET.task_epochs = struct();

addpath(genpath(fullfile(FOLDER.toolbox, 'letswave7-master')));

total_classified = 0;

%% ====================== BLOCK LOOP ======================
for block_idx = 1:length(DATASET.raw)

    fprintf('Processing Block %d:\n', block_idx);

    % Load continuous data
    lwdata.header = DATASET.raw(block_idx).header;
    lwdata.data   = DATASET.raw(block_idx).data;
    task_data     = DATASET.raw(block_idx).task;

    fs = 1 / lwdata.header.xstep;
    fprintf('  Sampling rate: %.0f Hz\n', fs);

 

    %% --- Segmentation (SAFE: no event modification before this) ---
    option = struct();
    option.event_labels = {params.trigger_code};
    option.x_start      = params.epoch_window(1);
    option.x_end        = params.epoch_window(2);
    option.x_duration   = diff(params.epoch_window);
    option.suffix       = params.suffix{2};
    option.is_save      = 0;
    lwdata = FLW_segmentation.get_lwdata(lwdata, option);

    %% --- DC + linear detrend (epoched) ---
    option = struct();
    option.linear_detrend = 1;
    option.suffix         = params.suffix{3};
    option.is_save        = 0;
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option);

    %% Store epoched data
    DATASET.raw(block_idx).epoched.header = lwdata.header;
    DATASET.raw(block_idx).epoched.data   = lwdata.data;

    %% === Trial UID alignment (SAFE) ===
    n_epochs = size(lwdata.data, 1);
    n_task   = length(task_data);
    n_trials = min(n_epochs, n_task);

    DATASET.raw(block_idx).epoched.trial_uid = {task_data(1:n_trials).trial_uid};

    fprintf('  Segmented epochs: %d | Task trials: %d | Using: %d\n', ...
        n_epochs, n_task, n_trials);

    %% === Classification ===
    block_classified = 0;

    for trial_idx = 1:n_trials
        try
            trial_info = task_data(trial_idx);
            label = classify_trial_your_rules(trial_info);

            if ~isempty(label)
                DATASET.epochs.(label) = ...
                    [DATASET.epochs.(label); block_idx, trial_idx];
                block_classified = block_classified + 1;
            else
                if ~isfield(DATASET.epochs, 'Unclassified')
                    DATASET.epochs.Unclassified = [];
                end
                DATASET.epochs.Unclassified = ...
                    [DATASET.epochs.Unclassified; block_idx, trial_idx];
            end
        catch ME
            fprintf('    Trial %d: %s\n', trial_idx, ME.message);
        end
    end

    fprintf('  Classified %d trials in this block\n\n', block_classified);
    total_classified = total_classified + block_classified;
end

%% SUMMARY
display_classification_summary(DATASET, total_classified);

%% OUTPUT DIRS
FOLDER = create_output_directories(FOLDER, subject_idx);

%% ====================== MERGING ======================
fprintf('\nMERGING + PREPROCESSING PER CATEGORY:\n');
fprintf('=====================================\n');

for c = 1:length(params.categories)

    category = params.categories{c};
    idx = DATASET.epochs.(category);

    if isempty(idx)
        fprintf('  Skipping %s\n', category);
        continue;
    end

    fprintf('  Category %s: %d epochs\n', category, size(idx,1));
    lwdataset = build_lwdataset_for_category(DATASET, category);

    option = struct();
    option.type    = 'epoch';
    option.suffix  = [params.suffix{4} '_' category];
    option.is_save = 0;
    lw_merged = FLW_merge.get_lwdata(lwdataset, option);

    lw_merged.header.name = sprintf('Subject_%03d_%s', subject_idx, category);
    DATASET.task_epochs.(category) = lw_merged;

    save_lwdata_pair(lw_merged, FOLDER.epochs, lw_merged.header.name);
end

fprintf('--- COMPLETED ---\n');
end
%% =======================================================================
%% HELPER FUNCTIONS (CLEAN – SINGLE DEFINITIONS)
%% =======================================================================

function params = setup_parameters()
% Default preprocessing + epoching parameters

    params.bandpass     = [0.1 45];   % Hz
    params.filter_order = 4;           % Butterworth
    params.trigger_code = 'start';     % trial start trigger

    params.epoch_window = [0 5.5];     % seconds

    params.baseline_window = [];       % disabled

    params.categories = {
        'Baseline', ...
        'DelaySelect', 'DelayNonSelect', ...
        'NoTMSSelect', 'NoTMSNonSelect', ...
        'CatchSelect', 'CatchNonSelect'
    };

    params.suffix = {'bandpass','ep','dc','m','bl','reref'};
end

% -----------------------------------------------------------------------

function params = validate_and_complete_params(params)
    if ~isfield(params,'baseline_window')
        params.baseline_window = [];
    end
    if ~isfield(params,'filter_order')
        params.filter_order = 4;
    end
    if ~isfield(params,'suffix') || numel(params.suffix) < 6
        params.suffix = {'bandpass','ep','dc','m','bl','reref'};
    end
end

% -----------------------------------------------------------------------

function DATASET = initialize_epoch_structure(DATASET, params)
% Initialise epoch index containers

    DATASET.epochs = struct();
    for i = 1:numel(params.categories)
        DATASET.epochs.(params.categories{i}) = [];
    end
end

% -----------------------------------------------------------------------

function lwdataset = build_lwdataset_for_category(DATASET, category)
% Build Letswave dataset per category

    lwdataset = struct('header',{},'data',{});
    idx_table = DATASET.epochs.(category);

    if isempty(idx_table)
        return;
    end

    blocks = unique(idx_table(:,1))';
    k = 0;

    for b = blocks
        ep_idx = idx_table(idx_table(:,1)==b,2);

        if ~isfield(DATASET.raw(b),'epoched')
            continue;
        end

        lw_block.header = DATASET.raw(b).epoched.header;
        lw_block.data   = DATASET.raw(b).epoched.data;

        ep_idx = ep_idx(ep_idx>=1 & ep_idx<=size(lw_block.data,1));
        if isempty(ep_idx)
            continue;
        end

        lw_block.data = lw_block.data(ep_idx,:,:,:,:,:);
        lw_block.header.datasize = size(lw_block.data);
        lw_block.header.name = sprintf('%s_block%02d',category,b);

        if isfield(lw_block.header,'events')
            lw_block.header.events = lw_block.header.events(ep_idx);
        end

        k = k + 1;
        lwdataset(k) = lw_block; %#ok<AGROW>
    end
end

% -----------------------------------------------------------------------

function FOLDER = create_output_directories(FOLDER, subject_idx)

    subject_dir = fullfile(FOLDER.processed, ...
        sprintf('Subject_%03d',subject_idx));

    if ~exist(subject_dir,'dir')
        mkdir(subject_dir);
    end

    epochs_dir = fullfile(subject_dir,'Classified_Epochs_MergedPreproc');
    if ~exist(epochs_dir,'dir')
        mkdir(epochs_dir);
    end

    FOLDER.subject = subject_dir;
    FOLDER.epochs  = epochs_dir;
end

% -----------------------------------------------------------------------

function success = save_lwdata_pair(lwdata, out_dir, base_name)

    success = false;
    try
        header = lwdata.header;
        data   = lwdata.data;

        save(fullfile(out_dir,[base_name '.lw6']),'header','-v7.3');
        save(fullfile(out_dir,[base_name '.mat']),'data','-v7.3');

        fprintf('    Saved: %s\n',base_name);
        success = true;
    catch ME
        fprintf('    Save error (%s): %s\n',base_name,ME.message);
    end
end

% -----------------------------------------------------------------------

function label = classify_trial_your_rules(trial)
% *** YOUR RULES – UNCHANGED ***

    label = '';

    required = {'linesTTL','timeTTL','Side','isGo','correct'};
    for f = 1:numel(required)
        if ~isfield(trial,required{f}); return; end
    end

    linesTTL = double(trial.linesTTL);
    timeTTL  = trial.timeTTL;
    Side     = trial.Side;
    isGo     = trial.isGo;
    correct  = trial.correct;

    if linesTTL==1 && strcmp(timeTTL,'Baseline') && isGo==1
        label = 'Baseline';

    elseif linesTTL==1 && strcmp(timeTTL,'Delay') && strcmp(Side,'R') ...
            && isGo==1 && correct==1
        label = 'DelaySelect';

    elseif linesTTL==1 && strcmp(timeTTL,'Delay') && strcmp(Side,'L') ...
            && isGo==1 && correct==1
        label = 'DelayNonSelect';

    elseif linesTTL==0 && strcmp(timeTTL,'Baseline') && strcmp(Side,'R') ...
            && isGo==1 && correct==1
        label = 'NoTMSSelect';

    elseif linesTTL==0 && strcmp(timeTTL,'Baseline') && strcmp(Side,'L') ...
            && isGo==1 && correct==1
        label = 'NoTMSNonSelect';

    elseif strcmp(Side,'R') && isGo==0 && isfield(trial,'result') ...
            && strcmp(trial.result,'good catch')
        label = 'CatchSelect';

    elseif strcmp(Side,'L') && isGo==0 && isfield(trial,'result') ...
            && strcmp(trial.result,'good catch')
        label = 'CatchNonSelect';
    end
end

% -----------------------------------------------------------------------

function display_classification_summary(DATASET,total)

    fprintf('\nCLASSIFICATION SUMMARY\n=======================\n');
    cats = fieldnames(DATASET.epochs);

    for i = 1:numel(cats)
        n = size(DATASET.epochs.(cats{i}),1);
        if n>0
            fprintf('%-20s: %3d\n',cats{i},n);
        end
    end
    fprintf('-----------------------\nTotal classified: %d\n',total);
end

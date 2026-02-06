function [DATASET, INFO] = NSMP_EEG_Import_task(FOLDER, INFO, DATASET, subject_idx)
    % Study parameters
    study = 'NSMP';
    output_file = fullfile(FOLDER.output, [study '_output.mat']);
    
    % Identify import folder
    session_folders = dir(sprintf('%s\\%s', FOLDER.raw, INFO(subject_idx).ID));
    for a = 1:length(session_folders)
        if contains(session_folders(a).name, 'task', 'IgnoreCase', true)
            params.folder = sprintf('%s\\%s', session_folders(a).folder, session_folders(a).name);
        end
    end

    % Identify recording blocks
    session_blocks = dir(params.folder);
    data_idx = logical([]);
    for b = 1:length(session_blocks)   
        if ~contains(session_blocks(b).name, 'NSMP')
            data_idx(b) = true;
        elseif contains(session_blocks(b).name, 'tsv')
            data_idx(b) = true;
        else
            data_idx(b) = false;
        end
    end
    session_blocks(data_idx) = [];

      if length(session_blocks) ~= length(DATASET.raw)
        error('ERROR: number of task blocks (%d) does not match EEG blocks (%d)', ...
            length(session_blocks), length(DATASET.raw));
      end

    % Extract task info and append it to DATASET
    for c = 1:length(session_blocks)
        % Load trial info
        load(sprintf('%s\\%s', params.folder, session_blocks(c).name), 'trials')
        fields = fieldnames(trials);
        n = numel(trials.(fields{1}));
        if n == INFO(subject_idx).EEG.dataset(c).trials            
            trials_final(1:n) = struct();
            for i = 1:n
                for f = 1:numel(fields)
                    trials_final(i).(fields{f}) = trials.(fields{f})(i);
                end
            end
        else
            error('ERROR: number of trials in block %d does not match with trials in task data (%d)!', c, n)
        end

        % Attach to DATASET
        DATASET.raw(c).task = trials_final;

        % === Assign persistent trial UID (for traceability) ===
        for t = 1:length(DATASET.raw(c).task)
            DATASET.raw(c).task(t).trial_uid = ...
                sprintf('S%03d_B%02d_T%03d', subject_idx, c, t);
        end

        clear trials_final
    end

    % Save INFO
    save(output_file, 'INFO','-append')
end
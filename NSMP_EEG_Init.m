function [FOLDER, INFO, subject_idx] = NSMP_EEG_Init()
    % Get directory paths
    FOLDER.toolbox = uigetdir(pwd, 'Choose the toolbox folder');
    FOLDER.raw = uigetdir(pwd, 'Choose the raw data folder');
    FOLDER.processed = uigetdir(pwd, 'Choose the processed data folder');
    FOLDER.output = uigetdir(pwd, 'Choose the output folder');
    
    % Add Letswave to path
    addpath(genpath(fullfile(FOLDER.toolbox, 'letswave7-master')));
    
    % Study parameters
    study = 'NSMP';
    output_file = fullfile(FOLDER.output, [study '_output.mat']);
    
    % Initialize output structures if they don't exist
    if exist(output_file, 'file')
        load(output_file);
    else
        INFO = struct;
        DATA = struct;
        MEASURES = struct;
        save(output_file, 'INFO', 'DATA', 'MEASURES');
    end
    
    % Get subject info
    prompt = {'Subject number:'};
    dlgtitle = 'Subject Information';
    dims = [1 35];
    definput = {''};
    info = inputdlg(prompt, dlgtitle, dims, definput);
    subject_idx = str2double(info{1});
%     if ~isfield(INFO, 'ID') || isempty(INFO(subject_idx).ID)
    if length(INFO) < subject_idx || ~isfield(INFO(subject_idx), 'ID') || isempty(INFO(subject_idx).ID)

        % Fill in the info
        prompt = {'Age:', 'Sex:', 'Handedness:', 'rMT:', 'Session:'};
        definput = {'', 'F/M', 'R/L', '', datestr(now, 'dd-mmm-yyyy')};
        info = inputdlg(prompt, dlgtitle, dims, definput);
        INFO(subject_idx).ID = sprintf('%s%s', study, num2str(subject_idx, '%03d'));
        INFO(subject_idx).session = info{5};
        INFO(subject_idx).age = str2double(info{1});
        INFO(subject_idx).sex = info{2};
        INFO(subject_idx).handedness = info{3};
        INFO(subject_idx).TMS.rMT = str2double(info{4});
        INFO(subject_idx).TMS.intensity = ceil(INFO(subject_idx).TMS.rMT * 1.15);

        % Save
        save(output_file, 'INFO','-append')
    end

    % Get recording info
%     if ~isfield(INFO(subject_idx).EEG, 'SR') || isempty(INFO(subject_idx).EEG.SR)
      if ~isfield(INFO(subject_idx), 'EEG') || ~isfield(INFO(subject_idx).EEG, 'SR') || isempty(INFO(subject_idx).EEG.SR)

        % Fill in the info
        prompt = {'sampling rate (Hz):', 'reference:', 'ground:', 'triggers:', 'blocks:'};
        dlgtitle = 'recording information';
        dims = [1 50];
        definput = {'20000', 'FCz', 'AFz', '1 - start, 2 - stimulation, 4 - preparatory, 8 - imperative', '1,2,3,4,5,6,7,8'};
%         recording_info = inputdlg(prompt,dlgtitle,dims,definput);
          recording_info = inputdlg(prompt, dlgtitle, dims, definput);

        INFO(subject_idx).EEG.SR = str2num(recording_info{1});
        INFO(subject_idx).EEG.ref = recording_info{2};
        INFO(subject_idx).EEG.ground = recording_info{3};
        triggers = split(recording_info{4}, ',');
        for t = 1:length(triggers)
            INFO(subject_idx).EEG.triggers(t).trigger = str2double(regexp(triggers{t}, '\d+\.?\d*', 'match', 'once'));
            pattern = sprintf('- ([a-zA-Z]+)');
            INFO(subject_idx).EEG.triggers(t).label = regexp(triggers{t}, pattern, 'tokens', 'once');
        end
        INFO(subject_idx).EEG.blocks = str2num(recording_info{5});

        % Save
        save(output_file, 'INFO','-append')
    end  
end
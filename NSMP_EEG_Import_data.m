function [DATASET, INFO] = NSMP_EEG_Import_data(FOLDER, INFO, subject_idx)
    % Study parameters
    study = 'NSMP';
    output_file = fullfile(FOLDER.output, [study '_output.mat']);

    % Define preprocessing parameters
    prompt = {'Crop margin (s):' 'Downsample ratio:'};
    dlgtitle = 'Continuous EEG preprocessing';
    dims = [1 60];
    definput = {'5' '20'};
    input = inputdlg(prompt, dlgtitle, dims, definput);
    params.suffix = {'crop' 'ds' 'notch' 'dc'};
    params.crop_margin = str2num(input{1});
    params.downsample = str2num(input{2});

    % Add letswave 6 to the top of search path
    addpath(genpath([FOLDER.toolbox '\letswave6-master']));

    % Identify import folder and subfolders
    session_folders = dir(sprintf('%s\\%s', FOLDER.raw, INFO(subject_idx).ID));
    session_date = datetime(INFO(subject_idx).session, 'InputFormat', 'dd-MMM-yyyy', 'Format', 'yyyy-MM-dd');
    session_date = char(session_date);
    for a = 1:length(session_folders)
        if contains(session_folders(a).name, session_date)
            params.folder = sprintf('%s\\%s', session_folders(a).folder, session_folders(a).name);
        end
    end
    params.blocks = INFO(subject_idx).EEG.blocks;

    % Check if identified subfolders are available 
    fprintf('Subject %d (%s): ', subject_idx, INFO(subject_idx).ID)
    data2import = dir(params.folder);
    if ~isempty(data2import)
        % Remove all datasets that are not labelled with numbers
        data_idx = logical([]);
        for b = 1:length(data2import)
            if isempty(str2num(data2import(b).name))
                data_idx(b) = true;
            else
                data2import(b).label = str2num(data2import(b).name);
                data_idx(b) = false;
            end
        end
        data2import(data_idx) = [];
        [~, file_idx] = sort([data2import.label]);
        data2import = data2import(file_idx);
    
        % Check for correct labels
        file_idx = false(1, length(data2import));
        for b = 1:length(data2import)
            for c = 1:length(params.blocks)
                if data2import(b).label == params.blocks(c)
                    file_idx(b) = true;
                end
            end
        end
        data2import = data2import(file_idx);
    
        % Check number of datasets
        fprintf('%d datasets found in the directory.\n', length(data2import))
        if length(data2import) ~= length(params.blocks)
            error('ERROR: This does not match with expected number of datasets (%d)\n.Please verify manually.\n', length(params.blocks))
        end
    else
        error('ERROR: no datasets found in the directory.\n')
    end

    % Load datasets
    fprintf('Loading datasets:\n')
    for d = 1:length(data2import)
        % provide update
        fprintf('--> block %d - ', d)
    
        % create the name
        params.name = sprintf('%s b%d', INFO(subject_idx).ID, d);
        
        % encode to the info structure
        if d == 1
            INFO(subject_idx).EEG.dataset(1).block = d;
        else
            INFO(subject_idx).EEG.dataset(end + 1).block = d;
        end
        INFO(subject_idx).EEG.dataset(end).subfolder = data2import(d).label;
        INFO(subject_idx).EEG.dataset(end).name = params.name;  
    
        % import the dataset
        DATASET.raw(d).block = d;
        [DATASET.raw(d).header, DATASET.raw(d).data, ~] = RLW_import_MEGA(data2import(d).folder, data2import(d).label);
    
        % rename in the header
        DATASET.raw(d).header.name = params.name;
    end  
    fprintf('Done.\n')

    % Check trigger and event numbers, ask for continuation
    for a = 1:length(DATASET.raw)
        % Check total number of triggers
        fprintf('block %d:\n%d triggers found\n', a, length(DATASET.raw(a).header.events));
        
        % Check trigger labels
        triggers = unique([DATASET.raw(a).header.events.code]');
        fprintf('triggers are labeled: ');
        for b = 1:length(triggers)
            fprintf('%s ', triggers(b))
        end
        fprintf('\n')
    
        % Check total number of events
        if ismember('1', triggers)
            events = 0;
            for b = 1:length(DATASET.raw(a).header.events)
                if str2double(DATASET.raw(a).header.events(b).code) == 1
                    events = events + 1;
                end
            end
        else
            error('ERROR: no trigger is labeled ''1''! Please check manually.')
        end
        fprintf('in total %d events found\n\n', events);
        pause(1)
    end

    % Add letswave 7 to the top of search path
    addpath(genpath([FOLDER.toolbox '\letswave7-master']));
    
    % Pre-process continuous data and save for letswave
    fprintf('Pre-processing:\n')
    for d = 1:length(DATASET.raw)
        % Provide update
        fprintf('--> %s\n', INFO(subject_idx).EEG.dataset(d).name)
    
        % Select data
        lwdata.header = DATASET.raw(d).header;
        lwdata.data = DATASET.raw(d).data; 
    
        % Assign electrode coordinates
        fprintf('Assigning electrode coordinates...\n')
        option = struct('filepath', sprintf('%s\\letswave7-master\\res\\electrodes\\spherical_locations\\Standard-10-20-Cap81.locs', FOLDER.toolbox), ...
            'suffix', '', 'is_save', 0);
        lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
        if d == 1
            INFO(subject_idx).EEG.processing(1).process = 'electrode coordinates assigned';
            INFO(subject_idx).EEG.processing(1).params.layout = 'standard 10-20-cap81';
            INFO(subject_idx).EEG.processing(1).suffix = [];
            INFO(subject_idx).EEG.processing(1).date = sprintf('%s', date);
        end
    
        % Re-label and count events
        fprintf('Checking events...\n') 
        event_idx = logical([]);
        for a = 1:length(lwdata.header.events)
            if isempty(str2num(lwdata.header.events(a).code))
                event_idx(a) = true;
            else
                event_idx(a) = false;
            end
        end
        lwdata.header.events(event_idx) = [];
        params.eventcodes = unique({lwdata.header.events.code});
        if length(params.eventcodes) == length(INFO(subject_idx).EEG.triggers)
            event_count = zeros(length(INFO(subject_idx).EEG.triggers), 1);
            for e = 1:length(lwdata.header.events)
                for a = 1:length(INFO(subject_idx).EEG.triggers)
                    if strcmp(lwdata.header.events(e).code, num2str(INFO(subject_idx).EEG.triggers(a).trigger))
                        lwdata.header.events(e).code = INFO(subject_idx).EEG.triggers(a).label{1};
                        event_count(a) = event_count(a) + 1;
                    end
                end
            end
        else
            error('ERROR: wrong number of triggers (%d) was found in the dataset!', length(params.eventcodes))
        end
        
    


        % Update & encode
        fprintf('%d events in total were found in the dataset:\n%s - %d events\n%s - %d events\n%s - %d events\n%s - %d events\n', ...
            length(lwdata.header.events), ...
            INFO(subject_idx).EEG.triggers(1).label{1}, event_count(1), ...
            INFO(subject_idx).EEG.triggers(2).label{1}, event_count(2), ...
            INFO(subject_idx).EEG.triggers(3).label{1}, event_count(3), ...
            INFO(subject_idx).EEG.triggers(4).label{1}, event_count(4));
        INFO(subject_idx).EEG.dataset(d).trials = event_count(1);
    
        % Crop data
        params.crop(1) = lwdata.header.events(1).latency - params.crop_margin;
        params.crop(2) = lwdata.header.events(end).latency + params.crop_margin;
        fprintf('Cropping ...\n')
        option = struct('xcrop_chk', 1, 'xstart', params.crop(1), 'xend', params.crop(2), ...
            'suffix', params.suffix{1}, 'is_save', 0);
        lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);
        if d == 1
            INFO(subject_idx).EEG.processing(2).process = 'continuous data cropped';
            INFO(subject_idx).EEG.processing(2).params.start = params.crop(1);
            INFO(subject_idx).EEG.processing(2).params.end = params.crop(2);
            INFO(subject_idx).EEG.processing(2).params.margin = params.crop_margin;
            INFO(subject_idx).EEG.processing(2).suffix = params.suffix{1};
            INFO(subject_idx).EEG.processing(2).date = sprintf('%s', date);
        end
    
        % Downsample 
        fprintf('Downsampling...\n')
        option = struct('x_dsratio', params.downsample, 'suffix', params.suffix{2}, 'is_save', 0);
        lwdata = FLW_downsample.get_lwdata(lwdata, option);
        if d == 1
            INFO(subject_idx).EEG.processing(3).process = sprintf('downsampled');
            INFO(subject_idx).EEG.processing(3).params.ratio = params.downsample;
            INFO(subject_idx).EEG.processing(3).params.fs_orig = 1/lwdata.header.xstep * params.downsample;
            INFO(subject_idx).EEG.processing(3).params.fs_final = 1/lwdata.header.xstep;
            INFO(subject_idx).EEG.processing(3).suffix = params.suffix{2};
            INFO(subject_idx).EEG.processing(3).date = sprintf('%s', date);
        end
    


        % Notch filter (50 Hz line noise)
        fprintf('Applying notch filter (50 Hz)...\n')
        option = struct();
        option.filter_type  = 'bandstop';
        option.low_cutoff   = 49;
        option.high_cutoff  = 51;
        option.filter_order = 4;
        option.suffix       = params.suffix{3};
        option.is_save      = 0;
        lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
        
        if d == 1
            INFO(subject_idx).EEG.processing(4).process = 'notch filtered';
            INFO(subject_idx).EEG.processing(4).params.type = 'bandstop';
            INFO(subject_idx).EEG.processing(4).params.limits = [49 51];
            INFO(subject_idx).EEG.processing(4).params.order = 4;
            INFO(subject_idx).EEG.processing(4).suffix = params.suffix{3};
            INFO(subject_idx).EEG.processing(4).date = sprintf('%s', date);
        end

        % Remove DC + linear detrend continuous data
        fprintf('Removing DC and applying linear detrend...\n')
        option = struct('linear_detrend', 1, 'suffix', params.suffix{4}, 'is_save', 0);
        lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
        if d == 1
            INFO(subject_idx).EEG.processing(5).process = sprintf('DC + linear detrend on continuous data');
            INFO(subject_idx).EEG.processing(5).suffix = params.suffix{4};
            INFO(subject_idx).EEG.processing(5).date = sprintf('%s', date);
        end
        fprintf('\n')

        % Save in letswave format to Data folder
        header = lwdata.header;
        data = lwdata.data;
        save(sprintf('%s\\%s.lw6', FOLDER.processed, header.name), 'header')
        save(sprintf('%s\\%s.mat', FOLDER.processed, header.name), 'data')
    
        % Update DATASET
        DATASET.raw(d).header = lwdata.header;
        DATASET.raw(d).data = lwdata.data; 
    end
    fprintf('Done.\n')

    % Save INFO
    save(output_file, 'INFO','-append')
end
function DATASET = NSMP_EEG_ReloadProcessed(FOLDER, DATASET, suffix)

fprintf('Reloading processed datasets (*%s*)...\n', suffix);

files = dir(fullfile(FOLDER.processed, ['*' suffix '*.lw6']));
assert(~isempty(files), 'No %s datasets found in processed folder', suffix);

DATASET.raw = struct([]);

for i = 1:numel(files)
    load(fullfile(files(i).folder, files(i).name), '-mat');   % loads header
    matfile = strrep(files(i).name, '.lw6', '.mat');
    load(fullfile(files(i).folder, matfile), 'data');

    DATASET.raw(i).header = header;
    DATASET.raw(i).data   = data;
end

fprintf('Reloaded %d rereferenced datasets.\n', numel(files));
end

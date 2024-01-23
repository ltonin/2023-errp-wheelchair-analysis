clearvars; clc;

subject = 'e2';

includepat  = {subject, 'discrete'};
excludepat  = {};
depthlevel  = 2;

rootpath      = '/mnt/data/Research/';
folder        = '2023_errp_wheelchair';
gdfpath       = [rootpath '/' folder '/'];

%% Processing parameters
eegchannels  = {'FP1', 'FP2', 'F1', 'FZ', 'F2', 'FC1', 'FCz', 'FC2', 'C1', 'CZ', 'C2', 'CP1', 'CP2'};
eogchannels  = {'EOG'};
filter.order = 4;
filter.bands = [1 45];
savedir      = ['analysis/amica/preprocessed/' num2str(length(eegchannels)) '/'];

%% Get datafiles
files = util_getfile3(gdfpath, '.gdf', 'include', includepat, 'exclude', excludepat, 'level', depthlevel);

nfiles = length(files);
if(nfiles > 0)
    util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
else
    error(['[io] - No files found with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
end

%% Create directory
util_mkdir(pwd, savedir);

%% Pre-processing files

for fId = 1:nfiles

    cfullname = files{fId};
    [cfilepath, cfilename, cfileext] = fileparts(cfullname);
    
    util_bdisp(['[io] + Loading file ' num2str(fId) '/' num2str(nfiles)]);
    disp(['     |- File: ' cfullname]);

    %% Loading data
    try
        [s, h] = sload(cfullname); 
    catch ME
        error(['[error] - Cannot load filename. Skipping it. ' Me.message]);
    end
    
    %% Extracting selected EEG channels
    util_bdisp('[proc] + Extracting selected EEG channels (in this ordered):');
    disp(['       |- Channels: ' strjoin(eegchannels, ', ')]);
   
    [~, eegchannelidx] = ismember(lower(eegchannels), lower(h.Label));
    s_eeg = s(:, eegchannelidx);
    
    %% Extracting selected EOG channels
    util_bdisp('[proc] + Extracting selected EOG channels (in this ordered):');
    disp(['       |- Channels: ' strjoin(eogchannels, ', ')]);
   
    [~, eogchannelidx] = ismember(lower(eogchannels), lower(h.Label));
    s_eog = s(:, eogchannelidx);

    %% Processing data
    util_bdisp('[proc] + Processing EEG the data');

    % Compute spatial filter
    disp('       |- Spatial filter: CAR');
    s_spatial = proc_car(s_eeg);

    % Compute bandpass filters
    disp(['       |- Bandpass filter order ' num2str(filter.order) ': [' strjoin(compose('%g', filter.bands), ' ') '] Hz']);
    s_bp = filt_bp(s_spatial, filter.order, filter.bands, h.SampleRate);

    %% Storing data and settings
    eeg = s_bp;
    eog = s_eog;

    settings.preprocess.spatial = 'car';
    settings.preprocess.filter  = filter;
    settings.channels.eeg       = eegchannels;
    settings.channels.eog       = eogchannels;
    settings.events.TYP         = h.EVENT.TYP;
    settings.events.POS         = h.EVENT.POS;
    settings.events.DUR         = h.EVENT.DUR;
    settings.samplerate         = h.SampleRate;

    %% Saving preprocessed data
    sfilename = fullfile(savedir, [cfilename '.mat']);
    util_bdisp(['[out] + Saving AMICA preprocessed data in: ' sfilename]);
    save(sfilename, 'eeg', 'eog', 'settings'); 
end
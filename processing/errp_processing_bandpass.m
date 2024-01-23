clearvars; clc;

subject = 'd6';

includepat  = {subject, 'discrete'};
excludepat  = {};
depthlevel  = 2;

rootpath      = '/mnt/data/Research/';
folder        = '2023_errp_wheelchair';
gdfpath       = [rootpath '/' folder '/'];
chanlocspath  = 'chanlocs64.mat';

%% Processing parameters
eegchannels       = {'FP1', 'FP2', 'F1', 'FZ', 'F2', 'FC1', 'FCz', 'FC2', 'C1', 'CZ', 'C2', 'CP1', 'CP2'};
eogchannels       = {'EOG'};
chanlocstr        = load(chanlocspath);
chanlocs          = errp_util_get_chanlocs(eegchannels, chanlocstr.chanlocs);
spatial.filter    = 'car'; 
filter.order      = 3;
filter.bands      = [2 8];
savedir           = ['analysis/none/' spatial.filter '/bandpass/' num2str(length(eegchannels)) '/'];

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

%% Processing files
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
    
    
    %% Get information from filename
    cinfo = errp_util_get_info(cfullname);
    
    %% Processing data
    util_bdisp('[proc] + Processing the data');

    % Compute Spatial filter
    disp(['       |- Spatial filter: ' spatial.filter]);

    switch(spatial.filter)
        case 'none'
            s_eeg_spatial = s_eeg;
        case 'car'
            s_eeg_spatial = proc_car(s_eeg, 'excluded', [1 2]);
        case 'laplacian'
            error('Laplacian spatial filter not implemented yet');
        otherwise
            error(['Unknown spatial filter selected ' spatial.filter]);
    end
    
    % Compute bandpass filters
    s_eeg_bp = filt_bp(s_eeg_spatial, filter.order, filter.bands, h.SampleRate);
    s_eog_bp = filt_bp(s_eog, filter.order, filter.bands, h.SampleRate);
    keyboard
    eeg = s_eeg_bp;
    eog = s_eog_bp;
    
    % Extracting events
    disp('       |-Extract events');
    cevents     = h.EVENT;
    events.TYP = cevents.TYP;
    events.POS = cevents.POS;
    events.DUR = cevents.DUR;
    
    % Task
    disp('       |-Extract task info');
    task = cinfo.task;
    
    % Extra
    disp('       |-Extract extra info');
    device  = cinfo.extra1;
    control = cinfo.extra2;

    
    %% Create settings structure
    settings.spatial           = spatial;
    settings.filter            = filter;
    settings.channels.eeg      = eegchannels;
    settings.channels.eog      = eogchannels;
    settings.channels.chanlocs = chanlocs;
    settings.events.TYP        = h.EVENT.TYP;
    settings.events.POS        = h.EVENT.POS;
    settings.events.DUR        = h.EVENT.DUR;
    settings.samplerate        = h.SampleRate;
    settings.task.name         = task;
    settings.device.name       = device;
    settings.control.name      = control;
    settings.control.legend    = {'discrete', 'continuous', 'unknown'};
    settings.info              = cinfo;
    
    sfilename = fullfile(savedir, [cfilename '.mat']);
    util_bdisp(['[out] - Saving bandpass in: ' sfilename]);
    save(sfilename, 'eeg', 'eog', 'events', 'settings'); 
   
end


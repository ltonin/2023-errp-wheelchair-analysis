clearvars; clc;

subject = 'e10';

includepat  = {subject};
excludepat  = {};
depthlevel  = 2;

rootpath    = '/mnt/data/Research/';
folder      = '2023_errp_wheelchair';
gdfpath     = [rootpath '/' folder '/'];

artifactrej       = 'none';
spatialfilter     = 'car';
savedir           = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];
recompute         = true;
lbchannels        = {'FP1', 'FP2', 'FZ', 'FC5', 'FC1', 'FC2', 'FC6', 'C3', 'CZ', 'C4', 'CP5', 'CP1', 'CP2', 'CP6', 'P3',  'Pz',  'P4', ...
                 'EOG', 'F1',  'F2', 'FC3', 'FCZ', 'FC4',  'C5',  'C1',  'C2',  'C6', 'CP3', 'CP4',  'P5',  'P1',  'P2',  'P6'};
eogchannels       = {'EOG'};
eogchannelidx     = find(ismember(lower(lbchannels), lower(eogchannels)));
eegchannels       = setdiff(lbchannels, eogchannels, 'stable');
eegchannelidx     = find(ismember(lower(lbchannels), lower(eegchannels)));


%% Processing parameters
nchannels          = length(eegchannels);
layout_channels    = 'eeg.antneuro.32.noeog.errp_mi';
montage_channels   = proc_get_montage(layout_channels);
included_car_chans = {'FCz', 'Cz', 'Pz', 'FC1', 'FC2', 'C1', 'C2', 'CP1', 'CP2', 'P1', 'P2'};     % Included channels for CAR
included_car_index = find(ismember(lower(eegchannels), lower(included_car_chans)));
laplacian_mask     = proc_laplacian_mask(montage_channels, nchannels);
bands              = [2 8];
filtorder          = 3;


%% Get datafiles
files = util_getfile3(gdfpath, '.gdf', 'include', includepat, 'exclude', excludepat, 'level', depthlevel);

NumFiles = length(files);
if(NumFiles > 0)
    util_bdisp(['[io] - Found ' num2str(NumFiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
else
    error(['[io] - No files found with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
end

%% Create/Check for savepath
util_mkdir(pwd, savedir);

%% Processing files
for fId = 1:NumFiles
    cfullname = files{fId};
    [cfilepath, cfilename, cfileext] = fileparts(cfullname);
    
    util_bdisp(['[io] + Loading file ' num2str(fId) '/' num2str(NumFiles)]);
    disp(['     |-File: ' cfullname]);
    
    %% Check if the file has been already processed
    [~, pfilename] = fileparts(cfullname);
    if (recompute == false) && exist([savedir pfilename '.mat'], 'file') == 2
        disp('     |-Processed bandpass already exists. Skipping the recomputing');
        continue;
    end
    
    %% Loading data
    disp('     |-Loading GDF data');
    try
        [s, h] = sload(cfullname);
        eeg = s(:, eegchannelidx);
        eog = s(:, eogchannelidx);
    catch ME
        warning('[warning] - Cannot load filename. Skipping it.');
        warning(['[warning] - Error: ' ME.message]);
        continue;
    end
    
    
    %% Get information from filename
    cinfo = errp_util_get_info(cfullname);
    
    %% Processing data
    util_bdisp('[proc] + Processing the data');

    % Compute Spatial filter
    disp(['       |-Spatial filter: ' spatialfilter ' (' layout_channels ')']);

    switch(spatialfilter)
        case 'none'
            eeg_filt = eeg;
        case 'car'
            eeg_filt = proc_car(eeg, 'included', included_car_index);
        case 'laplacian'
            eeg_filt = eeg*laplacian_mask;
        otherwise
            error(['Unknown spatial filter selected ' spatialfilter]);
    end
    
    % Compute bandpass filters
    eeg_bp = filt_bp(eeg_filt, filtorder, bands, h.SampleRate);
    eog_bp = filt_bp(eog, filtorder, bands, h.SampleRate);

    P = eeg_bp;
    E = eog_bp;
    
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
    settings.data.filename          = cfullname;
    settings.data.nsamples          = size(eeg, 1);
    settings.data.nchannels         = size(eeg, 2);
    settings.data.lchannels         = eegchannels;
    settings.data.samplerate        = h.SampleRate;
    settings.artifact.name          = artifactrej;
    settings.spatial.laplacian      = laplacian_mask;
    settings.spatial.included_car   = included_car_index;
    settings.spatial.filter         = spatialfilter;
    settings.bandpass.order         = filtorder;
    settings.bandpass.bands         = bands;
    settings.task.name              = task;
    settings.device.name            = device;
    settings.control.name           = control;
    settings.control.legend         = {'discrete', 'continuous', 'unknown'};
    settings.info                   = cinfo;
    
    
    sfilename = fullfile(savedir, [pfilename '.mat']);
    util_bdisp(['[out] - Saving bandpass in: ' sfilename]);
    save(sfilename, 'P', 'E', 'events', 'settings'); 

    
end


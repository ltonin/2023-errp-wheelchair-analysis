clearvars; clc;

subject = 'e6';

includepat  = {subject, 'discrete'};
excludepat  = {};
depthlevel  = 2;

artifactrej   = 'amica';
spatialfilter = 'car';
datapath      = ['analysis/' artifactrej '/' spatialfilter '/cleaned/'];
savedir       = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];
eeg_channels  = {'FP1', 'FP2', 'FZ', 'FC5', 'FC1', 'FC2', 'FC6', 'C3', 'CZ', 'C4', 'CP5', 'CP1', 'CP2', 'CP6', 'P3',  'Pz',  'P4', ...
                 'F1',  'F2', 'FC3', 'FCZ', 'FC4',  'C5',  'C1',  'C2',  'C6', 'CP3', 'CP4',  'P5',  'P1',  'P2',  'P6'};
eog_channels  = {'EOG'};
filtorder     = 3;
bands         = [2 8];

%% Get datafiles
files = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);

nfiles = length(files);
if(nfiles > 0)
    util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
else
    error(['[io] - No files found with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
end

%% Create directory
util_mkdir(pwd, savedir);

%% Get data files
files  = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
nfiles = length(files);

for fId = 1:nfiles

    cfullname = files{fId};
    [cfilepath, cfilename, cfileext] = fileparts(cfullname);

    util_bdisp(['[io] + Loading cleaned data ' num2str(fId) '/' num2str(nfiles)]);
    disp(['     |-File: ' cfullname]);
    
    data = load(cfullname);

    s          = double(data.s);
    eog        = double(data.eog);
    events     = data.events;
    settings   = data.settings;
    samplerate = data.settings.data.samplerate;

    %% Processing data
    util_bdisp('[proc] + Processing the data');

    disp(['       |-Bandpass filter order ' num2str(filtorder) ': [' num2str(bands(1)) ' ' num2str(bands(2)) '] Hz']);
   
    
    % Compute bandpass filters
    s_bp   = filt_bp(s, filtorder, bands, samplerate);
    eog_bp = filt_bp(eog, filtorder, bands, samplerate);
    
    P = s_bp;
    E = eog_bp;

    settings.bandpass.order = filtorder;
    settings.bandpass.bands = bands;

    sfilename = fullfile(savedir, [cfilename '.mat']);
    util_bdisp(['[out] - Saving bandpass in: ' sfilename]);
    save(sfilename, 'P', 'E', 'events', 'settings'); 

end
clearvars; clc;

subject = 'e6';

includepat  = {subject, 'discrete'};
excludepat  = {};
depthlevel  = 2;

artifactrej   = 'amica';
spatialfilter = 'car';
datapath      = ['analysis/' artifactrej '/' spatialfilter '/raw/'];
savedir       = ['analysis/' artifactrej '/' spatialfilter '/cleaned/'];
eeg_channels  = {'FP1', 'FP2', 'FZ', 'FC5', 'FC1', 'FC2', 'FC6', 'C3', 'CZ', 'C4', 'CP5', 'CP1', 'CP2', 'CP6', 'P3',  'Pz',  'P4', ...
                 'F1',  'F2', 'FC3', 'FCZ', 'FC4',  'C5',  'C1',  'C2',  'C6', 'CP3', 'CP4',  'P5',  'P1',  'P2',  'P6'};
eog_channels  = {'EOG'};

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

%% AMICA component selection
files  = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
nfiles = length(files);

for fId = 1:nfiles

    cfullname = files{fId};
    [cfilepath, cfilename, cfileext] = fileparts(cfullname);

    util_bdisp(['[io] + Loading amica components ' num2str(fId) '/' num2str(nfiles)]);
    disp(['     |-File: ' cfullname]);
    
    amica = load(cfullname);

    util_bdisp('[amica] + Running component selection');
    eeg        = amica.eeg;
    eeg.icaact = geticaact(eeg);
    eeg        = iclabel(eeg);
    [~, rmcomps] = errp_pop_viewprops(eeg, 0);
    
    util_bdisp('[amica] + Removing selected components');
    eeg = pop_subcomp(eeg, rmcomps);
    eeg.etc.removed_comps = rmcomps;
    
    s   = eeg.data';
    eog = amica.eog;

    %% Get information from filename
    util_bdisp('[proc] + Getting file information');
    cinfo = errp_util_get_info(cfullname);

    % Extracting events
    disp('       |-Extract events');
    cevents    = amica.h.EVENT;
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
    settings.data.nsamples          = size(s, 1);
    settings.data.nchannels         = size(s, 2);
    settings.data.lchannels         = eeg_channels;
    settings.data.samplerate        = amica.h.SampleRate;
    settings.artifact.name          = artifactrej;
    settings.spatial.filter         = spatialfilter;
    settings.task.name              = task;
    settings.device.name            = device;
    settings.control.name           = control;
    settings.control.legend         = {'discrete', 'continuous', 'unknown'};
    settings.info                   = cinfo;

    %% Save cleaned data
    sfilename = fullfile(savedir, [cfilename '.mat']);
    util_bdisp(['[out] - Saving cleaned data in: ' sfilename]);
    save(sfilename, 's', 'eog', 'events', 'settings'); 
end

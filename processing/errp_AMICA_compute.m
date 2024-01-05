clearvars; clc;

subject = 'e6';

includepat  = {subject, 'discrete'};
excludepat  = {};
depthlevel  = 2;

rootpath    = '/mnt/data/Research/';
folder      = '2023_errp_wheelchair';
gdfpath     = [rootpath '/' folder '/'];

artifactrej   = 'amica';
spatialfilter = 'car';
savedir       = ['analysis/' artifactrej '/' spatialfilter '/raw/'];
lbchannels    = {'FP1', 'FP2', 'FZ', 'FC5', 'FC1', 'FC2', 'FC6', 'C3', 'CZ', 'C4', 'CP5', 'CP1', 'CP2', 'CP6', 'P3',  'Pz',  'P4', ...
                 'EOG', 'F1',  'F2', 'FC3', 'FCZ', 'FC4',  'C5',  'C1',  'C2',  'C6', 'CP3', 'CP4',  'P5',  'P1',  'P2',  'P6'};

%% Processing parameters
eogchannel    = {'EOG'};
eogchannelidx = find(ismember(lower(lbchannels), lower(eogchannel)));
eegchannels   = setdiff(lbchannels, eogchannel, 'stable');
eegchannelidx = find(ismember(lower(lbchannels), lower(eegchannels)));

amica_filterorder  = 4;
amica_filterbands  = [2 80];

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

%% Get generic 64 channel locations
load('chanlocs64.mat');

%% Compute AMICA components for each file

for fId = 1:nfiles

    cfullname = files{fId};
    [cfilepath, cfilename, cfileext] = fileparts(cfullname);
    
    util_bdisp(['[io] + Loading file ' num2str(fId) '/' num2str(nfiles)]);
    disp(['     |-File: ' cfullname]);

    %% Loading data
    disp('     |-Loading GDF data');
    try
        [s, h] = sload(cfullname);
        s_eeg = s(:, eegchannelidx);
        s_eog = s(:, eogchannelidx);
    catch ME
        warning('[warning] - Cannot load filename. Skipping it.');
        warning(['[warning] - Error: ' ME.message]);
        continue;
    end
    
    %% Processing data
    util_bdisp('[proc] + Processing the data');

    % Compute spatial filter
    disp(['       |-Spatial filter: ' spatialfilter]);

    switch(spatialfilter)
        case 'none'
            s_eeg_filt = s_eeg;
        case 'car'           
            s_eeg_filt = proc_car(s_eeg);
        otherwise
            error(['Unknown spatial filter selected ' spatialfilter]);
    end
    
    % Compute bandpass filters
    disp(['       |-Bandpass filter order ' num2str(amica_filterorder) ': [' num2str(amica_filterbands(1)) ' ' num2str(amica_filterbands(2)) '] Hz']);
    s_eeg_bp = filt_bp(s_eeg_filt, amica_filterorder, amica_filterbands, h.SampleRate);
    s_eog_bp = filt_bp(s_eog, amica_filterorder, amica_filterbands, h.SampleRate);

    
    %% Compute AMICA
    util_bdisp('[proc] + Computing AMICA');
    eeg          = eeg_emptyset;
    eeg.srate    = h.SampleRate;
    eeg.data     = s_eeg_bp';
    eeg.chanlocs = errp_util_get_chanlocs(eegchannels, chanlocs);
    eog          = s_eog_bp;

    amicadrk = getrank(eeg.data);
    disp(['       |-Suggested rank: ' num2str(amicadrk)]);

    % AMICA arguments
    amica_arglist;
    arglist(4)  = {size(eeg.data, 1)};       % nchannels
    arglist(6)  = {amicadrk};               % PCA keep 
    arglist(8)  = {8};                       % Max number of threads
    arglist(12) = {1};                       % Number of processes
    
    % Run AMICA
    [W,S,mods] = runamica15(eeg.data(:,:), arglist{:});
    eeg.icaweights = W;
    eeg.icasphere  = S(1:size(W,1),:);
    eeg.icawinv    = mods.A(:,:,1);

    % Check coherence of data structure
    eeg = eeg_checkset(eeg, 'ica');
    
    % Save AMICA components
    sfilename = fullfile(savedir, [cfilename '.mat']);
    util_bdisp(['[out] - Saving AMICA components in: ' sfilename]);
    save(sfilename, 'eeg', 'eog', 'h'); 
end
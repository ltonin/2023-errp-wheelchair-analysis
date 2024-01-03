clearvars; clc;

subject = 'd6';

includepat  = {subject, 'discrete'};
excludepat  = {};
depthlevel  = 2;

rootpath    = '/mnt/data/Research/';
folder      = '2023_errp_wheelchair';
gdfpath     = [rootpath '/' folder '/'];

artifactrej   = 'amica';
spatialfilter = 'car';
savedir       = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];
lbchannels    = {'FP1', 'FP2', 'FZ', 'FC5', 'FC1', 'FC2', 'FC6', 'C3', 'CZ', 'C4', 'CP5', 'CP1', 'CP2', 'CP6', 'P3',  'Pz',  'P4', ...
                 'EOG', 'F1',  'F2', 'FC3', 'FCZ', 'FC4',  'C5',  'C1',  'C2',  'C6', 'CP3', 'CP4',  'P5',  'P1',  'P2',  'P6'};

%% Processing parameters
eogchannel    = {'EOG'};
eogchannelidx = find(ismember(lower(lbchannels), lower(eogchannel)));
eegchannels   = setdiff(lbchannels, eogchannel, 'stable');
eegchannelidx = find(ismember(lower(lbchannels), lower(eegchannels)));

amica_filterorder  = 4;
amica_filterbands  = [3 45];
bands              = [2 8];
filtorder          = 3;

timeborder = 3; % seconds


%% Get datafiles
files = util_getfile3(gdfpath, '.gdf', 'include', includepat, 'exclude', excludepat, 'level', depthlevel);

nfiles = length(files);
if(nfiles > 0)
    util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
else
    error(['[io] - No files found with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
end

%% Create/Check for savepath
util_mkdir(pwd, savedir);

%% Get data size info to allocate memory and speedup the concatenation
datasize     = get_gdf_datasize(files);
allnsamples  = sum(datasize(1, :));
allnchannels = unique(datasize(2, :));
fileseek = 1;

%% Processing files

s_eeg_bp_all = nan(allnsamples, length(eegchannelidx));
s_eog_bp_all = nan(allnsamples, 1);

for fId = 1:nfiles

    % Get current position 
    cstart   = fileseek;
    cstop    = cstart + datasize(1, fId) - 1;

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
    
    
    %% Get information from filename
    cinfo = errp_util_get_info(cfullname);
    
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

    % Remove borders
    disp(['       |-Remove borders (' num2str(timeborder) ' s)']) 
    samplesborder = floor(timeborder*h.SampleRate);
    s_eeg_bp(1:samplesborder, :)       = 0;
    s_eeg_bp(end-samplesborder:end, :) = 0;
    s_eog_bp(1:samplesborder, :)       = 0;
    s_eog_bp(end-samplesborder:end, :) = 0;

    s_eeg_bp_all(cstart:cstop, :) = s_eeg_bp;
    s_eog_bp_all(cstart:cstop, 1) = s_eog_bp;

    % Update the fileseek position
    fileseek = cstop + 1;

end

eeg = s_eeg_bp_all';

%%
% Compute AMICA
drk = getrank(eeg);
fprintf('[compute_ica] Suggested rank is %d \n', drk)

amica_arglist
arglist(4) = {size(eeg,1)};
arglist(6) = {drk};
[W,S,mods] = runamica15(eeg(:,:), arglist{:});
aEEG.data = eeg;
aEEG.icaweights = W;
aEEG.icasphere  = S(1:size(W,1),:);
aEEG.icawinv    = mods.A(:,:,1);

aEEG = eeg_checkset(eeg, 'ica');


%% IC selection

aEEG.icaact = geticaact(aEEG);
EEG = iclabel(EEG);
pop_viewprops(EEG, 0);

prompt = 'Enter the components to remove:\n';
rm_comps = (str2num(input(prompt, "s")))';

EEG = pop_subcomp(EEG, rm_comps);
EEG.etc.removed_comps = rm_comps;

data = EEG.data;

function dsizes = get_gdf_datasize(filepaths)

    nfiles = length(filepaths);
    ndimensions = 2;                            % samples x chans
    
    dsizes = zeros(ndimensions, nfiles);
    
    for fId = 1:nfiles
        cfilepath = filepaths{fId};
        h = sopen(cfilepath, 'r');
        dsizes(:, fId) = [h.NRec*h.SampleRate h.NS];    
    end

end

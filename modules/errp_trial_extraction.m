clearvars; clc;

subject = 'e5';
includepat  = {subject, 'discrete'};
excludepat  = {};
spatialfilter = 'laplacian';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath    = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];


datafiles = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
ndatafiles = length(datafiles);
util_bdisp(['[io] - Found ' num2str(ndatafiles) ' data files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

%% Event parameters
Stop        = 100;
CommandLx   = 101;
CommandFx   = 102;
CommandRx   = 103;
ErrorLx     = 5101;
ErrorRx     = 5103;
NoReleaseLx = 4101;
NoReleaseFx = 4102;
NoReleaseRx = 4103;

RefractoryTime = 3;
EpochTime      = [-0.5 2.0];


%% Import data
util_bdisp(['[io] - Importing ' num2str(ndatafiles) ' files from ' datapath ':']);
[P, events, labels, settings] = errp_concatenate_bandpass(datafiles);
nsamples  = size(P, 1);
nchannels = size(P, 2);
SampleRate = settings.data.samplerate;

F = P;

%% Preprocess event
% Check for unexptected events
unexpected = find(ismember(unique(events.TYP), [Stop CommandLx CommandRx CommandFx ErrorLx ErrorRx NoReleaseLx NoReleaseRx NoReleaseFx]) == 0);
if(isempty(unexpected) == false)
    for uId = 1:length(unexpected)
        warning(['Event ''' num2str(unexpected(uId)) ''' is in the event structure']);
    end
end

% Extract command events
cmdindex = events.TYP == CommandLx | events.TYP == CommandRx | events.TYP == ErrorLx | events.TYP == ErrorRx | events.TYP == NoReleaseLx | events.TYP == NoReleaseRx;
DUR = events.DUR(cmdindex);
TYP = events.TYP(cmdindex);
POS = events.POS(cmdindex);

%% Extract epochs
EpochSamples = floor(EpochTime*settings.data.samplerate);
sepoch  = sum(abs(EpochSamples));
ntrials = length(POS);
E = nan(sepoch, nchannels, ntrials);

for trId = 1:ntrials
    cstart = POS(trId) + EpochSamples(1);
    cstop  = POS(trId) + EpochSamples(2) - 1;
    E(:, :, trId) = F(cstart:cstop, :);
end


%% Create label vectors

% Not released
Nk = false(ntrials, 1);
noreleaseindex = TYP == NoReleaseLx | TYP == NoReleaseRx;
Nk(noreleaseindex) = true;

% Command correct
Ck = false(ntrials, 1);
correctindex = TYP == CommandLx | TYP == CommandRx;
Ck(correctindex) = true;

% Command error
Ek = false(ntrials, 1);
errorindex = TYP == ErrorLx | TYP == ErrorRx;
Ek(errorindex) = true;



%% Plotting
t = EpochTime(1):1/settings.data.samplerate:EpochTime(2) - 1/settings.data.samplerate;
refchannel = {'FCz'};

[~, refchannelidx] = ismember(upper(refchannel), settings.data.lchannels);
nrefchannel = length(refchannel);

index_err = Ek;
index_cor = Ck;
index_nor = Nk;
perc_err = 100*sum(index_err)./(sum(index_err) + sum(index_cor));

m_err = squeeze(mean(E(:, :, index_err), 3));
m_cor = squeeze(mean(E(:, :, index_cor), 3));
m_nor = squeeze(mean(E(:, :, index_nor), 3));

fig = figure;
nrows = 3;
ncols = 8;

subplot(nrows, ncols, [9 10 11 12]);
plot_errp(t, m_cor(:, refchannelidx), m_err(:, refchannelidx));
plot_vline(0, 'k');
plot_hline(0, 'k');
title(['subject: ' subject ' | channel: ' char(refchannel) ' | error vs. correct']);

subplot(nrows, ncols, [13 14 15 16]);
plot_errp(t, m_cor(:, refchannelidx), m_nor(:, refchannelidx));
plot_vline(0, 'k');
plot_hline(0, 'k');
title(['Subject: ' subject ' | channel: ' char(refchannel) ' | no-realease vs. correct']);

load('chanlocs64.mat');
topowins = [-0.5 -0.25; 0.25 0.5; 0.75 1.0; 1.25 1.5];
t_corr = nan(nchannels, size(topowins, 2));
t_err = nan(nchannels, size(topowins, 2));
t_nor = nan(nchannels, size(topowins, 2));
for tId = 1:size(topowins, 1)
    csupport = find(t >= topowins(tId, 1) & t <= topowins(tId, 2));
    t_corr(:, tId) = mean(m_cor(csupport, :));
    t_err(:, tId)  = mean(m_err(csupport, :));
    t_nor(:, tId)  = mean(m_nor(csupport, :));
end

maplimits = [-5 5];

[chanlocs32, index32] = plot_getchanlocs(settings.data.lchannels, chanlocs);
for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, tId);
    topoplot(t_corr(index32, tId), chanlocs32, 'maplimits', maplimits);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))])
end

for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, 16+tId);
    topoplot(t_err(index32, tId), chanlocs32, 'maplimits', maplimits);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))])
end

for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, 4+tId);
    topoplot(t_corr(index32, tId), chanlocs32, 'maplimits', maplimits);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))])
end

for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, 20+tId);
    topoplot(t_nor(index32, tId), chanlocs32, 'maplimits', maplimits);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))])
end




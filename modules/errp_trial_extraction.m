clearvars; clc;

subject = 'd6';
includepat  = {subject, 'discrete'};
excludepat  = {};
spatialfilter = 'car';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath      = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];

bagSubPath    = 'bags/aligned';
gdfSubPath    = [artifactrej '/' spatialfilter '/bandpass'];

datafiles     = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
bagfiles      = replace(datafiles, gdfSubPath, bagSubPath); 

ndatafiles = length(datafiles);
util_bdisp(['[io] - Found ' num2str(ndatafiles) ' data files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') ')']);

%% Parameters
Stop        = 100;
CommandLx   = 101;
CommandFx   = 102;
CommandRx   = 103;
ErrorLx     = 5101;
ErrorRx     = 5103;
NoReleaseLx = 4101;
NoReleaseFx = 4102;
NoReleaseRx = 4103;

EpochTime      = [-1.0 2.0];

%% Importing data
util_bdisp(['[io] - Importing ' num2str(ndatafiles) ' files from ' datapath ':']);
[eeg, eog, events, labels, settings] = errp_concatenate_bandpass(datafiles);
[pose, twist, cmdvel, bagevents, baglabels] = errp_concatenate_bag(bagfiles);

nsamples  = size(eeg, 1);
nchannels = size(eeg, 2);
SampleRate = settings.data.samplerate;

%% Extract trials
util_bdisp('[proc] + Extracting trials:');

% Extract events
disp('     |-Extract events');
evtidx = errp_util_get_event_type([CommandLx CommandRx ErrorLx ErrorRx NoReleaseFx NoReleaseRx NoReleaseLx], events.TYP);
DUR = events.DUR(evtidx);
TYP = events.TYP(evtidx);
POS = events.POS(evtidx);

bagevtidx = errp_util_get_event_type([CommandLx CommandRx ErrorLx ErrorRx NoReleaseFx NoReleaseRx NoReleaseLx], bagevents.TYP);
bagERR = bagevents.ERR(bagevtidx);
bagTYP = bagevents.TYP(bagevtidx);
bagPOS = bagevents.POS(bagevtidx);

bagDEL = errp_compute_delay(twist.z, bagPOS);

plot_delay(twist.z, POS, bagDEL)

POS = POS + bagDEL;

% Extract epochs
disp(['     |-Extract epochs [' num2str(EpochTime(1)) ' ' num2str(EpochTime(2)) '] seconds']);
EpochSamples = floor(EpochTime*settings.data.samplerate);
sepoch       = sum(abs(EpochSamples));
ntrials      = length(POS);
T = nan(sepoch, nchannels, ntrials);
E = nan(sepoch, 1, ntrials);

for trId = 1:ntrials
    cstart = POS(trId) + EpochSamples(1);
    cstop  = POS(trId) + EpochSamples(2) - 1;
    T(:, :, trId) = eeg(cstart:cstop, :);
    E(:, :, trId) = eog(cstart:cstop, :);
end

% Create label vectors
disp('     |-Create label vectors');

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

chanlocs64 = load('chanlocs64.mat');
chanlocs32 = errp_util_get_chanlocs(settings.data.lchannels, chanlocs64.chanlocs);

perc_err = 100*sum(Ek)./(sum(Ek) + sum(Ck));

m_eeg_err = squeeze(mean(T(:, :, Ek), 3));
m_eeg_cor = squeeze(mean(T(:, :, Ck), 3));
m_eog_err = squeeze(mean(E(:, :, Ek), 3));
m_eog_cor = squeeze(mean(E(:, :, Ck), 3));

topowins = [-0.5 -0.25; 0.25 0.5; 0.75 1.0; 1.25 1.5];
maplimits = [-2 2];

% Figure
fig = figure;
nrows = 4;
ncols = 8;

% Topoplot correct
for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, tId);
    cstart = find(t >= topowins(tId, 1), 1, 'first');
    cstop  = find(t <= topowins(tId, 2), 1, 'last');
    topoplot(mean(m_eeg_cor(cstart:cstop, :), 1), chanlocs32, 'maplimits', maplimits);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))])
end

% EOG signal
subplot(nrows, ncols, [9 10 11 12]);
hold on;
plot(t, m_eog_cor(:, 1), 'b');
plot(t, m_eog_err(:, 1), 'r');
hold off;
grid on;
plot_vline(0, 'k');
plot_hline(0, 'k');
title(['subject: ' subject ' | channel: EOG | error vs. correct']);

% EEG signal
subplot(nrows, ncols, [17 18 19 20]);
plot_errp(t, m_eeg_cor(:, refchannelidx), m_eeg_err(:, refchannelidx));
plot_vline(0, 'k');
plot_hline(0, 'k');
title(['Subject: ' subject ' | channel: ' char(refchannel) ' | no-realease vs. correct']);

% Topoplot error
for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, tId + 24);
    cstart = find(t >= topowins(tId, 1), 1, 'first');
    cstop  = find(t <= topowins(tId, 2), 1, 'last');
    topoplot(mean(m_eeg_err(cstart:cstop, :), 1), chanlocs32, 'maplimits', maplimits);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))])
end

% Imagesc error trials

subplot(nrows, ncols, [5 6 7 8 13 14 15 16]);
imagesc(t, 1:sum(Ek), squeeze(T(:, refchannelidx, Ek))', [-15 15]);
plot_vline(0, 'k');
title(['Subject: ' subject ' | channel: ' char(refchannel) ' | error trials']);
xlabel('time [s]')
ylabel('# trial')


subplot(nrows, ncols, [21 22 23 24 29 30 31 32]);
imagesc(t, 1:sum(Ck), squeeze(T(:, refchannelidx, Ck))', [-15 15]);
plot_vline(0, 'k');
title(['Subject: ' subject ' | channel: ' char(refchannel) ' | correct trials']);
xlabel('time [s]')
ylabel('# trial')




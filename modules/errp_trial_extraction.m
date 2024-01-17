clearvars; clc;

subject = 'd6';
includepat  = {subject, 'discrete'};
excludepat  = {};
spatialfilter = 'car';
artifactrej   = 'none'; % {'FORCe', 'none'}
datapath      = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];
bagspath      = 'analysis/bags/aligned/';

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

time_no_release = round(512 * 2.8);
compensation = 3; % 1: delay on the user perception, 3 delay on the "almost peak"

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

evtidx = errp_reduce_events(events, evtidx, time_no_release);

DUR = events.DUR(evtidx);
TYP = events.TYP(evtidx);
POS = events.POS(evtidx);

bagevtidx = errp_util_get_event_type([CommandLx CommandRx ErrorLx ErrorRx NoReleaseFx NoReleaseRx NoReleaseLx], bagevents.TYP);

bagevtidx = errp_reduce_events(bagevents, bagevtidx, time_no_release);

bagERR = bagevents.ERR(bagevtidx);
bagTYP = bagevents.TYP(bagevtidx);
bagPOS = bagevents.POS(bagevtidx);

% Compute and add the delay
bagDEL = errp_compute_delay(twist.z, bagPOS, compensation);

plot_delay(twist.z, bagPOS, bagDEL, true)

POS    = POS + bagDEL;
bagPOS = bagPOS + bagDEL;

% Extract epochs
disp(['     |-Extract epochs [' num2str(EpochTime(1)) ' ' num2str(EpochTime(2)) '] seconds']);
EpochSamples = floor(EpochTime*settings.data.samplerate);
sepoch       = sum(abs(EpochSamples));
ntrials      = length(POS);
T = nan(sepoch, nchannels, ntrials);
E = nan(sepoch, 1, ntrials);
velz = nan(sepoch, ntrials);

for trId = 1:ntrials
    cstart = POS(trId) + EpochSamples(1);
    cstop  = POS(trId) + EpochSamples(2) - 1;
    T(:, :, trId) = eeg(cstart:cstop, :);
    E(:, :, trId) = eog(cstart:cstop, :);
    velz(:, trId) = twist.z(cstart:cstop, :);
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
[~, refchannelidx] = ismember(upper(refchannel), (settings.data.lchannels));
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
ncols = 12;

eogtimeslot     = 4;
eegtimeslot     = 7;
eegerrimgslot   = [2 5];
eegcorrimgslot  = [8 11];
velxerrimgslot  = [3 6];
velxcorrimgslot = [9 12];

get_slot_layout = @(slot) sort(reshape(((slot-1).*4 + [1 2 3 4]'), 1, length(slot).*4));

% Topoplot error
for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, tId);
    cstart = find(t >= topowins(tId, 1), 1, 'first');
    cstop  = find(t <= topowins(tId, 2), 1, 'last');
    topoplot(mean(m_eeg_err(cstart:cstop, :), 1), chanlocs32, 'maplimits', maplimits);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))])
end

% EOG signal
subplot(nrows, ncols, get_slot_layout(4));
hold on;
plot(t, m_eog_cor(:, 1), 'b');
plot(t, m_eog_err(:, 1), 'r');
hold off;
grid on;
plot_vline(0, 'k');
plot_hline(0, 'k');
xlim([t(1) t(end)]);
title(['subject: ' subject ' | channel: EOG | error vs. correct']);

% EEG signal
subplot(nrows, ncols, get_slot_layout(7));
plot_errp(t, m_eeg_cor(:, refchannelidx), m_eeg_err(:, refchannelidx));
plot_vline(0, 'k');
plot_hline(0, 'k');
title(['Subject: ' subject ' | channel: ' char(refchannel) ' | error vs. correct']);

% Topoplot correct
for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, tId + 36);
    cstart = find(t >= topowins(tId, 1), 1, 'first');
    cstop  = find(t <= topowins(tId, 2), 1, 'last');
    topoplot(mean(m_eeg_cor(cstart:cstop, :), 1), chanlocs32, 'maplimits', maplimits);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))])
end

% Imagesc error trials
subplot(nrows, ncols, get_slot_layout([2 5]));
imagesc(t, 1:sum(Ek), squeeze(T(:, refchannelidx, Ek))', [-15 15]);
plot_vline(0, 'k');
set(gca, 'YDir', 'normal');
title(['Subject: ' subject ' | channel: ' char(refchannel) ' | error trials']);
xlabel('time [s]')
ylabel('# trial')

% Imagesc correct trials
subplot(nrows, ncols, get_slot_layout([8 11]));
imagesc(t, 1:sum(Ck), squeeze(T(:, refchannelidx, Ck))', [-15 15]);
plot_vline(0, 'k');
set(gca, 'YDir', 'normal');
title(['Subject: ' subject ' | channel: ' char(refchannel) ' | correct trials']);
xlabel('time [s]')
ylabel('# trial')

overvelz = find_over(abs(velz), 0.2);

% imagesc on velz
subplot(nrows, ncols, get_slot_layout([3 6]));
hold on;
imagesc(t, 1:sum(Ek), squeeze(abs(velz(:, Ek)))');
axis image
axis normal
h = gca;
plot(t(overvelz(Ek)), 1:sum(Ek), '.w');
hold off;
xlim(h.XLim);
ylim(h.YLim);
plot_vline(0, 'w');
title(['Subject: ' subject ' | abs angular velocity | error trials']);
xlabel('time [s]')
ylabel('# trial')

subplot(nrows, ncols, get_slot_layout([9 12]));
hold on;
imagesc(t, 1:sum(Ck), squeeze(abs(velz(:, Ck)))');
axis image
axis normal
h = gca;
plot(t(overvelz(Ck)), 1:sum(Ck), '.w');
hold off;
xlim(h.XLim);
ylim(h.YLim);
plot_vline(0, 'w');
title(['Subject: ' subject ' | abs angular velocity | correct trials']);
xlabel('time [s]')
ylabel('# trial')

function x= find_over(v, value)

    x = ones(size(v, 2), 1);

    for trId = 1:size(v, 2)
        tmp = find(v(:, trId) >= value, 1, 'first');
        if(isempty(tmp) == false)
            x(trId) = tmp;
        end
    end

end


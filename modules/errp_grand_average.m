clearvars; clc;

subject = 'd6';
includepat  = {subject, 'discrete'};
excludepat  = {};
spatialfilter = 'car';
artifactrej   = 'none'; % {'FORCe', 'none'}

with_delay    = true;

eegchannels   = {'FP1', 'FP2', 'F1', 'FZ', 'F2', 'FC1', 'FCz', 'FC2', 'C1', 'CZ', 'C2', 'CP1', 'CP2'};

datapath      = ['analysis/' artifactrej '/' spatialfilter '/bandpass/'];
bagspath      = 'analysis/bags/aligned/';

datafiles = util_getfile3(datapath, '.mat', 'include', includepat, 'exclude', excludepat);
bagsfiles = util_getfile3(bagspath, '.mat', 'include', includepat, 'exclude', excludepat);

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

epoch = [-1 2];

%% Importing data
util_bdisp(['[io] - Importing ' num2str(ndatafiles) ' files from ' datapath ':']);
[eeg, eog, events, labels, settings] = errp_concatenate_bandpass(datafiles);
[pose, twist, cmdvel, bagevents, baglabels] = errp_concatenate_bags(bagsfiles);

nchannels  = size(eeg, 2);
samplerate = settings.samplerate;
eTYP       = events.TYP;
ePOS       = events.POS;
bTYP       = bagevents.TYP;
bPOS       = bagevents.POS;

%% Analyze events
util_bdisp('[proc] + Analysizing events:');
% Removing events within refractory period to avoid overlap while computing AMICA
refractory = length(epoch(1):1/samplerate:epoch(2));
disp(['       |- Find events within the refractory period (' num2str(refractory/samplerate, '%3.2f') ' s)']);
rmevt_refractory = errp_util_events_refractory(ePOS, refractory);

% Find the first event and last event consistent with the epoch
disp(['       |- Find events not consistent with epoch [' strjoin(compose('%g', floor(epoch*samplerate)), ' ') ']']);
firstevt    = find(ePOS > floor(abs(epoch(1)*samplerate)), 1, 'first');
lastevt     = find(ePOS < size(eeg, 1) - floor(abs(epoch(2)*samplerate)), 1, 'last');
rmevt_epoch = setdiff(1:length(eTYP), firstevt:lastevt);

rm_evtidx = union(rmevt_refractory, rmevt_epoch);
disp(['       |- Excluded events index: [' strjoin(compose('%g', rm_evtidx), ' ') ']']);
ePOS(rm_evtidx) = [];
eTYP(rm_evtidx) = [];


%% Extract trials
util_bdisp('[proc] + Extracting trials:');

% Extract events
disp('     |- Extract events');
evtidx = errp_util_get_event_type([CommandLx CommandRx ErrorLx ErrorRx NoReleaseFx NoReleaseRx NoReleaseLx], eTYP);
eTYP = eTYP(evtidx);
ePOS = ePOS(evtidx);

% Compute delay based on wheelchair velocity
refvel = 0.2184;
rmevt_velocity = [];
disp(['     |- Compute delays based on reference velocity: ' num2str(refvel)]);
ntrials = length(ePOS);
delays  = zeros(ntrials, 1);
for trId = 1:ntrials
    cstart = ePOS(trId);
    cstop  = ePOS(trId) + floor(epoch(2)*samplerate);
    overidx = find(abs(twist.z(cstart:cstop)) >= refvel, 1, 'first');
    if isempty(overidx) == true
        rmevt_velocity = cat(1, rmevt_velocity, trId);
    else
        delays(trId) = overidx;
    end
end

if with_delay == true
    ePOS = ePOS + delays;
end

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
disp(['     |- Extract epochs [' strjoin(compose('%g', epoch), ' ') '] seconds']);
nsamples  = length(epoch(1):1/samplerate:epoch(2));
ntrials   = length(ePOS);
T    = nan(nsamples, nchannels, ntrials);
E    = nan(nsamples, 1, ntrials);
velz = nan(nsamples, ntrials);

for trId = 1:ntrials
    cstart = ePOS(trId) + floor(epoch(1)*samplerate);
    cstop  = ePOS(trId) + floor(epoch(2)*samplerate);
    T(:, :, trId) = eeg(cstart:cstop, :);
    E(:, :, trId) = eog(cstart:cstop, :);
    velz(:, trId) = twist.z(cstart:cstop, :);
end

% Create label vectors
disp('     |- Create label vectors');

% Command correct
Ck = false(ntrials, 1);
correctindex = eTYP == CommandLx | eTYP == CommandRx;
Ck(correctindex) = true;

% Command error
Ek = false(ntrials, 1);
errorindex = eTYP == ErrorLx | eTYP == ErrorRx;
Ek(errorindex) = true;

%% Plotting
t = epoch(1):1/samplerate:epoch(2);

refchannel = {'FCz'};
[~, refchannelidx] = ismember(upper(refchannel), upper(settings.channels.eeg));
nrefchannel = length(refchannel);

chanlocs = settings.channels.chanlocs;

perc_err = 100*sum(Ek)./(sum(Ek) + sum(Ck));

m_eeg_err = squeeze(mean(T(:, :, Ek), 3));
m_eeg_cor = squeeze(mean(T(:, :, Ck), 3));
m_eog_err = squeeze(mean(E(:, :, Ek), 3));
m_eog_cor = squeeze(mean(E(:, :, Ck), 3));

topowins = [-0.5 -0.25; 0.25 0.5; 0.75 1.0; 1.25 1.5];

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

% EOG signal
subplot(nrows, ncols, get_slot_layout(1));
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
subplot(nrows, ncols, get_slot_layout(4));
plot_errp(t, m_eeg_cor(:, refchannelidx), m_eeg_err(:, refchannelidx));
plot_vline(0, 'k');
plot_hline(0, 'k');
title(['Subject: ' subject ' | channel: ' char(refchannel) ' | error vs. correct']);

% Topoplot error
htop = [];
for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, 24 + tId);
    cstart = find(t >= topowins(tId, 1), 1, 'first');
    cstop  = find(t <= topowins(tId, 2), 1, 'last');
    h = topoplot(mean(m_eeg_err(cstart:cstop, :), 1), chanlocs);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))]);
    htop = cat(1, htop, h);
end

% Topoplot correct
for tId = 1:size(topowins, 1)
    subplot(nrows, ncols, 36 + tId);
    cstart = find(t >= topowins(tId, 1), 1, 'first');
    cstop  = find(t <= topowins(tId, 2), 1, 'last');
    h = topoplot(mean(m_eeg_cor(cstart:cstop, :), 1), chanlocs);
    title([num2str(topowins(tId, 1)) '-' num2str(topowins(tId, 2))])
    htop = cat(1, htop, h);
end
errp_set_clim(htop);

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

figtitle = [subject ' grand average']; 
if with_delay == true
    figtitle = [figtitle ' with delay | threshold = ' num2str(refvel) ' m/s'];
end
sgtitle(figtitle);

function x= find_over(v, value)

    x = ones(size(v, 2), 1);

    for trId = 1:size(v, 2)
        tmp = find(v(:, trId) >= value, 1, 'first');
        if(isempty(tmp) == false)
            x(trId) = tmp;
        end
    end

end


clearvars; clc;

includepat  = {'delay'};
excludepat  = {};
depthlevel  = 2;

rootpath    = '/mnt/data/Research';
folder      = '2023_wheelchair_bag_delay';
bagpath     = [rootpath '/' folder ''];
savedir     = 'analysis/bags/raw/';

%% Processing parameters
samplerate  = 512;

Stop        = 100;
CommandLx   = 101;
CommandFx   = 102;
CommandRx   = 103;
ErrorLx     = 5101;
ErrorRx     = 5103;
NoReleaseLx = 4101;
NoReleaseFx = 4102;
NoReleaseRx = 4103;

JoyStop = 3;
JoyLx   = 5;    
JoyRx   = 6;
JoyFx   = 2;

%% Get datafiles
files = util_getfile3(bagpath, '.bag', 'include', includepat, 'exclude', excludepat, 'level', depthlevel);

nfiles = length(files);
if(nfiles > 0)
    util_bdisp(['[io] - Found ' num2str(nfiles) ' files with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
else
    error(['[io] - No files found with the inclusion/exclusion criteria: (' strjoin(includepat, ', ') ') / (' strjoin(excludepat, ', ') '), depth: ' num2str(depthlevel)]);
end

%% Create/Check for savepath
util_mkdir(pwd, savedir);

%% Processing files
delays = [];
for fId = 1:nfiles

    cfullname = files{fId};
    [cfilepath, cfilename, cfileext] = fileparts(cfullname);
    
    util_bdisp(['[io] + Loading file ' num2str(fId) '/' num2str(nfiles)]);
    disp(['     |-File: ' cfullname]);
    
    % Read the bag
    bag = rosbagreader(cfullname);
    
    % Select the available topics
    disp('     |-Selecting topics from bag');
    seldel    = select(bag, 'Topic', '/delay');

    % Read message for delays
    disp('     |-Extracting button and neuro events');
    msgdel = readMessages(seldel, 'DataFormat', 'struct');

    tmpdel = [];
    for iMsg = 1:length(msgdel)
        tmpdel = [tmpdel; msgdel{iMsg,1}.Data];
    end

    % Save the messages
    delays = [delays; tmpdel];

end

% Clear the to large delays
delays = delays(delays<=0.85);
delays = delays(delays>=0.1);

util_bdisp(['[io] + Compute std and means ']);

delay_mean = median(delays);
delay_std  = mad(delays);

fprintf('     |-Mean: %d \n',delay_mean);
fprintf('     |-Std:  %d \n',delay_std);

figure()
plot(delays)
hold on
plot(ones(length(delays))*delay_mean)

save('mean_delay', "delay_mean", "delay_std")

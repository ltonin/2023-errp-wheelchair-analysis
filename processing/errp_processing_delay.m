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

time_no_release = round(512 * 3);

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
vels = [];
full_events.POS = [];
full_events.ERR = [];
full_events.TYP = [];

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
    selodom   = select(bag, 'Topic', '/odometry/filtered');
    seljoy    = select(bag, 'Topic', '/j0');

    % Read message for delays
    disp('     |-Extracting button and neuro events');
    msgdel   = readMessages(seldel, 'DataFormat', 'struct');
    msgjoy = readMessages(seljoy, 'DataFormat', 'struct');


    nuodom   = timeseries(selodom, 'Twist.Twist.Angular.Z');
    
    % Extract button press
    is_button_valid = cellfun(@(m) isempty(find(m.Buttons, 1)) == false, msgjoy);
    button_value    = cellfun(@(m) find(m.Buttons, 1), msgjoy, 'UniformOutput',false);
    button_value    = cell2mat(button_value(is_button_valid));
    button_time     = cellfun(@(m) double(m.Header.Stamp.Sec) + double(m.Header.Stamp.Nsec)*10^-9, msgjoy, 'UniformOutput', false);
    button_time     = cell2mat(button_time(is_button_valid));

    % Resample timeseries
    disp('     |-Resampling timeseries');
    t_common = nuodom.Time(1):1/samplerate:nuodom.Time(end)-1/samplerate;
    rodom   = resample(nuodom, t_common);

    % Align button and neuro event
    disp('     |-Aligning events to resampled time support');
    t_align = rodom.Time(1):1/samplerate:rodom.Time(end)-1/samplerate;
    [buttonPOS, buttonERR] = errp_find_closest(button_time, t_align);

    twist_z = rodom.Data(:, 1);

    % Create event structure
    buttonTYP = button_value;
    [events.POS, sindex] = sort(buttonPOS);
    events.ERR = buttonERR(sindex);
    events.TYP = buttonTYP(sindex);
    events.POS = events.POS' + length(vels);
    events.TYP = events.TYP;
    events.ERR = events.ERR';

    tmpdel = [];
    for iMsg = 1:length(msgdel)
        tmpdel = [tmpdel; msgdel{iMsg,1}.Data];
    end

    % Save the messages
    delays = [delays; tmpdel];

    % Save the messages
    vels = [vels; twist_z];
    % Save the joy
    full_events.POS = [full_events.POS; events.POS];
    full_events.TYP = [full_events.TYP; events.TYP];
    full_events.ERR = [full_events.ERR; events.ERR];


end

% Clear the to large delays
delays = delays(delays<=0.85);
delays = delays(delays>=0.1);

util_bdisp(['[io] + Compute std and means ']);

delay_mean = median(delays); % <-- IT IS A MEDIAN!!!
delay_std  = mad(delays);

fprintf('     |-Median:   %d \n',delay_mean);
fprintf('     |-Mad:      %d \n',delay_std);

figure()
plot(delays)
title("lenght of the delays")
hold on
plot(ones(length(delays))*delay_mean)

% Find the joy events that I need to compute the delay
evsidx = errp_util_get_event_type([ JoyLx JoyRx ], full_events.TYP);
evsidx = errp_reduce_events(full_events, evsidx, time_no_release);

delays = ones(length(full_events.POS),1);

% Move to the odom samplerate
frame = 1/samplerate;
pulse = delay_mean/frame;
delay_mean_pulsation = round(pulse);

delays = delays*delay_mean_pulsation; 

delayed_events.POS = full_events.POS + delays;
% Compute the mean velocity when the person press the button
abs_vels = abs(vels);
delay_vel = mean( abs_vels(delayed_events.POS(evsidx)) );
fprintf('     |-Vel Mean: %d \n', delay_vel);

% Save the data
save('mean_delay', "delay_mean", "delay_std", "delay_vel")

% Plot the data
plot_delay(vels,  full_events.POS(evsidx), delays(evsidx), false)






















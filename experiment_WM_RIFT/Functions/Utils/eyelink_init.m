function prm = eyelink_init(prm,id)

%% File and directory for the .edf file
prm.exp.eyeFile = sprintf('s%02d.edf',id.SubNo);
prm.exp.eyeDir = fullfile(id.path,'edf'); 
mkdir(prm.exp.eyeDir); % this folder to be created if it does not exist!

%% Initialize
ELsamplingrate = 1000;
dummymode = ~prm.exp.useEye;
EyelinkInit(dummymode); % 0 -> don't ask for dummymode, 1 -> dummymode

prm.eyelink = EyelinkInitDefaults(prm.w.Number); % sends pixel coordinates to the eyetracker

%% Colour choices for calibration
% We are changing calibration to a black background with almost white
% targets, no sound and smaller targets
prm.eyelink.backgroundcolour = 128;    % BlackIndex(prm.eyelink.expWin);
prm.eyelink.foregroundcolour = 250;  % not completely white % white;
prm.eyelink.msgfontcolour = 250;     %WhiteIndex(prm.eyelink.expWin);
prm.eyelink.imgtitlecolour = 250;    %WhiteIndex(prm.eyelink.expWin);
% set calibration sounds - 0 off, 1 on
prm.eyelink.targetbeep = 0;
prm.eyelink.feedbackbeep = 0;
prm.eyelink.calibration_target_beep = 0;
prm.eyelink.calibrationtargetcolour = prm.eyelink.foregroundcolour; % has to be same, because of bug in 'EyelinkDrawCalTarget.m'
% 'width' (inner circle) has to be smaller than 'size' (outer circle),
% otherwise calibration target has background colour, i.e. it is
% invisible, see >> open EyelinkDrawCalTarget
% prm.eyelink.calibrationtargetsize = .5;
% prm.eyelink.calibrationtargetwidth = .3; % property of "size"
prm.eyelink.calibrationtargetsize = .8;
prm.eyelink.calibrationtargetwidth = .3; % property of "size"

EyelinkUpdateDefaults(prm.eyelink); % update the defauls with new specifications

Eyelink('Command', 'set_idle_mode'); % Put tracker in idle/offline mode before recording
WaitSecs(0.05); % Allow some time for transition    

%% Set-up datafile -----------------------------------------------------
% Open an EDF file
failOpen = Eyelink('OpenFile', prm.exp.eyeFile); % 0 success, else error code
if failOpen ~= 0 % Abort if it fails to open
    error('Cannot create EDF file %s', prm.exp.eyeFile); % Print some text in Matlab's Command Window
end
preambleText = sprintf('Recorded by EyelinkToolbox, Experiment narratives, %s, subject-%02d''', prm.exp.dateString, id.SubNo);
Eyelink('command', 'add_file_preamble_text "%s"', preambleText);
WaitSecs(0.05);
Eyelink('message', '%s',prm.exp.dateString);
WaitSecs(0.05);

%% Set-up screen res ---------------------------------------------------
Eyelink('command','screen_distance = %ld %ld', prm.exp.edgeEyeDistTop, prm.exp.edgeEyeDistBottom);
Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, prm.monitor.width-1, prm.monitor.height-1);
% ... continued below!

%% Set calibration type ------------------------------------------------
Eyelink('command', 'calibration_type = HV9'); % use 9-point, since that's recommended for the long range setup as it is in the MEG lab (says Davide)
% manually specify location of calibration targets
Eyelink('command', 'generate_default_targets = NO');
% for 9-point calibration
% 5 1 6
% 3 0 4
% 7 2 8

%% here define the position of calibration point
% stimX = mean(prm.stim.DistanceInPxl);
% stimY = mean(prm.stim.DistanceInPxl);

stimX = 1200;%prm.stim.SizeInPxl;
stimY = 900;%prm.stim.SizeInPxl; 
% make the calibration targets be roughly the size of the image
Eyelink('command', 'calibration_targets = %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld', ...
    prm.monitor.width/2,prm.monitor.height/2, prm.monitor.width/2,prm.monitor.height/2-stimY/2, prm.monitor.width/2,prm.monitor.height/2+stimY/2,...
    prm.monitor.width/2-stimX/2,prm.monitor.height/2, prm.monitor.width/2+stimX/2,prm.monitor.height/2, prm.monitor.width/2-stimX/2,prm.monitor.height/2-stimY/2,...
    prm.monitor.width/2+stimX/2,prm.monitor.height/2-stimY/2, prm.monitor.width/2-stimX/2,prm.monitor.height/2+stimY/2, prm.monitor.width/2+stimX/2,prm.monitor.height/2+stimY/2);
Eyelink('command', 'validation_targets = %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld', ...
    prm.monitor.width/2,prm.monitor.height/2, prm.monitor.width/2,prm.monitor.height/2-stimY/2, prm.monitor.width/2,prm.monitor.height/2+stimY/2,...
    prm.monitor.width/2-stimX/2,prm.monitor.height/2, prm.monitor.width/2+stimX/2,prm.monitor.height/2, prm.monitor.width/2-stimX/2,prm.monitor.height/2-stimY/2,...
    prm.monitor.width/2+stimX/2,prm.monitor.height/2-stimY/2, prm.monitor.width/2-stimX/2,prm.monitor.height/2+stimY/2, prm.monitor.width/2+stimX/2,prm.monitor.height/2+stimY/2);
% Eyelink('command', 'calibration_targets = %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld', ...
%     prm.monitor.width/2,prm.monitor.height/2, prm.monitor.width/2,prm.monitor.height/4, prm.monitor.width/2,prm.monitor.height*3/4 ,...
%     prm.monitor.width/4,prm.monitor.height/2, prm.monitor.width*3/4,prm.monitor.height/2, prm.monitor.width/4,prm.monitor.height/4,...
%     prm.monitor.width*3/4,prm.monitor.height/4, prm.monitor.width/4,prm.monitor.height*3/4, prm.monitor.width*3/4,prm.monitor.height*3/4);
% Eyelink('command', 'validation_targets = %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld %ld,%ld  %ld,%ld  %ld,%ld  %ld,%ld', ...
%     prm.monitor.width/2,prm.monitor.height/2, prm.monitor.width/2,prm.monitor.height/4, prm.monitor.width/2,prm.monitor.height*3/4 ,...
%     prm.monitor.width/4,prm.monitor.height/2, prm.monitor.width*3/4,prm.monitor.height/2, prm.monitor.width/4,prm.monitor.height/4, ...
%     prm.monitor.width*3/4,prm.monitor.height/4, prm.monitor.width/4,prm.monitor.height*3/4, prm.monitor.width*3/4,prm.monitor.height*3/4);

% has to be executed after calibration targets were specified:
Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, prm.monitor.width-1, prm.monitor.height-1);

% only do monocular recording
Eyelink('command', 'binocular_enabled = YES'); WaitSecs(0.05); % see Eyelink Programmers Guide
% we want pupil area, not diameter
Eyelink('command', 'pupil_size_diameter = YES'); WaitSecs(0.05);

%% Set sampling rate ---------------------------------------------------
% set above.
Eyelink('command', 'sample_rate = %d', ELsamplingrate); WaitSecs(0.05);

%% Set movement thresholds (conservative) ------------------------------ 
Eyelink('command', 'saccade_velocity_threshold = 35'); WaitSecs(0.05);
Eyelink('command', 'saccade_acceleration_threshold = 9500'); WaitSecs(0.05);

%% Get tracker and software versions -----------------------------------
[~,vs] = Eyelink('GetTrackerVersion');
fprintf('Running experiment on a ''%s'' tracker.\n', vs);

%% Saving and linking data
% Select which files are saved to edf data ----------------------------------------------------
Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
WaitSecs(0.05);
Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
WaitSecs(0.05);

% Link data to Matlab -------------------------------------------------
Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
WaitSecs(0.05);
Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT,HTARGET');
WaitSecs(0.05);
Eyelink('command', 'link_event_data = GAZE,GAZERES,HREF,AREA,VELOCITY');
WaitSecs(0.05);

%% Set idle mode
Eyelink('command', 'set_idle_mode'); 
WaitSecs(0.05);
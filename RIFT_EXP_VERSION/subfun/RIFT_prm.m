function prm    = RIFT_prm(session)
% parameters for TW_SSVEP
%%=========================================================================
% session
default('session','');
%%=========================================================================
% filepath  
prm.main = 'D:\Users\sonbai\rift_pre_pilot-master\';
prm.path               = fileparts(mfilename('fullpath')); 
% cd(prm.path);cd('..\');
%%=========================================================================
% general exp prm
prm.exp.is_live        = 1;
prm.exp.eyelink_live   = 1;
prm.exp.view_distance  = 100;
prm.exp.edgeEyeDistTop = 650;
prm.exp.edgeEyeDistBottom = 665;
prm.exp.name           = 'RIFT';
prm.exp.dateString     = datestr(clock(), 'yyyy-mm-dd_HH-MM-SS');
prm.exp.path           = pwd;
fileID                 = fopen('instructions.txt','r');
prm.exp.instructions   = fscanf(fileID,'%c');fclose(fileID);
prm.exp.blocks         = 10;
prm.exp.ntrials        = 320;
prm.exp.ntrialsperblock = prm.exp.ntrials / prm.exp.blocks;
prm.exp.ndatap         = 10; % number of repetition in the whole experiment
% prm.exp.datetime       = datetime;
% for the practice session
if strcmp(session,'practice')
    prm.exp.is_live      = 1;
    prm.exp.eyelink_live = 1;
    prm.exp.name         = 'RIFT_practice';
    prm.exp.blocks       = 4;
    prm.exp.ntrials      = 128;
    prm.exp.ndatap       = 4;
end
%%========================================================================
% monitor prm

% prm.monitor                = InitPsychtoolbox(prm.exp.is_live);
prm.monitor            = monitor_info([],prm.exp.view_distance);
% add gamma correction table
% load(['/Users/gzy/Desktop/Experiments/Functions/monitor_cal/'   ...
%       'SmallLeft_monitor_calib_22_10_2019.mat'],'gammaTable1');
% prm.monitor.gamma      = gammaTable1;
    
%%=========================================================================
% stimuli prm
prm.stim.path          = fullfile(prm.main, '\stimuli_set\');
prm.stim.AllStim       = dir(append(prm.stim.path, '\**\*.png'));
prm.stim.AllStimPath   = [];
for i = 1:numel(prm.stim.AllStim)
    stim_path = append(prm.stim.AllStim(i).folder, '\', prm.stim.AllStim(i).name);
    prm.stim.AllStimPath   = [prm.stim.AllStimPath convertCharsToStrings(stim_path)];
end

% freqeuncy setting
prm.stim.refresh       = prm.w.RefreshRate * 12;
prm.stim.dt            = 1 / prm.stim.refresh; % in second
prm.stim.tag_frex      = [55,60,65,70];%[2 4 8 1];
prm.stim.duration      = 4; % in second
prm.stim.timax         = prm.stim.dt:prm.stim.dt:prm.stim.duration;
prm.stim.tag_sigs      = cos(2*pi*prm.stim.tag_frex(:)*prm.stim.timax) / 2 + 0.5;


% spatial paramter========================
prm.stim.SizeInDeg     = 3.5;
prm.stim.SizeInPxl     = round(prm.monitor.deg_to_pix * prm.stim.SizeInDeg); % vis deg
prm.stim.size          = [prm.stim.SizeInPxl prm.stim.SizeInPxl];
prm.stim.grating       = repmat(zeros(prm.stim.size),[1 1 4]);
for quardrant = 1:4 % center of each sub-screen
    prm.stim.subcenters(quardrant,:)    = MakeOffsetRect(prm.monitor, prm.stim, 0, 0,quardrant);
end
% define the position of stimuli in polar coordinate                      
prm.stim.DistanceInDeg = [4 6];
prm.stim.DistanceInPxl = round(prm.monitor.deg_to_pix * prm.stim.DistanceInDeg);% define the position of stimuli in polar coordinate
prm.stim.AngleRange    = [-15 15];
prm.stim.AngleInDeg    = repmat([225; 315; 135; 45],[1, 2]) + prm.stim.AngleRange; %

% prm.stim.distance      = [-prm.monitor.width/16, -prm.monitor.height/16; prm.monitor.width/16, -prm.monitor.height/16;...
%                           -prm.monitor.width/16, prm.monitor.height/16; prm.monitor.width/16, prm.monitor.height/16];
prm.stim.JitterInDeg   = [0 0];  % with jitter
prm.stim.JitterInInPxl = round(prm.monitor.deg_to_pix * prm.stim.JitterInDeg);

%%=========================================================================
% MEG triggers
if prm.exp.is_live
  prm.trigger.btsi = Bitsi('COM1'); % or whichever COM-port used by the PC
else

  prm.trigger.btsi.validResponses = KbName({'a' 'e'});
end
prm.trigger.ExpStart    = 9;
prm.trigger.FixStart    = 8; % start of a trial
prm.trigger.RiftStart   = 1; % stimuli onset
prm.trigger.FixEnd      = 88; % this is only necessary for fixation after sti onset
prm.trigger.RiftEnd     = 11; % stimuli offset
prm.trigger.MemorySti   = 4;
prm.trigger.response    = 6;
prm.trigger.ExpEnd      = 99;

%%=========================================================================
% fixation

%% Visual content==========================================================

% Fixation size

prm.fix.sizeAngle = 0.15;

prm.fix.size = visAng2pix(prm.fix.sizeAngle, 0, prm.monitor);
prm.fix.col =  [0 0 0];
prm.fix.col2 =  [200 0 0];



prm.fixation.color      = [0 0 0]; % fixation dot color
prm.fixation.size       = fix(.5*prm.monitor.deg_to_pix).*ones(1,2); 
% prm.fixation.pos        = position_shifter(prm.monitor,prm.fixation.size,[0 0]); 
% prm.fixation.PreDur     = prm.monitor.frame_ms(1500) * 12; % in physical frames
% prm.fixation.Dur        = prm.monitor.frame_ms(2000) * 12; % duration of fixation after stimuli onset in frames
prm.fixation.PreFixRange= 1000:1:1300; % range of fixation duration in ms
prm.fixation.FixRange   = 1000:1:1300; % duration of fixation after stimuli onset in frames
prm.fixation.rad_outer  = 0.4;
prm.fixation.rad_inner  = 0.2;

prm.fixation.maxDevDeg  = 2;% max deviation of fixation in visual degree 
prm.fixation.maxDevPxl  = round(prm.monitor.deg_to_pix * prm.fixation.maxDevDeg);
%%=========================================================================
% background
prm.bg.path          = fullfile(prm.main, '\background.jpg');
bg                   = imresize(imread(prm.bg.path), [prm.monitor.width/2, prm.monitor.height/2]); 
prm.bg.grating       = bg;
% initialize drawing rects (centered on physical screen)
prm.bg.positions     = zeros(4, 4);
for quadrant = 1:4
    prm.bg.positions(quadrant,:) = MakeOffsetRect(prm.monitor, prm.bg, 0, 0, quadrant);
end


%%=========================================================================
% diode tracking square
prm.diode_track.SizeInDeg = 1; % vis deg
prm.diode_track.SizeInPxl = round(prm.monitor.deg_to_pix * prm.diode_track.SizeInDeg);
prm.diode_track.freq      = 2; % which tagging freq ro check[1 2 3 4]
prm.diode_track.grating   = ones(prm.diode_track.SizeInPxl,prm.diode_track.SizeInPxl);
for quadrant = 1:4
    prm.diode_track.positions(quadrant,:) = MakeOffsetRect(prm.monitor, prm.diode_track, 0, 0, quadrant);
end

%%=========================================================================
% Memory Task
prm.task.nShown          = 2; 
prm.task.nUnShown        = 2;  

%%=========================================================================
% generate the exp design and conditions
CurvOrRect               = [1,2];  % whether distractors has a rect or curv outline: 1 for Curv, 2 for Rect
LargeOrSmall             = [1,2];  % whether distractors is large or small: 1 for large, 2 for small
OddballShapeOrSize       = [1,2];  % whether there is a shape oddball or Size oddball; 1 for shape oddball, 2 for size oddball
OddballPos               = [1 2 3 4]; % in which quardrant the oddball is placed; 1 upper left, 2 upper right, 3 lower left, 4 lower right
prm.exp.design           = create_design({OddballShapeOrSize,CurvOrRect,LargeOrSmall,OddballPos},prm.exp.ndatap,prm.exp.blocks);
prm.exp.design           = make_table({prm.exp.design,'OddballShapeOrSize','CurvOrRect','LargeOrSmall','OddballPos','block','trial'});
% add left/right attention label
prm.exp.design.OddCurvOrRect   = 2 - mod(prm.exp.design.CurvOrRect + prm.exp.design.OddballShapeOrSize,2);
prm.exp.design.OddLargeOrSmall = 2 - (prm.exp.design.LargeOrSmall == prm.exp.design.OddballShapeOrSize);
prm.exp.design.fixmode         = repmat([zeros(prm.exp.ntrialsperblock, 1); ones(prm.exp.ntrialsperblock, 1)], prm.exp.blocks/2, 1);

if RealRun
    prm.time.before = 12*1000;
else
    prm.time.before = 3*1000;
end

prm.time.countdown = 3;


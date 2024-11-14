function prm    = prm_RIFT(RunType, RealRun, Environment, useEye)

% parameters for TW_SSVEP


%%=========================================================================
% general exp prm
prm.exp.is_live        = RealRun;
prm.exp.eyelink_live   = useEye;
prm.exp.RunType        = RunType;
prm.exp.edgeEyeDistTop = 650;
prm.exp.edgeEyeDistBottom = 665;
prm.exp.dateString     = datestr(clock(), 'yyyy-mm-dd_HH-MM-SS');
prm.exp.path           = pwd;
prm.exp.datedir =      'Data_Lab';
% for the practice session
if strcmp(RunType,'pra')
    prm.N.praTrial       = 5;
    prm.exp.datedir =      'Data_Demo';

end
%%========================================================================
% monitor prm

% prm.monitor                = InitPsychtoolbox(prm.exp.is_live);
prm.monitor            = monitor_info([], Environment);
% add gamma correction table
% load(['/Users/gzy/Desktop/Experiments/Functions/monitor_cal/'   ...
%       'SmallLeft_monitor_calib_22_10_2019.mat'],'gammaTable1');
% prm.monitor.gamma      = gammaTable1;
    
%% temporal proterties ====================================================
prm.time.before    = 5;  % a resting period 
prm.time.countdown = 3;

% durations

prm.dur.fix    = 1200;
prm.dur.image1 = 1200;
prm.dur.mask   = 1000;
prm.dur.delay  = 6400;
prm.dur.wait   = 7000; % for response

prm.Nfr        = structfun(@(x) round(x/1000 * prm.monitor.hz), prm.dur,'UniformOutput',0); % need to
prm.SwitchFr   = cumsum(cell2mat(struct2cell(prm.Nfr)));
prm.dur.trial  = sum(cell2mat(struct2cell(prm.dur)));


%% tagging strips

if Environment==1
    prm.tag.refresh       = prm.monitor.hz ;
    prm.tag.tag_frex      = 30*ones(21,1);
elseif Environment==2 % MEG lab
    prm.tag.refresh       = prm.monitor.hz * 12;
    prm.tag.tag_frex      = 55:0.5:65;
end
prm.tag.dt            = 1 / prm.tag.refresh; % in second

prm.tag.duration      = prm.dur.delay; % in second
prm.tag.timax         = prm.tag.dt:prm.tag.dt:prm.tag.duration;
prm.tag.tag_sigs      = cos(2*pi*prm.tag.tag_frex(:)*prm.tag.timax) / 2 + 0.5;

prm.tag.angle         = 20:2.5:70;
prm.tag.nStrip        = length(prm.tag.angle);
prm.tag.stripWidth    = 2; % in angle (degree)

% for quardrant = 1:4 % center of each sub-screen
%     prm.tag.subcenters(quardrant,:)    = MakeOffsetRect(prm.monitor, prm.tag, 0, 0,quardrant);
% end
% % define the position of stimuli in polar coordinate                      
% prm.tag.DistanceInDeg = [4 6];
% prm.tag.DistanceInPxl = round(prm.monitor.deg_to_pix * prm.tag.DistanceInDeg);% define the position of stimuli in polar coordinate
% prm.tag.AngleRange    = [-15 15];
% prm.tag.AngleInDeg    = repmat([225; 315; 135; 45],[1, 2]) + prm.tag.AngleRange; %
% 
% % prm.tag.distance      = [-prm.monitor.width/16, -prm.monitor.height/16; prm.monitor.width/16, -prm.monitor.height/16;...
% %                           -prm.monitor.width/16, prm.monitor.height/16; prm.monitor.width/16, prm.monitor.height/16];
% prm.tag.JitterInDeg   = [0 0];  % with jitter
% prm.tag.JitterInInPxl = round(prm.monitor.deg_to_pix * prm.tag.JitterInDeg);

%%=========================================================================
% MEG triggers
if RealRun==1 && Environment==2
  prm.trigger.btsi = Bitsi('COM1'); % or whichever COM-port used by the PC
  prm.trigger.btsi.validResponses = ['a' 'e']; % index finger button on right (a)/left (e) button boxes
% else
%   prm.trigger.btsi = Bitsi('COM1');
%   prm.trigger.btsi.validResponses = KbName({'a' 'e'});
end

prm.trigger.ExpStart      = 9;
prm.trigger.FixStart      = 0; % start of a trial
prm.trigger.TargetStart   = 11; % stimuli onset
prm.trigger.TargetEnd     = 12; 
prm.trigger.DelayStart    = 21; 
prm.trigger.DelayEnd      = 22;
prm.trigger.Response      = 7;
prm.trigger.ExpEnd        = 99;

%%=========================================================================
% fixation
prm.fix.color      = [0 0 0]; % fixation dot color
prm.fix.color2     = [200, 0, 0];
prm.fix.size       = visAng2pix(0.15, 0, prm.monitor); 
% prm.fix.pos        = position_shifter(prm.monitor,prm.fix.size,[0 0]); 
% prm.fix.PreDur     = prm.monitor.frame_ms(1500) * 12; % in physical frames
% prm.fix.Dur        = prm.monitor.frame_ms(2000) * 12; % duration of fixation after stimuli onset in frames
prm.fix.PreFixRange= 1000:1:1300; % range of fixation duration in ms
prm.fix.FixRange   = 1000:1:1300; % duration of fixation after stimuli onset in frames
prm.fix.rad_outer  = 0.4;
prm.fix.rad_inner  = 0.2;

prm.fix.maxDevDeg  = 2;% max deviation of fixation in visual degree 
prm.fix.maxDevPxl  = round(prm.monitor.deg_to_pix * prm.fix.maxDevDeg);


%%=========================================================================
% diode tracking square
prm.diode_track.SizeInDeg = 1; % vis deg
prm.diode_track.SizeInPxl = round(prm.monitor.deg_to_pix * prm.diode_track.SizeInDeg);
prm.diode_track.freq      = 2; % which tagging freq ro check[1 2 3 4]
prm.diode_track.grating   = ones(prm.diode_track.SizeInPxl,prm.diode_track.SizeInPxl);
for quadrant = 1:4
    prm.diode_track.positions(quadrant,:) = MakeOffsetRect(prm.monitor, prm.diode_track, 0, 0, quadrant);
end



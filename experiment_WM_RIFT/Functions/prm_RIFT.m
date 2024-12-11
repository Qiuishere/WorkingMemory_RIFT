function prm    = prm_RIFT(SubNo, RunType, RealRun, Environment, useEye)

% parameters for tagging, stimuli, and conditions


%%=========================================================================
% general exp prm
prm.exp.SubNo          = SubNo;
prm.exp.RealRun        = RealRun;
prm.exp.useEye   = useEye;
prm.exp.Environment    = Environment;
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

%% Make directories for participant (if needed)============================
prm.subInfo.SubNo = SubNo;
prm.subInfo.path = fullfile(prm.exp.datedir, sprintf('Sub%02d', SubNo), RunType);

if ~exist(prm.subInfo.path, 'dir')
    mkdir(prm.subInfo.path);
end


%% determine order
if mod(SubNo,2)
    prm.exp.body1st = 1;
else
    prm.exp.body1st = 0;
end

if mod(str2double(prm.exp.dateString(end-1:end)),2)
    prm.tag.order = 1; % from small to large
else
    prm.tag.order = 0; % from large to small
end

%%========================================================================
% monitor prm
prm.monitor            = monitor_info([], Environment);
prm.monitor.mode       = [0,2,5];
prm.monitor.mode_fr     = [120, 480, 1440];
% add gamma correction table
% load(['/Users/gzy/Desktop/Experiments/Functions/monitor_cal/'   ...
%       'SmallLeft_monitor_calib_22_10_2019.mat'],'gammaTable1');
% prm.monitor.gamma      = gammaTable1;
    
%% temporal proterties ====================================================
prm.time.before    = 5;  % a resting period, in second
prm.time.countdown = 3;

% durations
 
prm.dur.image1 = 1200; % in ms
prm.dur.mask   = 1000;
prm.dur.delay  = 2400; % ms
prm.dur.wait   = 7000; % for response

prm.Nfr        = structfun(@(x) round(x/1000 * prm.monitor.hz), prm.dur,'UniformOutput',0); % need to
prm.SwitchFr   = cumsum(cell2mat(struct2cell(prm.Nfr)));
prm.dur.trial  = sum(cell2mat(struct2cell(prm.dur)));

prm.dur.fixRange = prm.monitor.frame_ms(1000:prm.monitor.dt*1000:1300);

prm.exp.Nmask = 6;
prm.Nfr.maskId = 1: prm.Nfr.mask/prm.exp.Nmask : prm.Nfr.mask;
%% tagging strips

if Environment==1
    prm.tag.refresh       = prm.monitor.hz *12 ;
    prm.tag.tag_frex      = linspace(3, 8, 11);
elseif Environment==2 % MEG lab
    prm.tag.refresh       = prm.monitor.hz * 12;
    prm.tag.tag_frex      = 55:10:65;%:0.5:65;
end

if prm.tag.order == 0
    prm.tag.tag_fex       = fliplr(prm.tag.tag_frex);
end


prm.tag.dt            = 1 / prm.tag.refresh; % in second

prm.tag.duration      = (prm.dur.image1 + prm.dur.mask + prm.dur.delay)/1000*3 ; % in second. *2 to ensure enough time
prm.tag.timax         = prm.tag.dt:prm.tag.dt:prm.tag.duration;
prm.tag.tag_sigs      = cos(2*pi*prm.tag.tag_frex(:)*prm.tag.timax) / 2 + 0.5;

prm.tag.angle         = 20:5:70;
prm.tag.nStrip        = length(prm.tag.tag_frex);
prm.tag.stripWidth    = 4; % in angle (degree)

% Spatial parameters
prm.img.WAng   = 8;
prm.img.WPix = visAng2pix(prm.img.WAng, 0, prm.monitor);


%%=========================================================================
% diode tracking square
prm.diode_track.SizeInDeg = 1; % vis deg
prm.diode_track.SizeInPxl = round(prm.monitor.deg_to_pix * prm.diode_track.SizeInDeg);
prm.diode_track.freq      = 2; % which tagging freq to check in the order of prm.tag.tag_frex
prm.diode_track.grating   = ones(prm.diode_track.SizeInPxl,prm.diode_track.SizeInPxl);



%%=========================================================================
% MEG triggers
if RealRun==1 && Environment==2
  prm.trigger.btsi = Bitsi('COM1'); % or whichever COM-port used by the PC
  prm.exp.key.up.press     = 99;
  prm.exp.key.up.release   = 67;
  prm.exp.key.down.press   = 100;
  prm.exp.key.down.release = 68;
  prm.exp.key.space        = 104;
  prm.trigger.btsi.validResponses = [prm.exp.key.up.press, prm.exp.key.down.press, prm.exp.key.up.release, prm.exp.key.down.release, prm.exp.key.space]; % R index finger, R middle, and L index

end

prm.trigger.ExpStart           = 1;
prm.trigger.FixStart           = 11; % start of a trial
prm.trigger.TargetStart.bar    = 21; % stimuli onset
prm.trigger.TargetStart.body   = 22; % stimuli onset
prm.trigger.TargetEnd          = 29; 
prm.trigger.DelayStart         = 31; 
prm.trigger.DelayEnd           = 39;
prm.trigger.Response           = 40;
prm.trigger.ExpEnd             = 99;

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
% prm.diode_track.SizeInDeg = 1; % vis deg
% prm.diode_track.SizeInPxl = round(prm.monitor.deg_to_pix * prm.diode_track.SizeInDeg);
% prm.diode_track.freq      = 2; % which tagging freq ro check[1 2 3 4]
% prm.diode_track.grating   = ones(prm.diode_track.SizeInPxl,prm.diode_track.SizeInPxl);
% for quadrant = 1:4
%     prm.diode_track.positions(quadrant,:) = MakeOffsetRect(prm.monitor, prm.diode_track, 0, 0, quadrant);
% end


%% Constructing all the trials=========================================
    prm.fac.adjustRange = (25:65)';
    prm.fac.targetRange = (30:5:60)';
    prm.fac.views = {'L', 'R'}; % 1: facing left (presented onthe right)
    prm.fac.figures = {'Fe'; 'Ma'};
    prm.fac.targetIds = [1;2];
    
    % Generate all combinations using ndgrid
    [t, d, v, f] = ndgrid(prm.fac.targetIds, prm.fac.targetRange, 1:length(prm.fac.views),  prm.fac.figures);
    
    % Combine into a cell array
    combinations = [num2cell(t(:)), num2cell(d(:)), num2cell(v(:)),  f(:)];
    
    % Create the table
    table = cell2table(combinations, ...
        'VariableNames', {'targetId','angleTarget', 'view1',  'figure1', });
   
   % table = repmat(table, 3, 1);
    
  
    table.angleProbe = randSamp(prm.fac.adjustRange, height(table), 'n');
       table = table(randperm(height(table)),:);
       
     % randomize and then make a copy for bars, counterbalance the order
      
    prm.fac.stimuli = ["Body","Bar"];
    table.stimulus = repmat("Body", height(table), 1);
    
    if  strcmp(prm.exp.RunType , 'pra')
        table = table(1: prm.N.praTrial,:);
    end
    
    table2 = table;
    table2.stimulus = repmat("Bar", height(table), 1);
    table2.figure1  = repmat([" "],height(table), 1);
    
    prm.N.trial = height(table)*2;
    prm.N.block = 4*2;
    if  strcmp(prm.exp.RunType , 'pra')
        prm.N.block = 1*2;
    end    
    prm.N.trialPerBlock = prm.N.trial/prm.N.block;
    
    T = [];
    for thebl = 1: 0.5*prm.N.block
        thetrialid = (1:prm.N.trialPerBlock) + prm.N.trialPerBlock*(thebl-1);
        if prm.exp.body1st==1
            T = [T; table(thetrialid,:); table2(thetrialid,:)];
        else
            T = [T; table2(thetrialid,:); table(thetrialid,:)];
        end
    end
     T.block = reshape(repmat(1:prm.N.block, prm.N.trialPerBlock,1),[],1);
   
    T.resp = zeros(height(T),1);
    T.rt   = NaN(height(T),1);
    T.subjectid = SubNo * ones(height(T),1);
    T = movevars(T, 'subjectid','Before', 1);
    
        
    prm.exp.T = T;
    
    prm.exp.inst = ['In each trial, you will see either a thin bar or a person lifting their arm\n\n',...
        'Please remember the position of the bar or the arm when you see them,\n\n later you will need to reproduce that position.\n\n',...
        'Then there will be some checkerboard images, followed by a blank screen, then the bar or the person will appear again.\n\n',...
        'And the fixation will turn into red. At this moment you should adjust the arm or the bar to the previous position.\n\n',...
        'You can press upkey or downkey to adjust, then press Space to confirm. You have 7 seconds to respond.\n\n',...
        'If you don''t confirm in 7s, the next trial will start regardless.\n\n',...
        'One important thing: Please make sure to always stare at the fixation point at the center. \n\n'...
        'If you don''t have any questions, Press Space to start!'];


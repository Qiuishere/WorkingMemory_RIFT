%function Body_RunRS(SubNo, RunInfo)

%{
the single-target design with flickering mask

This is the most updated version. can be easily adapted into an retro-cue
paradigm 
%}
sca; clear all; close all; clc; commandwindow;


if ~exist( 'SubNo', 'var') % this will run if the "function" line is commented
    
    addpath(genpath('Functions'));
    StartUp;
    SubNo     = 05;
    RunNo     = 1;
    NRuns     = 13;
    ThisRunNo = 1;
    Environment = 1;
    RealRun = 0;
    useEye = 0;
else
    RunNo = RunInfo(1);
    NRuns = RunInfo(2);
    ThisRunNo = RunInfo(3);
    global Environment
    global RealRun
end

tic;
RunType = 'pra';

try % open screen from here
    
    
    if  strcmp(RunType, 'task') && RealRun
        subInfo = GetParticipantInfo('age','gender','hand');
    end
    
    prm = prm_RIFT(RunType, Environment, useEye);
    Settings_General_MEG;

    %% Instruction
    
    % message (while loading images and waiting for scanner trigger)
    Inst{1} = ['In each trial, you will see two persons lifting their arm, one after another. \n\n',...
        'Please remember the position of their arms when you see them,\n\n later you will need to reproduce one of the two positions.\n\n',...
        'A cue of either "1" or "2" will follow to tell you which posture you need to reproduce.\n\n You can forget about the uncued one.\n\n',...
        'After a blank screen, the cued person will appear again.\n\n',...
        'You need to press upkey or downkey to move their arm to the previous position.\n\n',...
        'Then press Space to confirm. You have 7 seconds to respond.\n\n',...
        'If you don''t confirm in 7s, the next trial will start regardless',...
        'Please make sure to always stare at the fixation point at the center. \n\n'...
        'Press Space to start!'];
    Inst{2} = ['In each trial, you will see two bars, one after another, one on the left, one on the right \n\n',...
        'Please remember the position the bars when you see them,\n\n later you will need to reproduce one of the two positions.\n\n',...
        'A cue of either "1" or "2" will follow to tell you which bar you need to reproduce.\n\n You can forget about the uncued one.\n\n',...
        'After a blank screen, the cued bar will appear again.\n\n',...
        'You need to press upkey or downkey to move the bar to the previous position.\n\n',...
        'Then press Space to confirm. You have 7 seconds to respond\n\n',...
        'If you don''t confirm in 7s, the next trial will start regardless',...
        'Please make sure to always stare at the fixation point at the center. \n\n'...
        'Press Space to start!'];
    
    if ~mod(SubNo,2)
        Inst = fliplr(Inst);
    end
    DrawFormattedText(prm.w.Number, Inst{1}, 'center', 'center', prm.monitor.white);
    Screen('Flip', prm.w.Number);
    
    
    %% Set diary and record the current script
    
    LogDir = fullfile(RunDir, 'Logs');
    if ~exist(LogDir, 'dir')
        mkdir(LogDir);
    end
    
    RunStart = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
    diary(fullfile(LogDir, sprintf('Sub%02d_%s_%g_%s.txt', SubNo, RunType, ThisRunNo, RunStart)));
    
    
    ScriptName = {'Adjustment_withMask','Settings_General_MEG', 'RUN'};
    myCode = cell(length(ScriptName),1);
    for i = 1:length(ScriptName)
        fid = fopen([ScriptName{i} '.m'],'r'); % 'r': open file for reading
        myCode{i} = fscanf(fid,'%300c');
        fclose('all');
    end
    
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
    
    table = repmat(table, 3, 1);
    table.angle1 = table.angleTarget;
    table.angleProbe = randSamp(prm.fac.adjustRange, height(table), 'n');
    
    
    %% Get the texture for all the iamges
    
    imgTxts = zeros(length(prm.fac.figures), length(prm.fac.adjustRange), length(prm.fac.views));
    for i = 1:length(prm.fac.figures)
        for j = prm.fac.adjustRange(1): prm.fac.adjustRange(end)
            for k = 1:length(prm.fac.views)
                theimg = strcat('stimuli/Front-Upper-shoulder-middle/', prm.fac.figures{i}, '_Front_', num2str(j), prm.fac.views{k}, '.png');
                [X1, ~, alpha1] = imread(theimg);
                X1(:,:,4) = alpha1;
                imgTxts(i,j,k) =   Screen( 'MakeTexture', prm.w.Number, X1);
            end
        end
    end
    
    
    figureId = strcmp(table.figure1,prm.fac.figures{2})+1;
    table.img1Txt = arrayfun(@(x,y,z) imgTxts(x,y,z), figureId, table.angle1, table.view1);
    %table.img2Txt = arrayfun(@(x,y,z) imgTxts(x,y,z), 3-figureId, table.angle2, table.view2);
    
    nMask = 8;
    for theimg = 1:nMask
        maskFile = strcat('stimuli/masks/m', num2str(theimg), '.jpg');
        mask = imread(maskFile);
        maskTxts(theimg)  =   Screen( 'MakeTexture', prm.w.Number, mask);
    end
    
    %% randomize and then make a copy for bars, counterbalance the order
    table = table(randperm(height(table)),:);
    prm.fac.stimuli = ["Body","Bar"];
    table.stimulus = repmat("Body", height(table), 1);
    
    if  strcmp(RunType, 'pra')
        table = table(1: 5,:);
    end
    
    table2 = table;
    table2.stimulus = repmat("Bar", height(table), 1);
    
    if mod(SubNo,2)
        T = [table; table2];
    else
        T = [table2; table];
    end
    
    T.resp = zeros(height(T),1);
    T.rt   = NaN(height(T),1);
    T.subjectid = SubNo * ones(height(T),1);
    T = movevars(T, 'subjectid','Before', 1);
    
    %% Experimental precedures ============================================
    
    prm.N.trial = height(T);
    prm.N.block = 4*2;
    if  strcmp(RunType, 'pra')
        prm.N.block = 2;
    end
    prm.N.trialPerBlock = prm.N.trial/prm.N.block;
    T.block = reshape(repmat(1:prm.N.block, prm.N.trialPerBlock,1),[],1);
    if strcmp(RunType, 'pra')
        T.block = reshape(repmat([1,2],height(T)/2,1),[],1);
    end
    
    % temporal parameters
    
    prm.dur.fix    = 1200;
    prm.dur.image1 = 1200;
    prm.dur.mask   = 1000;
    prm.dur.delay  = 2400;
    prm.dur.wait   = 7000; % for response
    
    prm.Nfr        = structfun(@(x) round(x/1000 * prm.w.RefreshRate), prm.dur,'UniformOutput',0); % need to
    prm.SwitchFr  = cumsum(cell2mat(struct2cell(prm.Nfr)));
    prm.dur.trial = sum(cell2mat(struct2cell(prm.dur)));
    prm.time.task = height(T)*prm.dur.trial/1000/60;
    
    prm.time.theory = (height(T)*prm.dur.trial)/1000/60;
      
    
    % spatial parameters
    prm.img.W      = size(X1,2);
    prm.img.H      = size(X1,1);
    prm.img.WAng   = 8;
    prm.img.Scale  = visAng2pix(prm.img.WAng, 0, prm.monitor)/prm.img.W;
    prm.img.rect   = [0, 0, prm.img.W * prm.img.Scale, prm.img.H * prm.img.Scale];
    
    prm.img.offset = 139;  % vertical shift to center the shoulder at fixation
    prm.img.offPix = prm.img.offset * prm.img.Scale;
    
    prm.img.pos    = CenterRectOnPoint(prm.img.rect,  prm.w.Center(1), prm.w.Center(2)+prm.img.offPix); % for left-facing images: put on the right
    
    prm.bar.armLength = 216 * prm.img.Scale;
    prm.bar.color = 255 - 0.5 * prm.monitor.grey*[1,1,1];
    
    %% Scanner trigger ====================================================
    toc % get the time cost for code initiation before trigger
    if RealRun % if dummy scanner or lab, get trigger to start the run
        warning('Waiting for scanner trigger...');
        BitsiScanner.clearResponses();
        firstScan = 0;
        
        while firstScan == 0
            while BitsiScanner.numberOfResponses() == 0  % what does this do??
                WaitSecs(0.001);
            end
            [resp] = BitsiScanner.getResponse(0.001, true);
            if resp == 97
                firstScan = 1;
            end
        end
    else
        warning('Sham run: press Spacebar to continue.');
        waitforspace; waitfornokey;
    end
    
    % !!  GET THE TRIGGER TIME
    prm.time.trigger  = GetSecs;
    Screen('Flip', prm.w.Number);
    % Start (sham) scanning
    fprintf('\n SCANNER TRIGGER (run %g): %s \n', RunNo, datestr(now,'dd-mm-yyyy HH:MM:SS'));
    
    
    %% Delay period at the beginning
       
    WaitSecs(prm.time.before/1000 - prm.time.countdown); % during this period the text is still shown
    fprintf('Waiting for the task to begin in %g seconds...\n', prm.time.countdown);
    
    for n = 1:prm.time.countdown
        CountDownTxt = sprintf('Starting in %g seconds', prm.time.countdown - n);
        if n < prm.time.countdown
            DrawFormattedText(prm.w.Number, CountDownTxt, 'center', prm.w.Center(2) + 30, prm.monitor.white);
        else
            % show fixation at the last second
            Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.w.Center, 1);
        end
        Screen('Flip', prm.w.Number);
        WaitSecs(1);
    end
    
    
    
    %% real task ==========================================================
    
    vbl = Screen('Flip', prm.w.Number);
    
    for thetrial = 1: prm.N.trial
        %
        %         if T.targetId(thetrial)==1
        a = find(imgTxts==T.img1Txt(thetrial));
        lineAdjDir = (T.view1(thetrial)-1.5)*2;
        %         elseif T.targetId(thetrial)==2
        %             a = find(imgTxts==T.img2Txt(thetrial));
        %             lineAdjDir = (T.view2(thetrial)-1.5)*2;
        %         end
        [thefig, theang, theview] = ind2sub(size(imgTxts), a);
        theRange = squeeze(imgTxts(thefig, :, theview));
        
        line1Dir = (T.view1(thetrial)-1.5) *2;
        angle1Rad = T.angle1(thetrial)/180*pi;
        %        angle2Rad = T.angle2(thetrial)/180*pi;
        line1End = [prm.w.Center(1) + line1Dir *prm.bar.armLength* sin(angle1Rad), prm.w.Center(2) - prm.bar.armLength*cos(angle1Rad)];
        %        line2End = [prm.w.Center(1) - line1Dir *armLength* sin(angle2Rad), prm.w.Center(2) - armLength*cos(angle2Rad)];
        
        currentAng = T.angleProbe(thetrial);
        preRespT = vbl;
        history = currentAng;
        nf = 1;
        
        respFlag     =   0;
        if RealRun BitsiBB.clearResponses(); end
        
        while ~respFlag && nf <= prm.SwitchFr(end)
            if nf <= prm.SwitchFr(1) % draw dixation
                Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.w.Center, 1);
            elseif  nf <= prm.SwitchFr(2) % target 1
                if strcmp(T.stimulus(thetrial), 'Body')
                    Screen('DrawTexture', prm.w.Number, T.img1Txt(thetrial), [], prm.img.pos);
                else
                    Screen('DrawLine', prm.w.Number, prm.bar.color, prm.w.Center(1), prm.w.Center(2), line1End(1), line1End(2),10);
                end
                Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.w.Center, 1);
                
            elseif  nf <= prm.SwitchFr(3) % mask
                id = 1: prm.Nfr.mask/nMask : prm.Nfr.mask;
                theimg = max(find((nf- prm.SwitchFr(2) - id)>=0));
                Screen('DrawTexture', prm.w.Number, maskTxts(theimg), [], prm.img.pos);
                Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.w.Center, 1);
                
            elseif  nf <= prm.SwitchFr(4) % delay
                Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.w.Center, 1);
            elseif nf <= prm.SwitchFr(end) % adjustment
                if strcmp(T.stimulus(thetrial), 'Body')
                    Screen('DrawTexture', prm.w.Number, theRange(currentAng), [], prm.img.pos);
                else
                    lineAdjEnd = [prm.w.Center(1) + lineAdjDir *prm.bar.armLength* sin(currentAng/180*pi), prm.w.Center(2) - prm.bar.armLength*cos(currentAng/180*pi)];
                    Screen('DrawLine', prm.w.Number, prm.bar.color, prm.w.Center(1), prm.w.Center(2), lineAdjEnd(1), lineAdjEnd(2),10);
                end
                Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color2, prm.w.Center, 1);
            end
            
            Screen('DrawingFinished', prm.w.Number);
            vbl = Screen('Flip', prm.w.Number, vbl + .5 * prm.w.ifi);
            
            if nf == prm.SwitchFr(1) + 1
                T.T1Onset(thetrial) = vbl - prm.time.trigger;
                % elseif nf == prm.SwitchFr(3) + 1
                %  T.T2Onset(thetrial) = vbl - prm.time.trigger;
            elseif nf == prm.SwitchFr(end-1) + 1
                T.probeOnset(thetrial) = vbl - prm.time.trigger;
            end
            nf = nf +1; % after refreshing, count one more frame
            
            %% collect response
            if RealRun
                % retreive responses from button box
                [wkey, timeStamp] = BitsiBB.getResponse( timeout, true);
                if ismember( wkey, RespKeys) && ~respFlag && nf > prm.SwitchFr(7)  % only response after image 2 is presented is valid
                    respFlag = 1;
                    T.resp(thetrial)   =   find(wkey==RespKeys);          % 1 or 2
                    T.rt(thetrial)     =   timeStamp - T.T2Onset(thetrial) - prm.time.trigger;    % time from scanner trigger
                    fprintf('Responded %g, should be %g. \n', T.resp(thetrial), T.correctKey(thetrial));
                end
            end
            
            % also check keyboard, no matter realRun or not.
            [keydown,respT,keyCode] = KbCheck;
            if keydown
                if keyCode(KbName('ESCAPE')) %check quit
                    fprintf('\n\nExperiment terminated at %s, diary closed...\n', datestr(now));
                    if RealRun
                        close(BitsiScanner);
                        fprintf('All serial ports closed.\n');
                    end
                    
                    ShowCursor;
                    diary off
                    clear temprun
                    sca; return
                elseif ~respFlag && nf > prm.SwitchFr(end-1) %collect response only after delay finished
                    if keyCode(KbName('UpArrow')) && respT - preRespT > 1*prm.w.ifi % to control of speed of adjustment
                        currentAng = currentAng - 1;
                        currentAng = max(currentAng, min(prm.fac.adjustRange));
                        history = [history, currentAng];
                        if length(history)==2
                            T.adjustOnset(thetrial) = respT- prm.time.trigger;
                        end
                    elseif keyCode(KbName('DownArrow')) && respT - preRespT > 1*prm.w.ifi
                        currentAng = currentAng + 1;
                        currentAng = min(currentAng, max(prm.fac.adjustRange));
                        history = [history, currentAng];
                        if length(history)==2
                            T.adjustOnset(thetrial) = respT- prm.time.trigger;
                        end
                    elseif keyCode(KbName('Space'))
                        respFlag = 1;
                        T.resp(thetrial)   =   currentAng;
                        T.rt(thetrial)     =   respT - T.probeOnset(thetrial) - prm.time.trigger;
                        T.adjustHistory{thetrial} = history;
                        fprintf('Responded %g, should be %g.\n', T.resp(thetrial), T.angleTarget(thetrial));
                    end % end of kbcheck
                    preRespT = respT;
                end
            end
            
        end
        
        
        %% break
        if thetrial == prm.N.trial/2 % start the new stimulus
            msg = sprintf(['Well done! You completed half of the experiment!\n\n', ...
                'You can take a break. \n\n\n Press Space to start the next block']);
            DrawFormattedText(prm.w.Number, msg, 'center', 'center', prm.monitor.black);
            Screen('DrawingFinished', prm.w.Number);
            vbl = Screen('Flip', prm.w.Number, vbl + .5 * prm.w.ifi);
            WaitSecs(1); waitforspace; waitfornokey;
            
            
            DrawFormattedText(prm.w.Number, Inst{2}, 'center', prm.w.Center(2) - 30, prm.monitor.black);
            Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.w.Center, 1);
            Screen('DrawingFinished', prm.w.Number);
            vbl = Screen('Flip', prm.w.Number, vbl + .5 * prm.w.ifi);
            waitforspace; waitfornokey;
            
            % Otherwise, a normal break
        elseif thetrial>1 && mod(thetrial, prm.N.trialPerBlock)==0
            blockPrompt= sprintf('You finished %d blocks! \n\n Take a break. \n\n When you are ready, press SPACE to start!', T.block(thetrial)-1);
            DrawFormattedText(prm.w.Number, blockPrompt, 'center', 'center', prm.monitor.black);
            vbl = Screen('Flip', prm.w.Number, vbl + .5 * prm.w.ifi);
            waitforspace; waitfornokey;

        end
    end
  % end of task
    
    
    
    %% report participants'results
    T.Error = T.resp - T.angleTarget; % no reponse for MVPA trials are also correct
    figure;
    for thestim = 1:2
        trialid = strcmp(T.stimulus,prm.fac.stimuli{thestim}) ;
        meanError(thestim) = mean(T.Error(trialid)); % only calculate acc for RS trials
        STD(thestim) = std(T.Error(trialid));
    end
    bar([1,3], meanError); hold on
    xticklabels(prm.fac.stimuli);
    errorbar([1,3], meanError, STD)
    fprintf('((( mean Error for Body: %g, for Bar = %g )))\n', meanError(1), meanError(2)); % to print a % inside fprintf, you need %%
    
    %% Save 1. all data 2. resulttable=====================================
    
    fileName1 = fullfile(RunDir, sprintf('Sub%02d_%s_%g.mat', SubNo, RunType, ThisRunNo));
    if ~exist('subInfo', 'var')
        save(fileName1, 'prm','T')
    else
        save(fileName1,'prm','subInfo','T');
    end
    fileName2 = fullfile(RunDir, sprintf('ResultTable_Sub%02d_%s_%g.csv', SubNo, RunType, ThisRunNo));
    writetable(T, fileName2);
    
    
    %% EXIT and clean up===================================================
        % prm.monitor PTB performance    
    prm.dur.trialActual = 1000 * (T.probeOnset - T.T1Onset); % get the actual duraction for each trial
    prm.dur.trialActual = prm.dur.trialActual(prm.dur.trialActual> 0 & prm.dur.trialActual< prm.time.rest);
    TrialDur = round(nanmean(prm.dur.trialActual));
    TrialStd = round(nanstd( prm.dur.trialActual));
    fprintf('((( Trial Duration %g/%g ms, SD = %g ms )))\n', TrialDur, prm.dur.trial, TrialStd);
    %  figure;  scatter(ones(size(prm.dur.trialActual)), prm.dur.trialActual)
    prm.time.total = toc/60;
    
    
    
    
    
    if RealRun
        close(BitsiScanner);
        wmsg = 'scanner and response box';
    else
        wmsg = 'virtual response box';
    end
    fprintf('\nSerial ports (%s) CLOSED at %s. \n\n', wmsg, datestr(now));
    
    prm.time.actual  =  round( GetSecs - prm.time.trigger)/60;
    fprintf('RUN DURATION: %g mins (expected duration: %g mins). \n\n', prm.time.actual, prm.time.theory);
    

    % clean log
    fprintf('Diary closed (%s)\n\n', datestr(now));
    diary off
    
    % clean memory
    %  cleans up the Java heap from memory leaks, preventing the infamous
    %  Java OutOfMemory exception. https://nl.mathworks.com/matlabcentral/fileexchange/36757-java-heap-cleaner
    try  % use this to avoid a cessation from error
        jheapcl;
    end
    
    
    sca
    ShowCursor;
    Priority(0);
    toc;
    
catch ME
    sca
    ShowCursor;
    Priority(0);
    toc;
    rethrow(ME);
end
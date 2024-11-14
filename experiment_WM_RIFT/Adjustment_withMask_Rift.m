function Adjustment_withMask_Rift(SubNo, RunInfo)

%{
the single-target design with flickering mask

This is the most updated version. can be easily adapted into an retro-cue
paradigm 
%}


if ~exist( 'SubNo', 'var') % this will run if the "function" line is commented
    sca; clear all; close all; clc; commandwindow;

    addpath(genpath('Functions'));
    StartUp;
    SubNo     = 05;
    RunNo     = 1;
    NRuns     = 13;
    ThisRunNo = 1;
    RunType  = 'pra';
    Environment = 1;
    RealRun = 0;
    useEye = 0;
    prm = prm_RIFT(RunType, RealRun, Environment, useEye);
    
else
    RunNo = RunInfo(1);
    NRuns = RunInfo(2);
    ThisRunNo = RunInfo(3);

end

tic;

try % open screen from here

    
    Settings_General_MEG;

    %% Instruction
    
    % message (while loading images and waiting for scanner trigger)
    Inst{1} = ['In each trial, you will see two persons lifting their arm, one after another. \n\n',...
        'Please remember the position of their arms when you see them,\n\n later you will need to reproduce one of the two positions.\n\n',...
        'A cue of either "1" or "2" will follow to tell you which posture you need to reproduce.\n\n You can forget about the uncued one.\n\n',...
        'After a blank screen, the cued person will appear again.\n\n',...
        'You need to press upkey or downkey to move their arm to the previous position.\n\n',...
        'Then press Space to confirm. You have 7 seconds to respond.\n\n',...
        'If you don''t confirm in 7s, the next trial will start regardless\n\n',...
        'Please make sure to always stare at the fixation point at the center. \n\n'...
        'Press Space to start!'];
    Inst{2} = ['In each trial, you will see two bars, one after another, one on the left, one on the right \n\n',...
        'Please remember the position the bars when you see them,\n\n later you will need to reproduce one of the two positions.\n\n',...
        'A cue of either "1" or "2" will follow to tell you which bar you need to reproduce.\n\n You can forget about the uncued one.\n\n',...
        'After a blank screen, the cued bar will appear again.\n\n',...
        'You need to press upkey or downkey to move the bar to the previous position.\n\n',...
        'Then press Space to confirm. You have 7 seconds to respond\n\n',...
        'If you don''t confirm in 7s, the next trial will start regardless\n\n',...
        'Please make sure to always stare at the fixation point at the center. \n\n'...
        'Press Space to start!'];
    
    if ~mod(SubNo,2)
        Inst = fliplr(Inst);
    end
    DrawFormattedText(prm.w.Number, Inst{1}, 'center', 'center', prm.monitor.white);
    Screen('Flip', prm.w.Number);
    
    

    
    %% Constructing all the trials=========================================
    prm.fac.adjustRange = (25:65)';
    prm.fac.targetRange = (30:6:60)';
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
    
    if  strcmp(prm.exp.RunType , 'pra')
        table = table(1: prm.N.praTrial,:);
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
    if  strcmp(prm.exp.RunType , 'pra')
        prm.N.block = 2;
    end
    prm.N.trialPerBlock = prm.N.trial/prm.N.block;
    T.block = reshape(repmat(1:prm.N.block, prm.N.trialPerBlock,1),[],1);
    if strcmp(prm.exp.RunType , 'pra')
        T.block = reshape(repmat([1,2],height(T)/2,1),[],1);
    end
    

    % spatial parameters
    prm.img.W      = size(X1,2);
    prm.img.H      = size(X1,1);
    prm.img.WAng   = 8;
    prm.img.scale  = visAng2pix(prm.img.WAng, 0, prm.monitor)/prm.img.W;
    prm.img.rect   = [0, 0, prm.img.W * prm.img.scale, prm.img.H * prm.img.scale];
    
    prm.img.offset = 139;  % vertical shift to center the shoulder at fixation
    prm.img.offPix = prm.img.offset * prm.img.scale;
    
    prm.img.pos    = CenterRectOnPoint(prm.img.rect,  prm.w.Center(1), prm.w.Center(2)+prm.img.offPix); % for left-facing images: put on the right
    
    prm.bar.armLength = 216 * prm.img.scale;
    prm.bar.color = 255 - 0.5 * prm.monitor.grey*[1,1,1]; % a light grey
    
    % make tagging fans' texture 
    [blackTxts, whiteTxts] = make_tag_texture(prm);    
    
    %% press to start ====================================================
    
    warning('In instruction. Press Spacebar to continue.');
    waitforspace; waitfornokey;
    Screen('Flip', prm.w.Number);

    if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.ExpStart); end    
    prm.time.expStart  = GetSecs;
    % Start (sham) scanning
    fprintf('\n Starting the task: (run %g): %s \n', RunNo, datestr(now,'dd-mm-yyyy HH:MM:SS'));
    
    
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
        if prm.exp.eyelink_live
            Eyelink('Message', 'Trial %d', thetrial);
        end
    
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
                for thestrip = 1:prm.tag.nStrip
                  this_tag = prm.tag.tag_sigs(thestrip, nf);
                  Screen('DrawTexture', prm.w.Number, blackTxts(thestrip, theview), [],[],[],[], []); % draw the black line normally
                  Screen('DrawTexture', prm.w.Number, whiteTxts(thestrip, theview), [],[],[],[], this_tag); % use transparency to change
                end
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
                for thestrip = 1:prm.tag.nStrip
                  this_tag = prm.tag.tag_sigs(thestrip, nf);
                  Screen('DrawTexture', prm.w.Number, blackTxts(thestrip, theview), [],[],[],[], []); % draw the black line normally
                  Screen('DrawTexture', prm.w.Number, whiteTxts(thestrip, theview), [],[],[],[], this_tag); % use transparency to change
                end
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

            if nf== 1
                if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.FixStart); end
                
            elseif     nf == prm.SwitchFr(1) + 1
                T.T1Onset(thetrial) = vbl - prm.time.expStart;
                if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.TargetStart); end
            elseif     nf == prm.SwitchFr(1) + 1
                T.T1Onset(thetrial) = vbl - prm.time.expStart;
                if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.TargetStart); end                
            elseif nf == prm.SwitchFr(end-1) + 1
                T.probeOnset(thetrial) = vbl - prm.time.expStart;
                if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.probeStart); end
            end
            nf = nf +1; % after refreshing, count one more frame
            
            %% collect response
            if RealRun
                % retreive responses from button box
                [wkey, timeStamp] = BitsiBB.getResponse( timeout, true);
                if ismember( wkey, RespKeys) && ~respFlag && nf > prm.SwitchFr(7)  % only response after image 2 is presented is valid
                    respFlag = 1;
                    T.resp(thetrial)   =   find(wkey==RespKeys);          % 1 or 2
                    T.rt(thetrial)     =   timeStamp - T.T2Onset(thetrial) - prm.time.expStart;    % time from scanner trigger
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
                            T.adjustOnset(thetrial) = respT- prm.time.expStart;
                        end
                    elseif keyCode(KbName('DownArrow')) && respT - preRespT > 1*prm.w.ifi
                        currentAng = currentAng + 1;
                        currentAng = min(currentAng, max(prm.fac.adjustRange));
                        history = [history, currentAng];
                        if length(history)==2
                            T.adjustOnset(thetrial) = respT- prm.time.expStart;
                        end
                    elseif keyCode(KbName('Space'))
                        respFlag = 1;
                        if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.Response); end
                        T.resp(thetrial)   =   currentAng;
                        T.rt(thetrial)     =   respT - T.probeOnset(thetrial) - prm.time.expStart;
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
        elseif thetrial< prm.N.trial && mod(thetrial, prm.N.trialPerBlock)==0
            % first thing first, change the projector mode
            if RealRun && prm.exp.eyelink_live
                Datapixx('SetPropixxDlpSequenceProgram', 0); % 2 for 480, 5 for 1440 Hz, 0 for normal
                Datapixx('RegWrRd');
                
                blockPrompt= sprintf('You finished %d blocks! \n\n We are going to calibrate the eye-tracker again', T.block(thetrial)-1);
            else                        
            blockPrompt= sprintf('You finished %d blocks! \n\n Take a break. \n\n When you are ready, press SPACE to start!', T.block(thetrial)-1);
            end
            DrawFormattedText(prm.w.Number, blockPrompt, 'center', 'center', prm.monitor.black);
            vbl = Screen('Flip', prm.w.Number, vbl + .5 * prm.w.ifi);
            waitforspace; waitfornokey;
            if prm.exp.eyelink_live
                Eyelink('Message', 'ENDEXP');
                Eyelink('StopRecording'); % Stop tracker recording
                Eyelink('command', 'set_idle_mode');
                WaitSecs(0.5);
                Eyelink('CloseFile'); % Close EDF file on Host PC
                %% Transfer the edf file
                transfer_edf(prm); % try transfering the edf file
                Eyelink('Shutdown');
            end
        end
    end
    % end of task
    if RealRun
        prm.trigger.btsi.sendTrigger(prm.trigger.ExpEnd);
        Datapixx('SetPropixxDlpSequenceProgram', 0); % set the projector to normal mode
        Datapixx('RegWrRd')
    end

    
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
    
    fileName1 = fullfile(RunDir, sprintf('Sub%02d_%s_%g.mat', SubNo, prm.exp.RunType , ThisRunNo));
    if ~exist('subInfo', 'var')
        save(fileName1, 'prm','T')
    else
        save(fileName1,'prm','subInfo','T');
    end
    fileName2 = fullfile(RunDir, sprintf('ResultTable_Sub%02d_%s_%g.csv', SubNo, prm.exp.RunType , ThisRunNo));
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
    
    prm.time.actual  =  round( GetSecs - prm.time.expStart)/60;
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
%function Body_RunRS(SubNo, RunInfo, DataDir)

%{
the retro-cue paradigm


%}
clear

if ~exist( 'SubNo', 'var') % this will run if the "function" line is commented
    
    addpath(genpath('Functions'));
    StartUp;
    SubNo     = 05;
    RunNo     = 1;
    NRuns     = 13;
    ThisRunNo = 1;
    DataDir   = 'Data_Demo';
    Environment = 1;
    RealRun = 0;
else
    RunNo = RunInfo(1);
    NRuns = RunInfo(2);
    ThisRunNo = RunInfo(3);
    global Environment
    global RealRun
end

tic;
RunType = 'Test';

try % open screen from here
    

    if  strcmp(RunType, 'task') &&RealRun
        info = GetParticipantInfo('age','gender','hand');
    end
    
    Settings_General;
    
    %% Instruction
    
    % message (while loading images and waiting for scanner trigger)
    Inst{1} = ['In each trial, you will see two persons lifting their arm, one after another. \n\n',...
        'Please remember the position of their arms when you see them,\n\n later you will need to reproduce one of the two positions.\n\n',...
        'A cue of either "1" or "2" will follow to tell you which posture you need to reproduce.\n\n You can forget about the uncued one.\n\n',...
        'After a blank screen, the cued person will appear again.\n\n',...
        'You need to press upkey or downkey to move their arm to the previous position.\n\n',...
        'Then press Space to confirm. You have 7 seconds to respond\n\n',...
        'Please make sure to always stare at the fixation point at the center. \n\n'...
        'Press Space to start!'];
   Inst{2} = ['In each trial, you will see two bars, one after another, one on the left, one on the right \n\n',...
        'Please remember the position the bars when you see them,\n\n later you will need to reproduce one of the two positions.\n\n',...
        'A cue of either "1" or "2" will follow to tell you which bar you need to reproduce.\n\n You can forget about the uncued one.\n\n',...
        'After a blank screen, the cued bar will appear again.\n\n',...
        'You need to press upkey or downkey to move the bar to the previous position.\n\n',...
        'Then press Space to confirm. You have 7 seconds to respond\n\n',...
        'Please make sure to always stare at the fixation point at the center. \n\n'...
        'Press Space to start!'];
    if ~mod(SubNo,2)
        Inst = fliplr(Inst);
    end
    
    DrawFormattedText(w.Number, Inst{1}, 'center', 'center', White);
    Screen('Flip', w.Number);
    
    
    %% Set diary and record the current script
    
    LogDir = fullfile(RunDir, 'Logs');
    if ~exist(LogDir, 'dir')
        mkdir(LogDir);
    end
    
    RunStart = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
    diary(fullfile(LogDir, sprintf('Sub%02d_%s_%g_%s.txt', SubNo, RunType, ThisRunNo, RunStart)));
    
    
    ScriptName = {'UpDownTask_retroCue_flip','Settings_General', 'RUN'};
    myCode = cell(length(ScriptName),1);
    for i = 1:length(ScriptName)
        fid = fopen([ScriptName{i} '.m'],'r'); % 'r': open file for reading
        myCode{i} = fscanf(fid,'%300c');
        fclose('all');
    end
    
    %% Constructing all the trials=========================================
    adjustRange = (25:65)';
    targetRange = (30:5:60)';
    views = {'L', 'R'}; % 1: facing left (presented onthe right)
    figures = {'Fe'; 'Ma'};
    targetIds = [1;2];
    
    % Generate all combinations using ndgrid
    [t, d, v, f] = ndgrid(targetIds, targetRange, 1:length(views),  figures);
    
    % Combine into a cell array
    combinations = [num2cell(t(:)), num2cell(d(:)), num2cell(v(:)),  f(:)];
    
    % Create the table
    table = cell2table(combinations, ...
        'VariableNames', {'targetId','angleTarget', 'view1',  'figure1', });
    table.view2 = 3 - table.view1;
    table.figure2 = figures(strcmp(table.figure1, figures{1}) + 1);
    
    table = repmat(table, 3, 1);
    
    % angle 1 ang angle 2 have to be fully independent to ensure they can be analyzed separately
    % first column: cued, second column: uncued
    target_nontarget = [table.angleTarget, randSamp(targetRange, height(table), 'n')];
    table.angle1 = target_nontarget(sub2ind(size(target_nontarget), (1:height(table))', table.targetId));
    table.angle2 = target_nontarget(sub2ind(size(target_nontarget), (1:height(table))', 3-table.targetId));
    table.angleProbe = randSamp(adjustRange, height(table), 'n');
    
    scatter(table.angleTarget, target_nontarget(:,2))
    
    %% Get the texture for all the iamges

    imgTxts = zeros(length(figures), length(adjustRange), length(views));
    for i = 1:length(figures)
        for j = adjustRange(1): adjustRange(end)
            for k = 1:length(views)
                theimg = strcat('Stimuli/Front-Upper-shoulder-middle/', figures{i}, '_Front_', num2str(j), views{k}, '.png');
                [X1, ~, alpha1] = imread(theimg);
                X1(:,:,4) = alpha1;
                imgTxts(i,j,k) =   Screen( 'MakeTexture', w.Number, X1);
            end
        end
    end
    
    
    figureId = strcmp(table.figure1,figures{2})+1;
    table.img1Txt = arrayfun(@(x,y,z) imgTxts(x,y,z), figureId, table.angle1, table.view1);
    table.img2Txt = arrayfun(@(x,y,z) imgTxts(x,y,z), 3-figureId, table.angle2, table.view2);
    
    
    maskFile = strcat('Stimuli/masks/m', num2str(2), '.jpg');
    mask = imread(maskFile);
    maskTxt  =   Screen( 'MakeTexture', w.Number, mask);
    
    % randomize and then make a copy 
    table = table(randperm(height(table)),:);
    stimuli = {"Body","Bar"};
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
    
    % temporal parameters
    %iti   = 1800 + [-200, 0, 200];
    %T.iti = randSamp(iti, height(T), 'n');
    dur.fix    = 1200;
    dur.image1 = 400;
    dur.blank1 = 500;
    dur.image2 = 400;
    dur.blank2 = 500;
    dur.cue    = 400;
    dur.delay  = 1000;
    dur.wait   = 7000; % for response 
    
    Nfr = structfun(@(x) round(x/1000 * w.RefreshRate), dur,'UniformOutput',0); % need to
    SwitchFr = cumsum(cell2mat(struct2cell(Nfr)));
    dur.trial = sum(cell2mat(struct2cell(dur)));
    time.task = height(T)*dur.trial/1000/60;
    time.rest = 12 * 1000; % try to use ms for all time variables
    
    Nfr.rest = round(time.rest/1000 * w.RefreshRate);
    time.theory = (height(T)*dur.trial + time.before *2 + time. rest)/1000/60;
    
    NTrial = height(T);
    
    
    % spatial parameters
    img.W      = size(X1,2);
    img.H      = size(X1,1);
    img.WAng   = 8;
    img.Scale = visAng2pix(img.WAng, 0, monitor)/img.W;
    img.rect = [0, 0, img.W * img.Scale, img.H * img.Scale];
    
    img.offset = 139;  % vertical shift to center the shoulder at fixation
    img.offPix = img.offset * img.Scale;
    
    img.pos = CenterRectOnPoint( img.rect,  w.Center(1), w.Center(2)+img.offPix); % for left-facing images: put on the right
    
    armLength = 216 * img.Scale;
    
    bar.color = 255 - 0.5 * Grey*[1,1,1];
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
    time.trigger  = GetSecs;
    Screen('Flip', w.Number);
    % Start (sham) scanning
    fprintf('\n SCANNER TRIGGER (run %g): %s \n', RunNo, datestr(now,'dd-mm-yyyy HH:MM:SS'));
    
    
    %% Delay period at the beginning
    
    
    WaitSecs(time.before/1000 - time.countdown); % during this period the text is still shown
    fprintf('Waiting for the task to begin in %g seconds...\n', time.countdown);
    
    for n = 1:time.countdown
        CountDownTxt = sprintf('Starting in %g seconds', time.countdown - n);
        if n < time.countdown
            DrawFormattedText(w.Number, CountDownTxt, 'center', w.Center(2) + 30, White);
        else
            % show fixation at the last second
            Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
        end
        Screen('Flip', w.Number);
        WaitSecs(1);
    end
    
    
    
    %% real task ==========================================================
    
    vbl = Screen('Flip', w.Number);
    
    for thetrial = 1: NTrial
        
        if T.targetId(thetrial)==1
            a = find(imgTxts==T.img1Txt(thetrial));
            lineAdjDir = (T.view1(thetrial)-1.5)*2;
        elseif T.targetId(thetrial)==2
            a = find(imgTxts==T.img2Txt(thetrial));
            lineAdjDir = (T.view2(thetrial)-1.5)*2;
        end
        [thefig, theang, theview] = ind2sub(size(imgTxts), a);
        theRange = squeeze(imgTxts(thefig, :, theview));
        
        line1Dir = (T.view1(thetrial)-1.5) *2;
        angle1Rad = T.angle1(thetrial)/180*pi;
        angle2Rad = T.angle2(thetrial)/180*pi;
        line1End = [w.Center(1) + line1Dir *armLength* sin(angle1Rad), w.Center(2) - armLength*cos(angle1Rad)];
        line2End = [w.Center(1) - line1Dir *armLength* sin(angle2Rad), w.Center(2) - armLength*cos(angle2Rad)];
        
        currentAng = T.angleProbe(thetrial);
        preRespT = vbl;
        history = currentAng;
        nf = 1;
        
        respFlag     =   0;
        if RealRun BitsiBB.clearResponses(); end
        
        while ~respFlag && nf <= SwitchFr(end)
            if nf <= SwitchFr(1) % draw dixation
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
            elseif  nf <= SwitchFr(2) % target 1
                if strcmp(T.stimulus(thetrial), 'Body')
                    Screen('DrawTexture', w.Number, T.img1Txt(thetrial), [], img.pos);
                else
                    Screen('DrawLine', w.Number, bar.color, w.Center(1), w.Center(2), line1End(1), line1End(2),10);
                end
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
                
            elseif  nf <= SwitchFr(3) % blank 1
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
            elseif  nf <= SwitchFr(4) % target 2
                if strcmp(T.stimulus(thetrial), 'Body')                    
                    Screen( 'DrawTexture', w.Number, T.img2Txt(thetrial), [], img.pos);
                else
                    Screen('DrawLine', w.Number, bar.color, w.Center(1), w.Center(2), line2End(1), line2End(2),10);
                end
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
            elseif  nf <= SwitchFr(5) % blank 2
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
            elseif nf <= SwitchFr(6) % cue
                DrawFormattedText(w.Number, num2str(T.targetId(thetrial)), 'center', 'center', text.Color);
            elseif  nf <= SwitchFr(7) % blank 3
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
            elseif nf <= SwitchFr(8) % adjustment
                if strcmp(T.stimulus(thetrial), 'Body')                    
                    Screen('DrawTexture', w.Number, theRange(currentAng), [], img.pos);
                else
                    lineAdjEnd = [w.Center(1) + lineAdjDir *armLength* sin(currentAng/180*pi), w.Center(2) - armLength*cos(currentAng/180*pi)];
                    Screen('DrawLine', w.Number, bar.color, w.Center(1), w.Center(2), lineAdjEnd(1), lineAdjEnd(2),10);
                end
                Screen('DrawDots', w.Number, [0; 0], fix.Size, [0,0,0], w.Center, 1);
            end
            
            Screen('DrawingFinished', w.Number);
            vbl = Screen('Flip', w.Number, vbl + .5 * w.ifi);
            
            if nf == SwitchFr(1) + 1
                T.T1Onset(thetrial) = vbl - time.trigger;
            elseif nf == SwitchFr(3) + 1
                T.T2Onset(thetrial) = vbl - time.trigger;
            elseif nf == SwitchFr(end) + 1
                T.probeOnset(thetrial) = vbl - time.trigger;
            end
            nf = nf +1; % after refreshing, count one more frame
            
            %% collect response
            if RealRun
                % retreive responses from button box
                [wkey, timeStamp] = BitsiBB.getResponse( timeout, true);
                if ismember( wkey, RespKeys) && ~respFlag && nf > SwitchFr(7)  % only response after image 2 is presented is valid
                    respFlag = 1;
                    T.resp(thetrial)   =   find(wkey==RespKeys);          % 1 or 2
                    T.rt(thetrial)     =   timeStamp - T.T2Onset(thetrial) - time.trigger;    % time from scanner trigger
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
                    close(BitsiBB);
                    ShowCursor;
                    diary off
                    clear temprun
                    sca; return
                elseif ~respFlag && nf > SwitchFr(7) %collect response only after delay finished
                    if keyCode(KbName('UpArrow')) && respT - preRespT > 1*w.ifi % to control of speed of adjustment
                        currentAng = currentAng - 1;                     
                        currentAng = max(currentAng, min(adjustRange));
                        history = [history, currentAng];
                        if length(history)==2
                            T.adjustOnset(thetrial) = respT- time.trigger;
                        end
                    elseif keyCode(KbName('DownArrow')) && respT - preRespT > 1*w.ifi 
                        currentAng = currentAng + 1;
                        currentAng = min(currentAng, max(adjustRange));
                        history = [history, currentAng];
                        if length(history)==2
                            T.adjustOnset(thetrial) = respT- time.trigger;
                        end
                    elseif keyCode(KbName('Space'))
                        respFlag = 1;
                        T.resp(thetrial)   =   currentAng;
                        T.rt(thetrial)     =   respT - T.probeOnset(thetrial) - time.trigger;
                        T.adjustHistory{thetrial} = history;
                        fprintf('Responded %g, should be %g.\n', T.resp(thetrial), T.angleTarget(thetrial));
                    end % end of kbcheck
                    preRespT = respT;
                end
            end
            
        end
        
        if thetrial == NTrial/2
                msg = sprintf(['Well done! You completed this block\n\n', ...
        'You can take a break. \n\n\n Press Space to start the next block']);
            DrawFormattedText(w.Number, msg, 'center', 'center', Black);
            Screen('DrawingFinished', w.Number);
            vbl = Screen('Flip', w.Number, vbl + .5 * w.ifi);
            WaitSecs(3); waitforspace; waitfornokey;
            
                
                DrawFormattedText(w.Number, Inst{2}, 'center', w.Center(2) - 30, Black);
                Screen('DrawDots', w.Number, [0; 0], fix.Size, fix.Col, w.Center, 1);
                Screen('DrawingFinished', w.Number);
                vbl = Screen('Flip', w.Number, vbl + .5 * w.ifi);
                waitforspace; waitfornokey;
            
        end
    end
    
    %% report participants'results
    T.Error = T.resp - T.angleTarget; % no reponse for MVPA trials are also correct
    figure;
    for thestim = 1:2
            trialid = strcmp(T.stimulus,stimuli{thestim}) ;
            meanError(thestim) = mean(T.Error(trialid)); % only calculate acc for RS trials
            STD(thestim) = std(T.Error(trialid));
    end
    bar([1,3], meanError); hold on
xticklabels(stimuli);
errorbar([1,3], meanError, STD)
   fprintf('((( mean Error for Body: %g, for Bar = %g )))\n', meanError(1), meanError(2)); % to print a % inside fprintf, you need %%
    
    
    % Monitor PTB performance
    
    dur.trialActual = 1000 * (T.probeOnset - T.T1Onset); % get the actual duraction for each trial
    dur.trialActual = dur.trialActual(dur.trialActual> 0 & dur.trialActual< time.rest);
    TrialDur = round(nanmean(dur.trialActual));
    TrialStd = round(nanstd( dur.trialActual));
    fprintf('((( Trial Duration %g/%g ms, SD = %g ms )))\n', TrialDur, dur.trial, TrialStd);
     %  figure;  scatter(ones(size(dur.trialActual)), dur.trialActual)
    time.total = toc/60;
    
    
    %% EXIT and clean up===================================================

    close(BitsiBB);
    
    
    if RealRun
        close(BitsiScanner);
        wmsg = 'scanner and response box';
    else
        wmsg = 'virtual response box';
    end
    fprintf('\nSerial ports (%s) CLOSED at %s. \n\n', wmsg, datestr(now));
    
    time.actual  =  round( GetSecs - time.trigger)/60;
    fprintf('RUN DURATION: %g mins (expected duration: %g mins). \n\n', time.actual, time.theory);
    
    %% Save 1. all data 2. resulttable=====================================
    
    fileName1 = fullfile(RunDir, sprintf('Sub%02d_%s_%g.mat', SubNo,RunType, ThisRunNo));
    save(fileName1);
    fileName2 = fullfile(RunDir, sprintf('ResultTable_Sub%02d_%s_%g.csv', SubNo,RunType, ThisRunNo));
    writetable(T, fileName2);
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
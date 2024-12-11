%function [prm, T] = Adjustment_withMask_RIFT(prm)

%{
the single-target design with flickering mask

1. make sure the tags are boardering instead of overlapping with the

%}


if ~exist( 'prm', 'var') % this will run if the "function" line is commented
    sca; clear; close all; clc; commandwindow;

    addpath(genpath('Functions'));
    StartUp;
    SubNo       = 98;
    Environment = 2;
    RealRun     = 1;
    useEye      = 0;
    RunType     = 'pra';
    prm = prm_RIFT(SubNo, RunType, RealRun, Environment, useEye);

else
    RealRun      = prm.exp.RealRun;
    useEye       = prm.exp.useEye;
    Environment  = prm.exp.Environment;
    SubNo        = prm.exp.SubNo;
    RunType      = prm.exp.RunType;
end

T = prm.exp.T;
tic;

try % open screen from here


    Settings_General_MEG;

    %% setup eye-tracker

    % eye-link calibaration
    if useEye
        Datapixx('Open');
        Datapixx('SetPropixxDlpSequenceProgram', 0); % 2 for 480, 5 for 1440 Hz, 0 for normal
        Datapixx('RegWrRd');

        prm = eyelink_init(prm, prm.subInfo); % initialize the eye tracker parameters
        EyelinkEnterSetup(prm.eyelink);% enter the set up for the eyelink
        Screen('FillRect',prm.w.Number,127); % make it grey again
    end


    % start the eye tracker
    if useEye
        Eyelink('StartRecording');
        WaitSecs(0.050);
        Eyelink('Message', 'STARTEXP');
    end

    % message (while loading images and waiting for scanner trigger)
    instruction_screen(prm, prm.exp.inst);


    %% Get the texture for all the iamges

    [imgTxts, maskTxts, prm] = make_stimuli_texture(prm);
    % make tagging fans' texture
    [blackTxts, whiteTxts, prm] = make_tag_texture(prm);

    figureId = strcmp(T.figure1, prm.fac.figures{2})+1;
    T.img1Txt = arrayfun(@(x,y,z) imgTxts(x,y,z), figureId, T.angleTarget, T.view1);


    %diode tracking square
    sz = prm.diode_track.SizeInPxl;
    prm.diode_track.texture = Screen('OpenOffscreenWindow', prm.w.Number, [0 0 0 0], [0 0 prm.w.Center(1) prm.w.Center(2)]);
    Screen('FillRect', prm.diode_track.texture, 255, [0 prm.w.Center(2)-sz sz prm.w.Center(2)]); % on the lower-left corner



    % determine the positions in four quadrants

    for thequa = 1:4
        subcenter = MakeOffsetRect(prm.w, [0,0], 0, 0, thequa);
        prm.tag.subcenters(thequa,:)  = subcenter(1:2);

        prm.img.pos(thequa, 1:4)    = MakeOffsetRect(prm.w, prm.img.presentedSize, 0, 0, thequa);
        prm.tag.pos(thequa, 1:4, 1)    = MakeOffsetRect(prm.w, prm.tag.size, -prm.tag.size(1)/2, -prm.tag.size(1)/2, thequa);
        rightRect    = MakeOffsetRect(prm.w, prm.tag.size, prm.tag.size(1)/2,  -prm.tag.size(1)/2, thequa);
        prm.tag.pos(thequa, 1:4, 2)    = rightRect([3, 2, 1, 4]); %flipping of tag images to the right is achieved by flipping destinationRect
        prm.diode_track.positions(thequa,:) = MakeOffsetRect(prm.w, prm.w.Center, 0, 0, thequa);

        for j = prm.fac.adjustRange(1): prm.fac.adjustRange(end)
            for k = 1:length(prm.fac.views)
                lineDir = (k-1.5) *2;
                lineEnds(j,k,thequa,:) = [prm.tag.subcenters(thequa, 1) + lineDir *prm.bar.armLength* sind(j), prm.tag.subcenters(thequa, 2) - prm.bar.armLength*cosd(j)];

            end
        end
    end




    %%  start =============================================================

    vbl = Screen('Flip', prm.w.Number);
    prm.time.expStart  = vbl;
    if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.ExpStart); end

    % Start scanning
    fprintf('\n Starting the task: (%s): %s \n', RunType, datestr(now,'dd-mm-yyyy HH:MM:SS'));


    %% Delay period at the beginning

    fprintf('Waiting for the task to begin in %g seconds...\n', prm.time.countdown);

    if RealRun
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
    end


    %% real task ==========================================================

    if Environment==2 % switch to 1440Hz mode
        Datapixx('Open');
        Datapixx('SetPropixxDlpSequenceProgram', 5); % 2 for 480, 5 for 1440 Hz, 0 for normal
        Datapixx('RegWr');
        fprintf('Switch projector mode to %g hz', prm.monitor.mode_fr(prm.monitor.mode==5))
        Datapixx('SetPropixxLedCurrent', 0, 8);
        Datapixx('SetPropixxLedCurrent', 1, 10);
        Datapixx('SetPropixxLedCurrent', 2, 9);

        Datapixx('RegWrRd');
    end

    vbl = Screen('Flip', prm.w.Number);

    for thetrial = 1: prm.N.trial
        if useEye
            Eyelink('Message', 'Trial %d', thetrial);
        end

        a = find(imgTxts==T.img1Txt(thetrial));
        [thefig, theang, theview] = ind2sub(size(imgTxts), a);
        theRange = squeeze(imgTxts(thefig, :, theview));
        theBarRange = squeeze(lineEnds(:, T.view1(thetrial), :, :));

        currentAng = T.angleProbe(thetrial);
        preRespT = vbl;
        history = currentAng;
        fixFr = randSamp(prm.dur.fixRange, 1, 'n'); % jitter
        theSwitchFr = [fixFr; (prm.SwitchFr +fixFr)];

        nf = 1;  % number of flipped frame (on the monitor). reset to 1 for every trial
        nf_proj = 1;
        fliptimes = nan(1, prm.SwitchFr(end) + max(prm.dur.fixRange));

        holdUp         = 0;
        holdDown       = 0;
        respFlag       = 0;
        timeforDrawing = 0;


        while ~respFlag && nf <= theSwitchFr(end)
            % select the proper colour channel to draw into
            % colorchan will increase by 1 every 4 physical frame and reset after
            % 12 physical frames
            t1 = GetSecs;
            for colorchan = 1:3
                colmask = zeros(4,1);
                colmask(colorchan) = 1;
                Screen('BlendFunction', prm.w.Number, [], [], colmask);

                % quadrant will increase by 1 every physical frame and reset after
                % 4 physical frames
                for quadrant = 1:4

                    if nf <= fixFr % draw dixation

                        Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.tag.subcenters(quadrant,:), 1);
                    elseif  nf <= theSwitchFr(2) % target 1

                        % always draw the tracking sqaure
                        Screen('DrawTexture', prm.w.Number, prm.diode_track.texture, [],...
                            prm.diode_track.positions(quadrant,:), [], [], [], 255 * prm.tag.tag_sigs(prm.diode_track.freq, nf_proj));

                        % draw tag
                        this_tag = prm.tag.tag_sigs(:, nf_proj);
%                         Screen('DrawTextures', prm.w.Number, blackTxts, [0,0,prm.tag.size], prm.tag.pos(quadrant, :, theview),[],[], []); % draw the black line normally
%                         Screen('DrawTextures', prm.w.Number, whiteTxts, [0,0,prm.tag.size], prm.tag.pos(quadrant, :, theview),[],[], this_tag); % use transparency to change
% 

                        % draw stimulus
                        if strcmp(T.stimulus(thetrial), 'Body')
                            Screen('DrawTexture', prm.w.Number, T.img1Txt(thetrial), prm.img.sourceRect, prm.img.pos(quadrant, :));
                        else
                            Screen('DrawLine', prm.w.Number, prm.bar.color, prm.tag.subcenters(quadrant, 1), prm.tag.subcenters(quadrant,2), theBarRange(T.angleTarget(thetrial),quadrant, 1), theBarRange(T.angleTarget(thetrial),quadrant, 2),10);
                        end

                        Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.tag.subcenters(quadrant, :), 1);

                    elseif  nf <= theSwitchFr(3) % mask

                        theimg = max(find((nf- theSwitchFr(2) - prm.Nfr.maskId)>=0));
                        Screen('DrawTexture', prm.w.Number, maskTxts(theimg), [], prm.img.pos(quadrant, :));
                        Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.tag.subcenters(quadrant, :), 1);

                    elseif  nf <= theSwitchFr(4) % delay
                        this_tag = prm.tag.tag_sigs(:, nf_proj);
                        % always draw the tracking sqaure
                        Screen('DrawTexture', prm.w.Number, prm.diode_track.texture, [],...
                            prm.diode_track.positions(quadrant,:), [], [], [], 255 * this_tag(prm.diode_track.freq));


%                         Screen('DrawTextures', prm.w.Number, blackTxts, [0,0,prm.tag.size], prm.tag.pos(quadrant, :, theview),[],[], []); % draw the black line normally
%                         Screen('DrawTextures', prm.w.Number, whiteTxts, [0,0,prm.tag.size], prm.tag.pos(quadrant, :, theview),[],[], this_tag); % use transparency to change
                        %                        Screen('DrawTextures', prm.w.Number, blackTxts(:, theview), [], prm.tag.pos(quadrant, :), [],[],[], ones(4,11));
                        %                         Screen('DrawTexture', prm.w.Number, blackTxts(this_tag,theview));
                        Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color, prm.tag.subcenters(quadrant, :), 1);

                    elseif nf <= theSwitchFr(end) % adjustment
                        if strcmp(T.stimulus(thetrial), 'Body')
                            Screen('DrawTexture', prm.w.Number, theRange(currentAng), prm.img.sourceRect, prm.img.pos(quadrant, :));
                        else
                            Screen('DrawLine', prm.w.Number, prm.bar.color, prm.tag.subcenters(quadrant, 1), prm.tag.subcenters(quadrant,2), theBarRange(currentAng,quadrant,1), theBarRange(currentAng,quadrant,2),10);
                        end
                        Screen('DrawDots', prm.w.Number, [0; 0], prm.fix.size, prm.fix.color2, prm.tag.subcenters(quadrant, :), 1);
                    end

                    nf_proj = nf_proj + 1;
                end
            end
            t2 = GetSecs - t1;
            timeforDrawing = [timeforDrawing; t2];

            % now flip. flip happens every monitor frame, not projector frame
            Screen('DrawingFinished', prm.w.Number);
            vbl = Screen('Flip', prm.w.Number, vbl + .5 * prm.w.ifi);
            fliptimes(nf) = vbl - prm.time.expStart;


            if nf== 1
                if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.FixStart); end
            elseif nf == fixFr + 1
                T.T1Onset(thetrial) = vbl - prm.time.expStart;
                if RealRun
                    if strcmp(T.stimulus(thetrial), 'Body')
                        prm.trigger.btsi.sendTrigger(prm.trigger.TargetStart.body);
                    else prm.trigger.btsi.sendTrigger(prm.trigger.TargetStart.bar);
                    end
                end
            elseif nf == theSwitchFr(3) + 1
                T.delayOnset(thetrial) = vbl - prm.time.expStart;
                if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.DelayStart); end
            elseif nf == theSwitchFr(4) + 1
                T.probeOnset(thetrial) = vbl - prm.time.expStart;
                if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.DelayEnd); end
            end

            nf = nf + 1; % after refreshing, count one more frame

            %% collect response
            if RealRun

                % retreive responses from button box

                if  ~respFlag && nf > theSwitchFr(end-1) % only listen to button box when
                    [keydown,respT] = prm.trigger.btsi.waitResponse(0.0001, true);
                    if keydown == prm.exp.key.up.press %&& respT - preRespT > 1*prm.w.ifi % to control of speed of adjustment
                        holdUp = 1;
                    elseif keydown == prm.exp.key.down.press %&& respT - preRespT > 1*prm.w.ifi
                        holdDown = 1;
                    elseif ismember(keydown, [prm.exp.key.up.release, prm.exp.key.down.release])
                        holdUp = 0; % release the key to stop adjusting
                        holdDown = 0;
                    elseif keydown == prm.exp.key.space
                        respFlag = 1;
                        if RealRun; prm.trigger.btsi.sendTrigger(prm.trigger.Response); end
                        T.resp(thetrial)   =   currentAng;
                        T.rt(thetrial)     =   respT - T.probeOnset(thetrial) - prm.time.expStart;
                        T.adjustHistory{thetrial} = history;
                        fprintf('Responded %g, should be %g.\n', T.resp(thetrial), T.angleTarget(thetrial));
                    end % end of kbcheck
                    preRespT = respT;
                end

                if holdUp ==1 && mod(nf,3) % add one degree every three frames to make it less flexible
                    currentAng = currentAng - 1
                    currentAng = max(currentAng, min(prm.fac.adjustRange));
                    history = [history, currentAng];
                    if length(history)==2
                        T.adjustOnset(thetrial) = respT- prm.time.expStart;
                    end
                elseif holdDown==1 && mod(nf,3)
                    currentAng = currentAng + 1
                    currentAng = min(currentAng, max(prm.fac.adjustRange));
                    history = [history, currentAng];
                    if length(history)==2
                        T.adjustOnset(thetrial) = respT- prm.time.expStart;
                    end
                end
            end

            % also check keyboard, no matter realRun or not.
            [keydown,respT,keyCode] = KbCheck;
            if keydown
                if keyCode(KbName('ESCAPE')) %check quit
                    fprintf('\n\nExperiment terminated at %s, diary closed...\n', datestr(now));
                    save('errorState.mat');
                    ShowCursor;
                    diary off
                    sca; return

                    if useEye %#ok<*UNRCH>
                        % stop
                        Eyelink('Message', 'ENDEXP');
                        Eyelink('StopRecording'); % Stop tracker recording
                        Eyelink('command', 'set_idle_mode');
                        WaitSecs(0.5);
                        Eyelink('CloseFile'); % Close EDF file on Host PC
                        % Transfer the edf file
                        transfer_edf(prm); % try transfering the edf file
                        Eyelink('Shutdown');

                    end

                elseif ~respFlag && nf > theSwitchFr(end-1) %collect response only after the probe appeared


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


        %% break===========================================================


        % a normal break
        if thetrial< prm.N.trial && mod(thetrial, prm.N.trialPerBlock)==0

            if RealRun

                % first thing first, change the projector mode
                Datapixx('SetPropixxDlpSequenceProgram', 0); % 2 for 480, 5 for 1440 Hz, 0 for normal
                Datapixx('RegWrRd');
                fprintf('Switch projector mode to %g hz', prm.monitor.mode_fr(prm.monitor.mode==0))
                WaitSecs(0.2)

                blockPrompt=  sprintf('You finished %d/%d blocks! \n\n Take a rest and wait for the experimenter to continue', T.block(thetrial), prm.N.block);
                instruction_screen(prm, blockPrompt);

                if useEye
                    % recalibrate eye tracker
                    EyelinkEnterSetup(prm.eyelink);WaitSecs(0.50);

                    Eyelink('StartRecording'); WaitSecs(0.050);

                    Eyelink('Message', 'resume');
                end

                % -----------------------------------------------------------------
                % set the projector mode back to 1440Hz
                Datapixx('SetPropixxDlpSequenceProgram', 5); % set the projector back to 1440Hz
                Datapixx('RegWr');
                fprintf('Switch projector mode to %g hz', prm.monitor.mode_fr(prm.monitor.mode==5))

                Datapixx('SetPropixxLedCurrent', 0, 8);
                Datapixx('SetPropixxLedCurrent', 1, 10);
                Datapixx('SetPropixxLedCurrent', 2, 9);

                Datapixx('RegWrRd');
                WaitSecs(0.2)


            end

            WaitSecs(0.2)
            blockPrompt= sprintf('You are going to start another block! \n\n When you are ready, press SPACE to start!');
            instruction_screen(prm, blockPrompt);
            WaitSecs(0.5)
            % make sure vbl is available for the next trial
            vbl = Screen('Flip', prm.w.Number);

        end
        prm.time.fliptimes(thetrial, :) = fliptimes;
    end % end of task


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

    fileName1 = fullfile(RunDir, sprintf('Sub%02d_%s_%g.mat', prm.exp.SubNo, prm.exp.RunType ));
    if ~exist('subInfo', 'var')
        save(fileName1, 'prm','T')
    else
        save(fileName1,'prm','T');
    end
    fileName2 = fullfile(RunDir, sprintf('ResultTable_Sub%02d_%s.csv', prm.exp.SubNo, prm.exp.RunType ));
    writetable(T, fileName2);


    %% EXIT and clean up===================================================
    % PTB performance
    plot_flip_time(prm, 3)

    %% If eyetracking experiment - close the eyetracker
    if useEye
        Eyelink('Message', 'ENDEXP');
        Eyelink('StopRecording'); % Stop tracker recording
        Eyelink('command', 'set_idle_mode');
        WaitSecs(0.5);
        Eyelink('CloseFile'); % Close EDF file on Host PC
        % Transfer the edf file
        transfer_edf(prm); % try transfering the edf file
        Eyelink('Shutdown');
    end

    if Environment==2
        % change the projector mode to normal
        Datapixx('SetPropixxDlpSequenceProgram', 0); % 2 for 480, 5 for 1440 Hz, 0 for normal
        Datapixx('RegWrRd');
        fprintf('Switch projector mode to %g hz', prm.monitor.mode_fr(prm.monitor.mode==0))
    end


    % clean log
    fprintf('Diary closed (%s)\n\n', datestr(now));
    diary off

    % clean memory
    %  cleans up the Java heap from memory leaks, preventing the infamous
    %  Java OutOfMemory exception. https://nl.mathworks.com/matlabcentral/fileexchange/36757-java-heap-cleaner
    try  % use this to avoid a cessation from error
        jheapcl;
    end

    LoadIdentityClut(prm.w.Number)
    sca
    ShowCursor;
    Priority(0);
    toc;

catch ME
    if Environment==2
        % change the projector mode to normal
        Datapixx('SetPropixxDlpSequenceProgram', 0); % 2 for 480, 5 for 1440 Hz, 0 for normal
        Datapixx('RegWrRd');
        fprintf('Switch projector mode to %g hz', prm.monitor.mode_fr(prm.monitor.mode==0))
    end
    save('errorState.mat');
    sca
    ShowCursor;
    Priority(0);
    toc;
    rethrow(ME);
end
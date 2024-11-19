function RIFT_main(prm,id)

% cd(id.path);
%% =========================================================================
% open screen and show instructions
% if isfield(prm.monitor,'gamma')
% prm.monitor    = open_screen(prm.monitor,[],[],prm.monitor.gamma);
% else
% prm.monitor    = open_screen(prm.monitor,[],[]);
% end
prm.monitor      = InitPsychtoolbox(prm.monitor, prm.exp.is_live);

Screen('BlendFunction',prm.monitor.window,'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA');
instruction_screen(prm)
        
%% setup propixx mode
 
WaitSecs(1);

% eye-link calibaration
if prm.exp.eyelink_live
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram', 0); % 2 for 480, 5 for 1440 Hz, 0 for normal
    Datapixx('RegWrRd');
    
    prm = eyelink_init(prm,id); % initialize the eye tracker parameters
    EyelinkEnterSetup(prm.eyelink);% enter the set up for the eyelink
    Screen('FillRect',prm.monitor.window,127); % make it grey again
end

if prm.exp.is_live % switch to 1440Hz mode
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram', 5); % 2 for 480, 5 for 1440 Hz, 0 for normal
    Datapixx('RegWr');

    Datapixx('SetPropixxLedCurrent', 0, 8);
    Datapixx('SetPropixxLedCurrent', 1, 10);
    Datapixx('SetPropixxLedCurrent', 2, 9);

    Datapixx('RegWrRd');
end
%%=========================================================================

% make fixation texture
w_px = prm.monitor.width/2;
h_px = prm.monitor.height/2;
buf    = Screen('OpenOffscreenWindow', prm.monitor.window, [0 0 0 0], [0 0 w_px h_px], [], [], 4);
fixOut = prm.monitor.deg_to_pix*prm.fixation.rad_outer;
fixIn  = prm.monitor.deg_to_pix*prm.fixation.rad_inner;

rect = [w_px/2-fixOut, h_px/2-fixOut, w_px/2+fixOut, h_px/2+fixOut];
Screen('FillOval', buf, 255, rect);
rect = [w_px/2-fixIn, h_px/2-fixIn, w_px/2+fixIn, h_px/2+fixIn];
Screen('FillOval', buf, 0, rect);
prm.fixation.tex = Screen('OpenOffscreenWindow', prm.monitor.window, [0 0 0 0], [0 0 w_px h_px]);
Screen('CopyWindow', buf, prm.fixation.tex);

% diode tracking square
sz = prm.diode_track.SizeInPxl;
prm.diode_track.texture = Screen('OpenOffscreenWindow', prm.monitor.window, [0 0 0 0], [0 0 w_px h_px]);
Screen('FillRect', prm.diode_track.texture, 255, [0 h_px-sz sz h_px]);

%%=========================================================================
% start the eye tracker
if prm.exp.eyelink_live
    Eyelink('StartRecording');
    WaitSecs(0.050);
    Eyelink('Message', 'STARTEXP');
end

%%=========================================================================
% initialize data structure for saving
data           = struct('id',id,'prm',prm);

%%=========================================================================
% main trial loop: 
data.results   = make_table();
block          = 1;
prm.trigger.btsi.sendTrigger(prm.trigger.ExpStart); % trigger for the start of EXP
for i = 1:prm.exp.ntrials
    if prm.exp.eyelink_live
        Eyelink('Message', 'Trial %d', i); 
    end 
    %----------------------------------------------------------------------
    % ITI
    WaitSecs(randsample(.8:.001:1.5,1));    
    %----------------------------------------------------------------------    
    % trial routine
    trial_in         = prm.exp.design(i,:);
    
    % Here the mode is 1440Hz
    % target and non-target change orientation simultaneously
    [trial_out,data] = RIFT_trial(trial_in,data); 
    
%     WaitForResp = true; % whether tester gives a response
%     while WaitForResp == true
%         instruction_screen(prm,'press space to continue');
%         [~, ~, keyCode] = KbCheck;  
%         if keyCode(KbName('ESCAPE'))
%             sca;
%             cd(data.id.path);
%             save([data.id.number '_' prm.exp.name '_exit.mat'],'data');
% %             LoadIdentityClut(prm.monitor.window);
%             sca;ShowCursor;fclose all;
%         elseif keyCode(KbName('SPACE'))
%             WaitForResp = false;
%         end
%     end
    %----------------------------------------------------------------------
    % save log
    data.results     = vertcat(data.results,trial_out);
    %----------------------------------------------------------------------
    % block breaks; until here the projector is still in 1440 mode
    if (i>1 && mod(i,prm.exp.ntrials/prm.exp.blocks)==0)
        %------------------------------------------------------------------
        % first thing first, change the projector mode
        Datapixx('SetPropixxDlpSequenceProgram', 0); % 2 for 480, 5 for 1440 Hz, 0 for normal
        Datapixx('RegWrRd');
        
        %------------------------------------------------------------------
        % End of block memory task
        SbjResp = MemoryTask(data,block);
        data.results.EndOfBlockTask{i} = SbjResp;


        WaitForResp = true; % whether tester gives a response
        while WaitForResp == true && i~=prm.exp.ntrials % don't show if the last trial is finished
  
            instruction_screen(prm,['End of block ' num2str(block) ', we are going to do eye-tracker calibration' ]);
            [~, ~, keyCode] = KbCheck;  
            if keyCode(KbName('Q'))
                sca;
                cd(data.id.path);
                save([data.id.number '_' prm.exp.name '_exit.mat'],'data');

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
                
                LoadIdentityClut(prm.monitor.window);
                sca;ShowCursor;fclose all;
            elseif keyCode(KbName('SPACE'))
                WaitForResp = false;
            end
        end
        
        %% re-calibrate
        if prm.exp.eyelink_live && i~=prm.exp.ntrials % here the projector is still in normal mode
     
            EyelinkEnterSetup(prm.eyelink);
%             Screen('FillRect',prm.monitor.window,127); % make it grey again;
            WaitSecs(0.50);
            
            Eyelink('StartRecording');
            WaitSecs(0.050);
            Eyelink('Message', 'resume');
        
            %% finish re-calibration and show hint for next block's condition

            if prm.exp.design(i,:).fixmode == 1 && i ~= prm.exp.ntrials
                instruction4fix = [];
            else 
                instruction4fix = ', please keep fixating until fixation dot disappears';
            end
            
            DrawFormattedText(prm.monitor.window, ['block ' num2str(block+1) ' will start in 3 seconds' instruction4fix], ...
               'center','center', prm.monitor.white, [], [], [], 1.5);
            Screen('Flip', prm.monitor.window);    

            WaitSecs(3);
        else % if not eyetracker just wait for some time
            WaitSecs(1);
        end
        block  = block+1;

        % ----------------------------------------------------------------- 
        % set the projector mode back to 1440Hz
        Datapixx('SetPropixxDlpSequenceProgram', 5); % set the projector back to 1440Hz
        Datapixx('RegWr');

        Datapixx('SetPropixxLedCurrent', 0, 8);
        Datapixx('SetPropixxLedCurrent', 1, 10);
        Datapixx('SetPropixxLedCurrent', 2, 9);
        
        Datapixx('RegWrRd');
        WaitSecs(1)
    end
   
end

%----------------------------------------------------------------------
prm.trigger.btsi.sendTrigger(prm.trigger.ExpEnd); % trigger for the End of EXP
Datapixx('SetPropixxDlpSequenceProgram', 0); % set the projector to normal mode
Datapixx('RegWrRd');
instruction_screen(prm,'End');

%% If eyetracking experiment - close the eyetracker
if prm.exp.eyelink_live
    Eyelink('Message', 'ENDEXP');
    Eyelink('StopRecording'); % Stop tracker recording
    Eyelink('command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile'); % Close EDF file on Host PC
    % Transfer the edf file
    transfer_edf(prm); % try transfering the edf file
    Eyelink('Shutdown');
end


%%=========================================================================
% save and close
cd(data.id.path);
save([data.id.number '_' prm.exp.name '.mat'],'data');
LoadIdentityClut(prm.monitor.window);
Screen('CloseAll');sca;ShowCursor;fclose all;
function [trial_out,data] = RIFT_trial(trial_in,data)
   
prm                 = data.prm;

%% stimuli prm
cd(prm.stim.path)
% first determine the neutral stimuli
shape         = {'curv','rect'}; 
RWsize        = {'large','small'}; 

% determine the stimuli category
distract_set  = convertCharsToStrings([shape{trial_in.CurvOrRect} '_' RWsize{trial_in.LargeOrSmall} '\']);
oddball_set   = convertCharsToStrings([shape{trial_in.OddCurvOrRect} '_' RWsize{trial_in.OddLargeOrSmall} '\']);

neutral_struct= indir(strcat(prm.stim.path, distract_set));
odd_struct    = indir(strcat(prm.stim.path, oddball_set));

neutral_img   = {neutral_struct.name};
odd_img       = {odd_struct.name};

idistractor   = randsample(3,3); % select 3 out of 4 stimuli categories
iodd          = randsample(3,1); % select 1 out of 4 stimuli categories

stimuli_path  = cell(4,1);

pos_distractor  = Shuffle(setdiff(1:4,trial_in.OddballPos));

stimuli_path{trial_in.OddballPos}  = strcat(prm.stim.path, oddball_set, odd_img{iodd});
for i =1:3
    pos = pos_distractor(i);
    stimuli_path{pos} = strcat(prm.stim.path, distract_set, neutral_img{idistractor(i)});
end

idx_stim = nan(1,4);
Texture_img = zeros(1,4);
Texture_mask = zeros(1,4);
sti         = repmat(zeros(prm.stim.size),[1 1 4]);
for i = 1:numel(stimuli_path)
    cd(stimuli_path{i});
    stim = indir;
    stim_fullpath = strcat(stimuli_path{i}, '\', stim(randsample(numel(stim),1)).name);
    stimuli_path{i}  = stim_fullpath;
    idx_stim(i)      = find(strcmp(prm.stim.AllStimPath, stim_fullpath) == 1);
    [img, ~, alpha]  = imread(stimuli_path{i}); 
    img              = imresize(img,prm.stim.size);
    alpha            = imresize(alpha,prm.stim.size);
    mask             = img;
    mask(alpha~=0)     = 255;
    mask_top(:, :, 1:3)   = repmat(mask,[1 1 3]);
    mask_top(:, :, 4)     = alpha;
    sti(:, :, 1:3)   = repmat(img,[1 1 3]);
    sti(:, :, 4)     = alpha;
    Texture_img(i)   = Screen('MakeTexture', prm.monitor.window, sti);
    Texture_mask(i)  = Screen('MakeTexture', prm.monitor.window, mask_top);
end

% here we define the position of stimuli
Jitters =  prm.stim.JitterInInPxl(1) + rand(1,4)*(diff(prm.stim.JitterInInPxl));
Dis     =  prm.stim.DistanceInPxl(1) + rand(1,4)*(diff(prm.stim.DistanceInPxl));
Angle   = zeros(1,4);
for quadrant = 1:4
    Angle(quadrant)   = prm.stim.AngleInDeg(quadrant,1) + rand*(diff(prm.stim.AngleInDeg(quadrant,:)));
end
% add the additional position coordinate to each center of sub-screen
subscreen_rects = cell(1,4);
for icenter = 1:4
    rects = zeros(4, 4);
    for quadrant = 1:4 % quadrant within in each sub-screen
        theta  = deg2rad(Angle(quadrant));
        R      = Dis(quadrant) + Jitters(quadrant);         
        rects(quadrant,:) = prm.stim.subcenters(icenter,:) + [R * cos(theta) R * sin(theta) R * cos(theta) R * sin(theta)];
    end    
    subscreen_rects{icenter} = rects;
end

% make the background texture
bg_scrambled = phase_scramble(prm.bg.grating);
prm.bg.Texture_bg    = Screen('MakeTexture', prm.monitor.window, bg_scrambled);

%% ========================================================================
% pre-fixation & fixation after stimuli onset
% for i = 1:4
%     Screen('DrawTexture', prm.monitor.window, prm.bg.Texture_bg,[],prm.bg.positions(i ,:));
% end
% Screen('DrawTextures', prm.monitor.window, prm.fixation.tex, [], prm.bg.positions');
% Screen('Flip', prm.monitor.window);
% WaitSecs(prm.fixation.PreDur);
% Screen('Flip', prm.monitor.window, fix_start+ (prm.fixation.duration-.5)*prm.monitor.ifi,0);
 
%% ========================================================================
% tagged presentation loop
PreDur = prm.monitor.frame_ms(Sample(prm.fixation.PreFixRange)) * 12;
Dur    = prm.monitor.frame_ms(Sample(prm.fixation.FixRange)) * 12;
RiftStartSend = 0; % whether the trigger for stimuli onset is send or not
FixEndSend = 0; % whether the trigger for stimuli onset is send or not


phys_frame_ind = 1;
n_restart = 0;
vbl = []; % the flip vbl timestamp
fliptimes = nan(1, size(prm.stim.tag_sigs,2)/12);
flipind = 1;
tStart              = GetSecs;
totol_frame = PreDur + size(prm.stim.tag_sigs,2);

if trial_in.fixmode 
    fixation_frame = PreDur + Dur;
else
    fixation_frame = PreDur;
end

%% main loop
while phys_frame_ind <= totol_frame
    % where are we going to draw it?
%     % quadrant will increase by 1 every 3 physical frames and reset after
%     % 12 physical frames
%     quadrant = mod(floor((phys_frame_ind-1)/3), 4) + 1;
% 
%     % select the proper colour channel to draw into
%     % colorchan will increase by 1 every physical frame and reset after 3
%     % physical frames
%     colorchan = mod(phys_frame_ind-1, 3) + 1;

    % quadrant will increase by 1 every physical frame and reset after
    % 4 physical frames
    quadrant = mod(phys_frame_ind-1, 4) + 1;

    % select the proper colour channel to draw into
    % colorchan will increase by 1 every 4 physical frame and reset after
    % 12 physical frames
    colorchan = mod(floor((phys_frame_ind-1)/4), 3) + 1;

    colmask = zeros(4,1);
    colmask(colorchan) = 1;
    Screen('BlendFunction', prm.monitor.window, [], [], colmask);

    % draw the stimulus with the specified tagging type

    % luminance tagging between 0 and 100%: multiply intensity with tagging signal
    % (black stays black, white goes between 0-100%)

    Screen('DrawTexture', prm.monitor.window, prm.bg.Texture_bg,[],prm.bg.positions(quadrant,:),[],[],0.5);

    
    % fixation after stimuli onset
    if phys_frame_ind <= fixation_frame
        Screen('DrawTextures', prm.monitor.window, prm.fixation.tex, [], prm.bg.positions');
    end
    
    % draw the real stimuli after prefixation
    if phys_frame_ind > PreDur
        Screen('DrawTexture', prm.monitor.window, prm.diode_track.texture, [],...
            prm.bg.positions(quadrant,:), [], [], [], 255 * prm.stim.tag_sigs(prm.diode_track.freq, phys_frame_ind-PreDur));
        for isti = 1:4
            this_tag = prm.stim.tag_sigs(isti, phys_frame_ind-PreDur);
            %color_mod = ones(1, 3) * 255 * this_tag;
            sti_pos = subscreen_rects{quadrant};
            % full color to black
            % Screen('DrawTexture',prm.monitor.window, Texture_img(isti),[],sti_pos(isti,:),[],[],[], color_mod);
            
            % always draw image 
            Screen('DrawTexture',prm.monitor.window, Texture_img(isti),[],sti_pos(isti,:),[],[]);
            
            
            % draw the white mask on top of imgae
            Screen('DrawTexture',prm.monitor.window, Texture_mask(isti),[],sti_pos(isti,:),[],[], this_tag*0.8+0.2);
        end
    end
    % flip if necessary (every 12 physical frames)
    if phys_frame_ind > 1 && mod(phys_frame_ind, 12) == 0
        % draw fixation dot in every quadrant
%         Screen('BlendFunction', prm.monitor.window, [], [], mask_old);
%         Screen('DrawTextures', prm.monitor.window, stim.tex.fix, [], rects');

        if ~isempty(vbl)
            vbl = Screen('Flip', prm.monitor.window, vbl+prm.monitor.ifi/2);
        else % the first frame
            % the first time we're flipping we don't have a previous
            % timestamp yet
            vbl = Screen('Flip', prm.monitor.window);
            prm.trigger.btsi.sendTrigger(prm.trigger.FixStart); % after the first flip, we send the trigger for fixation onset
        end
        
        if phys_frame_ind > PreDur && RiftStartSend == 0
            prm.trigger.btsi.sendTrigger(prm.trigger.RiftStart); % we send the trigger for stimuli onset (where RIFT begin)
            RiftStartSend = 1;
        end
        
        if trial_in.fixmode && phys_frame_ind > fixation_frame && FixEndSend == 0 % only sen this trigger when fixater after stimuli onset
            prm.trigger.btsi.sendTrigger(prm.trigger.FixEnd); % we send the trigger for end of fixation, which is onset of free explore
            FixEndSend = 1;           
        end
        
%         %% use eye-tracker to check the fixation
%         if prm.exp.eyelink_live && phys_frame_ind <= fixation_frame
%             while 1
%                 % THIS Procedure misses really long enduring blinks, but such
%                 % blinks are very very rare (and they are visible from online
%                 % inspection of the data on the eye-tracker screen).
%                 NextDataType = Eyelink('GetNextDataType'); % takes 0.1 ms [sic] to call that function (tested separately).
%                 if NextDataType >= 3 && NextDataType <= 8 % some eye event (el.FIXUPDATE event, which is 9, excluded)
%                     fc_start_time = GetSecs(); % restart timer? Honestly not sure what should be done
%                 elseif NextDataType == 9 % el.FIXUPDATE, usually comes every 50 ms!
%                     FixUpItem = Eyelink('GetFloatData', NextDataType);
%                     % check distance from screen center
%                     if sqrt( (FixUpItem.gavx - prm.monitor.center(1)).^2 + (FixUpItem.gavy - prm.monitor.center(2)).^2 ) > prm.fixation.maxDevPxl
%                         % too far away
%                         Screen('Flip', prm.monitor.window);
%                         Screen('BlendFunction', prm.monitor.window, [], [], [1,1,1,0]);
%                         instruction_screen(prm,'Sorry, this trial will restart. Please fixate! \n Press space to continue');
%                         Screen('Flip', prm.monitor.window);
%                         WaitSecs(1)
%                         n_restart = n_restart+1;
%                         phys_frame_ind = 0;
%                         fliptimes = nan(1, size(prm.stim.tag_sigs,2)/12);
%                         flipind = 1;
%                     end
%                 end
%             end % while 1
%         end
 
        
        %% for debugging when there is no eye-trakcer
        if ~prm.exp.eyelink_live && phys_frame_ind <= fixation_frame
            [~, ~, keyCode] = KbCheck;  
            if keyCode(KbName('b')) % b for break
                Screen('Flip', prm.monitor.window);
                Screen('BlendFunction', prm.monitor.window, [], [], [1,1,1,0]);
                instruction_screen(prm,'Sorry, this trial will restart. Please fixate! \n Press space to continue');
                Screen('Flip', prm.monitor.window);
                WaitSecs(1)
                n_restart = n_restart+1;
                phys_frame_ind = 0;
                fliptimes = nan(1, size(prm.stim.tag_sigs,2)/12);
                flipind = 1;
            end
        end  
        
        fliptimes(flipind) = vbl;
        flipind = flipind + 1;       
    end
    phys_frame_ind = phys_frame_ind + 1;
    %     % for debugging 
    %     doclear = phys_frame_ind > 1 && mod(phys_frame_ind, 12) == 0;
    %     fprintf('flipping, tag=%.4g %%...', tag_sig(phys_frame_ind)*100);
    %     Screen('Flip', prm.monitor.window, 0, double(~doclear));
    %     fprintf('flipped.\n');
    %     
    %     % store image
    %     allimgs{end+1} = Screen('GetImage', prm.monitor.window );
    %     WaitSecs(0.1);
    % %     KbWait();
end

t_end = Screen('Flip',prm.monitor.window);
prm.trigger.btsi.sendTrigger(prm.trigger.RiftEnd);
tduration           = t_end-tStart;
WaitSecs(0.5);

tmp                 = [{stimuli_path}, idx_stim, Jitters, Dis, Angle, PreDur, Dur, tStart, tduration, fliptimes, flipind,n_restart,{[]}];
tmp_tab             = make_table({tmp, 'stimuli_path', 'stim_idx', 'Jitters','Distance','Angle','PreDur','Dur'...
                                'tStart','tduration','fliptimes','flipind','n_restart','EndOfBlockTask'});

trial_out           = horzcat(trial_in,tmp_tab);

Screen('BlendFunction', prm.monitor.window, [], [], [1,1,1,0]);

for i = 1:4
    Screen('Close', Texture_img(i)) 
end
Screen('Close', prm.bg.Texture_bg) 
% Screen('Close', stim.tex.fix) 
% Screen('Close', stim.tex.diode_track) 


end
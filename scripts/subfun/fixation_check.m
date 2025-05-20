function fixation = fixation_check(fixation_each_trial,tri_keep)
% 29.10.2024 Songyun at DCC
% check:
% 1. whether the duration of oddball is longer than distractor
% 2. whether the oddball is fixated more frequent than distractor
% 3. in the condition when sbj keep fixating after stimuli onset, do they
% fixated at oddball more frequent?

    % Predefine trial counts for prefix/postfix shape/size conditions
    N_trial_PreFix = struct('shape', 0, 'size', 0);
    N_trial_PostFix = struct('shape', 0, 'size', 0);
    
    % Predefine first-fixation counts for oddball conditions
    first_fix_oddPre = struct('shape', 0, 'size', 0);
    first_fix_oddPost = struct('shape', 0, 'size', 0);
    
    % Predefine fixation durations
    OddDur = struct('shape', [], 'size', []);
    DistrDur = struct('shape', [], 'size', []);
    AccDur_fix_oddPost = struct('shape', [], 'size', []);
    AccDur_fix_stdPost = struct('shape', [], 'size', []);
    AccDur_fix_oddPre  = struct('shape', [], 'size', []);
    AccDur_fix_stdPre  = struct('shape', [], 'size', []);

    n_AllFix = 0;

    
    for i = 1:numel(tri_keep)
        itrl = tri_keep(i);
        fix_dur = fixation_each_trial.Duration{itrl};
        fix_label = fixation_each_trial.tagged_pos{itrl};
        whether_fix_odd = fixation_each_trial.whether_fix_odd{itrl};
        block_idx = fixation_each_trial.Block_number(itrl);
        oddball_type = fixation_each_trial.OddballShapeOrSize{itrl};  % 1 for shape, 2 for size
    
        % Remove NaN entries
        valid_idx = ~isnan(fix_label);
        fix_label = fix_label(valid_idx);
        whether_fix_odd = whether_fix_odd(valid_idx);
        fix_dur = fix_dur(valid_idx);
    
        % Determine trial type (pre-fix or post-fix) and check first fixation
        if ~isempty(fix_label) % check whethere this is a valid trial
            is_post_fix = mod(block_idx, 2) == 0; % if post_fix
            oddball_key = "shape";
            if oddball_type == 2 % size oddball
                oddball_key = "size";
            end
            
            % Update trial and first-fixation counts
            if is_post_fix % if post fix 
                N_trial_PostFix.(oddball_key) = N_trial_PostFix.(oddball_key) + 1; %coun the trial
                if whether_fix_odd(1) 
                    first_fix_oddPost.(oddball_key) = first_fix_oddPost.(oddball_key) + 1; % count the oddball
                end
                
                % accumulated duration
                AccDur_fix_oddPost.(oddball_key) = [AccDur_fix_oddPost.(oddball_key), sum(fix_dur(whether_fix_odd))];
                AccDur_fix_stdPost.(oddball_key) = [AccDur_fix_stdPost.(oddball_key), sum(fix_dur(~whether_fix_odd))];

            else % pre-fix
                N_trial_PreFix.(oddball_key) = N_trial_PreFix.(oddball_key) + 1;
                if whether_fix_odd(1)
                    first_fix_oddPre.(oddball_key) = first_fix_oddPre.(oddball_key) + 1;
                end

                % accumulated duration
                AccDur_fix_oddPre.(oddball_key) = [AccDur_fix_oddPre.(oddball_key), sum(fix_dur(whether_fix_odd))];
                AccDur_fix_stdPre.(oddball_key) = [AccDur_fix_stdPre.(oddball_key), sum(fix_dur(~whether_fix_odd))];
            end
        end
    
        % Count fixations and durations
        for ifix = 1:numel(fix_label)
            n_AllFix = n_AllFix + 1; % total number of fixation
            if whether_fix_odd(ifix) % whether current fixation is on the oddball 
                OddDur.(oddball_key) = [OddDur.(oddball_key), fix_dur(ifix)]; 
            else
                DistrDur.(oddball_key) = [DistrDur.(oddball_key), fix_dur(ifix)];
            end
        end

    end

    fixation.fixation = n_AllFix;
    fixation.N_trial_PreFix = N_trial_PreFix;
    fixation.N_trial_PostFix = N_trial_PostFix; 

    fixation.first_fix_oddPre = first_fix_oddPre;  % number of trials which sbj made first saccades towards oddball
    fixation.first_fix_oddPost = first_fix_oddPost;


    fixation.OddDur = OddDur;
    fixation.DistrDur = DistrDur;

    fixation.AccDur_fix_oddPost = AccDur_fix_oddPost;
    fixation.AccDur_fix_stdPost = AccDur_fix_stdPost;
    fixation.AccDur_fix_oddPre  = AccDur_fix_oddPre;
    fixation.AccDur_fix_stdPre  = AccDur_fix_stdPre;
end
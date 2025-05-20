function classification(subj_id)
%subj_id = 1
try

    addpath(pwd)
    subjects = Start_up;
    mkdir(fullfile(subjects(subj_id).results, '03_decoding'))

    fprintf('*** SUBJECT %02d : epoching ***\n', subj_id);

    %% load data

    % read trl_keep which contain index of trials kept after pre-proc
    load(fullfile(subjects(subj_id).dir, 'preproc-artifacts-rejectvisual-muscle.mat'),...
        'tri_keep');

    data_MEG = load_clean_data(subj_id);

    cfg = [];
    cfg.trials = find(data_MEG.trialinfo(:,2)>10);
    data_MEG = ft_selectdata(cfg, data_MEG);


    % load behavior data
    load(subjects(subj_id).behav); % variable named as data
    prm.tag.tag_frex

    Behav = T; % this is only for subj 2
    trial_keep_this = tri_keep(11:end) - 10;
    Behav = Behav(trial_keep_this, :);
    assert(height(Behav)==size(data_MEG.trialinfo, 1), 'Trial number does not match behavior!')
    stimuli = {'Body', 'Bar'};

    trial_stim = strcmp(Behav.stimulus,'Body') +1; % 1 for line, 2 for body

    %% classification for body v.s. line as a sanity check
    avg_trialN = 6;

    cfg = [] ;
    cfg.method          = 'mvpa';
    cfg.latency         = [-2 2.4];
    cfg.avgovertime     = 'no';
    cfg.design          = trial_stim; %data_MEG.trialinfo(:,1);
    cfg.features        = 'chan';
    cfg.mvpa            = [];
    cfg.mvpa.classifier = 'multiclass_lda';
    cfg.mvpa.metric     = 'accuracy';
    cfg.mvpa.k          = 3;

    cfg.preprocess      = 'average_samples';
    cfg.preprocess_param = [];
    cfg.preprocess_param.group_size =avg_trialN;

    thestat = ft_timelockstatistics(cfg, data_MEG);


    % plot
    figure; hold on;
    fg = plot(thestat.time, thestat.accuracy');
    fill([-1.5 -1 -1 -1.5], [0.4, 0.4,  1, 1],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
    fill([0, 2.4, 2.4, 0], [0.4, 0.4,  1, 1],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
    yline(1/length(unique(cfg.design)))
    xline(0, '--')
    xlabel('time (s)')
    ylabel('decoding accuracy')
    title(['Stimuli decoding'])
    saveas(fg,fullfile(subjects(subj_id).results, ['/03_decoding/decode_stimuli' '.jpg']));


                %% time genralization of stimuli to check the pinging effect


                cfg = [] ;
                cfg.method          = 'mvpa';
                cfg.generalize      = 'time';
                cfg.latency         = [-2 2.4];
                cfg.avgovertime     = 'no';
                cfg.design          = trial_stim;
                cfg.features        = 'chan';
                cfg.mvpa            = [];
                cfg.mvpa.classifier = 'multiclass_lda';
                cfg.mvpa.metric     = 'accuracy';
                cfg.mvpa.k          = 3;

                thestat = ft_timelockstatistics(cfg, data_MEG);
                stat_temp_generalization= thestat;


                save(fullfile(subjects(subj_id).results, '03_decoding', 'decode_generalization_stimuli.mat'),...
                    'stat_temp_generalization','cfg');

                % plot
                mv_plot_result(stat_temp_generalization.mvpa, stat_temp_generalization.time, stat_temp_generalization.time)


    %% classification of angles
avg_trialN = 6;

    for thestim = 1:2
        for theview = 1:2

            % use only the same facing direction the same stimulus for
            % decoding
            trials = find(trial_stim==thestim& Behav.view1==theview); % 1 for line, 2 for body
    
            cfg = [];
            cfg.trials = trials;
            thedata = ft_selectdata(cfg, data_MEG);

            cfg_angle = [] ;
            cfg_angle.method          = 'mvpa';
            cfg_angle.latency         = [-2 2.4];
            cfg_angle.avgovertime     = 'no';
            cfg_angle.design          = Behav.angleTarget(trials);
            cfg_angle.features        = 'chan';
            cfg_angle.mvpa            = [];
            cfg_angle.mvpa.classifier = 'multiclass_lda';
            cfg_angle.mvpa.metric     = 'accuracy';
            cfg_angle.mvpa.k          = 3;
            cfg_angle.preprocess      = 'average_samples';
            cfg_angle.preprocess_param = [];
            cfg_angle.preprocess_param.group_size =avg_trialN;

            thestat = ft_timelockstatistics(cfg_angle, thedata);
            stat_angle{thestim,theview} = thestat;

            % plot
            figure; hold on;
            fg = plot(thestat.time, thestat.accuracy');
            fill([-1.5 -1 -1 -1.5], [0, 0,  0.4, 0.4],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
            fill([0, 2.4, 2.4, 0], [0, 0,  0.4, 0.4],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
            yline(1/length(unique(cfg_angle.design)))
            xline(0, '--')
            xlabel('time (s)')
            ylabel('decoding accuracy')
            title(['Angle decoding for stimulus: ' stimuli{thestim}])
            saveas(fg,fullfile(subjects(subj_id).results, ['/03_decoding/decode_angle_' stimuli{thestim} num2str(theview) '.jpg']));
        end
    end

    save(fullfile(subjects(subj_id).results, '03_decoding', 'decode_angle_avg_trial.mat'),...
        'stat_angle','cfg_angle');

    %% classification of angles, binary
avg_trialN = 6;

Behav.angleUP = Behav.angleTarget;
Behav.angleUP(Behav.angleTarget==45) =  0;
Behav.angleUP(Behav.angleTarget<45)  =  1;
Behav.angleUP(Behav.angleTarget>45)  =  2;


    for thestim = 1:2
        for theview = 1:2

            % use only the same facing direction the same stimulus for
            % decoding
            trials = find(trial_stim==thestim & Behav.view1==theview & Behav.angleUP~=0); % 1 for line, 2 for body
    
            cfg = [];
            cfg.trials = trials;
            thedata = ft_selectdata(cfg, data_MEG);

            cfg_angle = [] ;
            cfg_angle.method          = 'mvpa';
            cfg_angle.latency         = [-2 2.4];
            cfg_angle.avgovertime     = 'no';
            cfg_angle.design          = Behav.angleUP(trials);
            cfg_angle.features        = 'chan';
            cfg_angle.mvpa            = [];
            cfg_angle.mvpa.classifier = 'multiclass_lda';
            cfg_angle.mvpa.metric     = 'accuracy';
            cfg_angle.mvpa.k          = 3;

            cfg_angle.preprocess      = 'average_samples';
            cfg_angle.preprocess_param = [];
            cfg_angle.preprocess_param.group_size = avg_trialN;

            thestat = ft_timelockstatistics(cfg_angle, thedata);
            stat_angle{thestim,theview} = thestat;

            % plot
            figure; hold on;
            fg = plot(thestat.time, thestat.accuracy');
            fill([-1.5 -1 -1 -1.5], [0, 0,  1, 1],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
            fill([0, 2.4, 2.4, 0], [0, 0,  1, 1],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
            yline(1/length(unique(cfg_angle.design)))
            xline(0, '--')
            xlabel('time (s)')
            ylabel('decoding accuracy')
            title(['Angle decoding for stimulus: ' stimuli{thestim}])
            saveas(fg,fullfile(subjects(subj_id).results, ['/03_decoding/decode_angle_binary_avg' stimuli{thestim} num2str(theview) '.jpg']));
        end
    end

    save(fullfile(subjects(subj_id).results, '03_decoding', 'decode_angle_binary_avg_trial.mat'),...
        'stat_angle','cfg_angle');


    %% time genralization of angle

    for thestim = 1:2
        for theview = 1:2
            trials = find(trial_stim==thestim & Behav.view1==theview); % 1 for line, 2 for body

           cfg_select = [];
            cfg_select.trials = trials;
            thedata = ft_selectdata(cfg_select, data_MEG);
            
            cfg = [] ;
            cfg.method          = 'mvpa';
            cfg.generalize      = 'time';
            cfg.latency         = [-2 2.4];
            cfg.avgovertime     = 'no';
            cfg.design          = Behav.angleTarget(trials);
            cfg.features        = 'chan';
            cfg.mvpa            = [];
            cfg.mvpa.classifier = 'multiclass_lda';
            cfg.mvpa.metric     = 'accuracy';
            cfg.mvpa.k          = 3;

            thestat = ft_timelockstatistics(cfg, thedata);
            stat_temp_generalization{thestim, theview} = thestat;
        end
    end

    save(fullfile(subjects(subj_id).results, '03_decoding', 'decode_generalization.mat'),...
        'stat_temp_generalization','cfg');

    %% plot
    mv_plot_result(stat_temp_generalization{1,1}.mvpa, stat_temp_generalization{1,1}.time, stat_temp_generalization{1}.time)

catch ME
    save(fullfile(subjects(subj_id).results, '03_decoding','errorState.mat'));
    rethrow(ME);
end
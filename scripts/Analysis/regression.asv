function regression(subj_id)
%subj_id = 1
try

Start_up;



    subjects = datainfo();
    mkdir(fullfile(subjects(subj_id).results, '03_decoding'))

    fprintf('*** SUBJECT %02d : epoching ***\n', subj_id);

    %% load data

    % read trl_keep which contain index of trials kept after pre-proc
    load(fullfile(subjects(subj_id).dir, 'preproc-artifacts-rejectvisual-muscle.mat'),...
        'tri_keep');

    data_MEG = load_clean_data(subj_id);

    % load behavior data
    load(subjects(subj_id).behav); % variable named as data
    prm.tag.tag_frex

    Behav = [T;T];
    Behav = Behav(tri_keep, :);
    assert(height(Behav)==size(data_MEG.trialinfo, 1), 'Trial number does not match behavior!')

%% sort the data in the format that MVPA-LIGHT needs

for thetrial = 1:length(data_MEG.trial)
    thesig = data_MEG.trial{thetrial};
    allData(thetrial, :,:)  = thesig;
end

    %% classification of angles
    stimuli = {'Body', 'Line'};

    trial_stim = strcmp(Behav.stimulus,'Body') +1; % 1 for line, 2 for body
angles = unique(Behav.angleTarget);

    for thestim = 1:2
    % use only the same facing direction the same stimulus for
            % decoding
        for theview = 1:2


            % for theang = 1:length(angles)
            % 
            % trials = find(trial_stim==thestim& Behav.view1==theview & Behav.angleTarget==angles(theang)); % 1 for line, 2 for body
            %  occi_chan = find(cellfun(@(x) contains(x, 'O', 'IgnoreCase', false), data_MEG.label));
            % 
            % meanERP = squeeze(mean(mean(allData(trials,occi_chan,:),1),2));
            % plot(data_MEG.time{1}(1:800), meanERP(1:800),'Color', theang*[0.1,0.1,0.11])
            % hold on 
            % 
            % end

            trials = find(trial_stim==thestim& Behav.view1==theview); % 1 for line, 2 for body
            thedata = allData(trials,:,:);
            labels = Behav.angleTarget(trials);

            
            cfg = [];
            cfg.metric      = 'mae';

            thestat = mv_regress(cfg, thedata, labels);

            close all
            plot(data_MEG.time{1}, thestat)
            xlabel('Time'), ylabel('R sqaure')

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

    save(fullfile(subjects(subj_id).results, '03_decoding', 'decode_angle.mat'),...
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
            cfg.generalize       = 'time';
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
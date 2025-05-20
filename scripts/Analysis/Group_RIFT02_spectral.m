close all; clear all
addpath('/project/3018085.01/scripts/subfun')
rootdir = '/project/3018085.01';
subjects = Start_up(rootdir);

plotchans = {'MLO11', 'MLO12', 'MLO13', 'MLO14', 'MLO21', 'MLO22', 'MLO23',...
    'MLO24', 'MLO31', 'MLO32', 'MLO34', 'MLO41', 'MLO42', 'MLO43', 'MLO44',...
    'MLO53', 'MLP31', 'MLP41', 'MLP42', 'MLP51', 'MLP52', 'MLP53', 'MLP54',...
    'MLP55', 'MLT16', 'MLT27', 'MLT47', 'MLT57', 'MRO11', 'MRO12', 'MRO13',...
    'MRO14', 'MRO21', 'MRO22', 'MRO23', 'MRO24', 'MRO31', 'MRO32', 'MRO34',...
    'MRO41', 'MRO42', 'MRO43', 'MRO44', 'MRO53', 'MRP31', 'MRP41', 'MRP42',...
    'MRP51', 'MRP52', 'MRP53', 'MRP55', 'MRT16', 'MRT27', 'MRT47', 'MRT57',...
    'MZO01', 'MZO02', 'MZP01'};

subs = [1 2 3 6];
Nsub = length(subs);
for thesub = 1:Nsub
    
    subj_id = subs(thesub);
    
    fprintf('*** SUBJECT %02d : loading ***\n', subj_id);
    
    % read trl_keep which contain index of trials kept after pre-proc
    load(fullfile(subjects(subj_id).dir, 'preproc-artifacts-rejectvisual-muscle.mat'),...
      'tri_keep');
         
    
    %    load behavorial data
    load(subjects(subj_id).behav); % variable named as data

    allBehav{thesub} = T(tri_keep, :);
    
    allFreqs(thesub, :) = prm.tag.tag_frex; % order: first one is the upmost
    allFreqs


    % exclude unresponded trial from both behavior and MEG
    validtrial = ~isnan(allBehav{thesub}.rt);
    allBehav{thesub} = allBehav{thesub}(validtrial,:);

    load(fullfile(subjects(subj_id).results, '02_spectral', 'freq_trials.mat'));

    cfg = [];
    cfg.trials   = validtrial;
    freq_valid = ft_selectdata(cfg, freq_trials);

    allFourier{thesub} = freq_valid;
    assert(height(allBehav{thesub})==size(allFourier{thesub}.trialinfo,1), 'Trial numbers don''t match!')

end

angles = sort(unique(allBehav{1}.angleTarget));
angles
stimuli = {'Body', 'Bar'};

chanids = match_str(allFourier{thesub}.label, plotchans);
n_added_chan = 4;

%% Behavior averaging ================================================

for thesub = 1:Nsub
for thestim = 1:2
    for theang = 1:length(angles)

        trialid = strcmp(allBehav{thesub} .stimulus, stimuli{thestim}) & allBehav{thesub}.angleTarget==angles(theang) ;
        meanError(theang,thestim,thesub) = mean(allBehav{thesub}.Error(trialid)); % only calculate acc for RS trials
        meanResp(theang,thestim,thesub) = mean(allBehav{thesub}.resp(trialid)); % only calculate acc for RS trials
        STD(theang,thestim,thesub) = std(allBehav{thesub}.Error(trialid));
    end
end
end

%% coherence compute and plot


for thesub = 1:Nsub
    figure
    for i = 1:length(stimuli)

        thesti = stimuli{i};
        cfg = [];
        cfg.trials   = find(strcmp(allBehav{thesub}.stimulus,thesti) );
        freq_sti  = ft_selectdata(cfg, allFourier{thesub});


        cfg = [];
        cfg.method = 'coh';
        nchan = numel(allFourier{thesub}.label)-n_added_chan;
        cfg.channelcmb = [allFourier{thesub}.label(1:end-n_added_chan) repmat({ 'tag_all'}, nchan, 1)]; % compute coherence between all chan with 'tagall'
        coh_stim = ft_connectivityanalysis(cfg, freq_sti);

        subplot(1, length(stimuli), i);hold on
        %occi_chan = find(cellfun(@(x) contains(x, 'O', 'IgnoreCase', false), coh_trial.labelcmb(:,1)));
        plot(coh_stim.freq, mean(coh_stim.cohspctrm(chanids,:),1))
        xlabel('frequencies')
        ylabel('coherence')
        title(['mean coherence ' thesti])
        ylim([0,0.8])

         for k = 1: size(allFreqs,2)       % Freqs is already in the order of upper-> lower
                coh_sti(k, i, thesub) = mean(coh_stim.cohspctrm(chanids, coh_stim.freq==allFreqs(thesub, k)));
            end

    end

end

%% now separate by angles

for thesub = 1:Nsub
    for i = 1:length(stimuli)
        figure;
        for j = 1:length(angles)

            thesti = stimuli{i};
            cfg = [];
            cfg.trials   = find(strcmp(allBehav{thesub}.stimulus,thesti) & allBehav{thesub}.angleTarget==angles(j));
            freq_sti  = ft_selectdata(cfg, allFourier{thesub});


            cfg = [];
            cfg.method = 'coh';
            nchan = numel(allFourier{thesub}.label)-n_added_chan;
            cfg.channelcmb = [allFourier{thesub}.label(1:end-n_added_chan) repmat({ 'tag_all'}, nchan, 1)]; % compute coherence between all chan with 'tagall'
            coh_stim = ft_connectivityanalysis(cfg, freq_sti);

            subplot(1, length(angles), j);hold on
            temp = mean(coh_stim.cohspctrm(chanids,:));
            plot(coh_stim.freq(55:70),temp(55:70))
            xlabel('frequencies')
            ylabel('coherence')
            yline(1)
            title([thesti '-' num2str(angles(j))])
            %    saveas(cfg,fullfile(subjects(subj_id).results, '02_spectral', 'cohByAng.jpg'));

            for k = 1: size(allFreqs,2)       % Freqs is already in the order of upper-> lower
                coh_ang(k, j, i, thesub) = mean(coh_stim.cohspctrm(chanids, coh_stim.freq==allFreqs(thesub, k)));
            end

        end
    end
end

save(fullfile(subjects(subj_id).results, '02_spectral', 'coh_tagged'),...
      'coh_ang');

%% split by behavior response
nbin = 5;

for thesub = 1:Nsub
    for i = 1:length(stimuli)

        thesti = stimuli{i};
        stimtrials = strcmp(allBehav{thesub}.stimulus,thesti);

        edges(i,:, thesub) = quantile(allBehav{thesub}.resp(stimtrials), 0:1/nbin:1);
        quantileGroup = discretize(allBehav{thesub}.resp, edges(i,:, thesub) , 'IncludedEdge','right');

        figure;
        for q = 1:nbin

            cfg = [];
            cfg.trials   = find(quantileGroup == q & stimtrials);
            size(cfg.trials,1)
            freq_q  = ft_selectdata(cfg, allFourier{thesub});


            cfg = [];
            cfg.method = 'coh';
            nchan = numel(allFourier{thesub}.label)-n_added_chan;
            cfg.channelcmb = [allFourier{thesub}.label(1:end-n_added_chan) repmat({ 'tag_all'}, nchan, 1)]; % compute coherence between all chan with 'tagall'
            coh_resp = ft_connectivityanalysis(cfg, freq_q);
            % 
            % subplot(1, length(edges)-1, q);hold on
            % temp = mean(coh_resp.cohspctrm(chanids,:));
            % plot(coh_resp.freq(55:70),temp(55:70))
            % xlabel('frequencies')
            % ylabel('coherence')
            % yline(1)
            % title([thesti '-' num2str(edges(q))])

            for k = 1: size(allFreqs,2)       % Freqs is already in the order of upper-> lower
                coh_res(k, q, i, thesub) = mean(coh_resp.cohspctrm(chanids, coh_resp.freq==allFreqs(thesub, k)));
            end

        end
    end
end

%% plot group level ratio

% draw which one? by angle, by resp, by stimuli
metric2draw = coh_ang(:,:,:,:);
x = angles;
x = 1:nbin;

Down2up_Body = metric2draw(2,:,1,:)./metric2draw(1,:,1,:);
Down2up_Line = metric2draw(2,:,2,:)./metric2draw(1,:,2,:);


figure;
for thesub = 1:size(metric2draw, 4)
    subplot(1,2,1)
    plot(x,metric2draw(1,:,1,thesub),'b-','LineWidth',2); hold on;
    plot(x,metric2draw(2,:,1,thesub),'r-','LineWidth',2);
    legend({'upper','lower'})
    title(stimuli{1})


    subplot(1,2,2)
    plot(x,metric2draw(1,:,2,thesub),'b-','LineWidth',2); hold on;
    plot(x,metric2draw(2,:,2,thesub),'r-','LineWidth',2); legend({'upper','lower'})
    title(stimuli{2})
end


for thesub = 1:size(metric2draw, 4)
    figure;
    plot(x, Down2up_Body(:,:,:,thesub),'b-','LineWidth',2)
    hold on;
    plot(x, Down2up_Line(:,:,:,thesub),'r-','LineWidth',2)
end

figure;
for thesub = 1:Nsub
    plot(x, Down2up_Body(:,:,:,thesub),'b-','LineWidth',2)
    hold on;
    plot(x, Down2up_Line(:,:,:,thesub),'r-','LineWidth',2)
end
legend(stimuli)
yline(1)
title('ratiot lower/upper')


%% formal plot============================================================

MarkerSize = 5;
Col.pink = [255 153 153] /255;
Col.green = [91 153 153] /255;
Col.black = [0 0 0];
Col.grey = [.4 .4 .4];
jitter = 0.1;

figure;
groupB = plot(x+jitter,mean(Down2up_Body,4),'o-','MarkerFaceColor',Col.pink,'MarkerEdgeColor',Col.black,'Color',Col.pink); hold on
groupL = plot(x,mean(Down2up_Line,4),'o-','MarkerFaceColor',Col.green,'MarkerEdgeColor',Col.black,'Color', Col.green);
lg = legend([groupB, groupL], stimuli);
set(lg, 'AutoUpdate', 'off');

for ii = 1:numel(x)
    todraw = squeeze(Down2up_Body(:,ii,1,:));
    [Group_mu(ii),Group_sigma(ii), Group_MUCI(ii,:), Group_SigmaCI(ii,:)] = normfit(todraw,0.05);
    cilength(ii) = Group_MUCI(ii,2) - Group_mu(ii);

    scathd = plot(x(ii)+jitter,todraw,'o','MarkerFaceColor',Col.pink,'MarkerEdgeColor',Col.grey,'MarkerSize',MarkerSize-1);
    errorbar(x(ii)+jitter, Group_mu(ii), cilength(ii),  'k','linestyle','none','linewidth',1.2,'alignvertexcenters','on');

    todraw = squeeze(Down2up_Line(:,ii,1,:));
    [Group_mu(ii),Group_sigma(ii), Group_MUCI(ii,:), Group_SigmaCI(ii,:)] = normfit(todraw,0.05);
    cilength(ii) = Group_MUCI(ii,2) - Group_mu(ii);

    scathd = plot(x(ii),todraw,'o','MarkerFaceColor',Col.green,'MarkerEdgeColor',Col.grey,'MarkerSize',MarkerSize-1);
    errorbar(x(ii), Group_mu(ii), cilength(ii),  'k','linestyle','none','linewidth',1.2,'alignvertexcenters','on');

end
yline(1,':')
set(gca,'XTick',x, 'XTickLabel', x,'fontsize',12);

ylabel('Lower <- Coherence ratio -> Higher','fontsize',12);
title('MEG')
box off
axis ij
xlim([0.8,5.2])


%% plot behavior
plotError = 0;
x = angles;
if plotError
    Y = meanError;
else
    Y = meanResp;
end
figure;
groupB = plot(x,mean(Y(:,1,:),3),'o-','MarkerFaceColor',Col.pink,'MarkerEdgeColor',Col.black,'Color',Col.pink); hold on
groupL = plot(x+1,mean(Y(:,2,:),3),'o-','MarkerFaceColor',Col.green,'MarkerEdgeColor',Col.black,'Color', Col.green)
lg = legend([groupB, groupL], stimuli,'Location','northeast')
set(lg, 'AutoUpdate', 'off');

for ii = 1:numel(x)
    todraw = squeeze(Y(ii,1,:));
    [Group_mu(ii),Group_sigma(ii), Group_MUCI(ii,:), Group_SigmaCI(ii,:)] = normfit(todraw,0.05);
    cilength(ii) = Group_MUCI(ii,2) - Group_mu(ii);

    scathd = plot(x(ii),todraw,'o','MarkerFaceColor',Col.pink,'MarkerEdgeColor',Col.grey,'MarkerSize',MarkerSize-1);
    errorbar(x(ii), Group_mu(ii), cilength(ii),  'k','linestyle','none','linewidth',1.2,'alignvertexcenters','on');

    todraw = squeeze(Y(ii, 2,:));
    [Group_mu(ii),Group_sigma(ii), Group_MUCI(ii,:), Group_SigmaCI(ii,:)] = normfit(todraw,0.05);
    cilength(ii) = Group_MUCI(ii,2) - Group_mu(ii);

    scathd = plot(x(ii)+1,todraw,'o','MarkerFaceColor',Col.green,'MarkerEdgeColor',Col.grey,'MarkerSize',MarkerSize-1);
    errorbar(x(ii)+1, Group_mu(ii), cilength(ii),  'k','linestyle','none','linewidth',1.2,'alignvertexcenters','on');

end
yline(0,':')
set(gca,'XTick',x, 'XTickLabel', x,'fontsize',12);
if plotError
    ylabel('Lower <- Error (degree) -> Higher','fontsize',12);
else
    ylabel('Lower <- Response (degree) -> Higher','fontsize',12);
    ylim([25, 65])
end
title('Behavior')
box off
axis equal
axis ij

xlim([29, 62])

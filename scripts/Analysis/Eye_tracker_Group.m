  clear  
rootdir = '/project/3018085.01';
addpath(rootdir)
addpath('/project/3018085.01/scripts/subfun')

subjects = Start_up(rootdir);

subs = [1 2 3 6];
Nsub = numel(subs);
for thesub = 1:Nsub
    subj_id = subs(thesub);

    fprintf('*** SUBJECT %02d : loading ***\n', subj_id);

    % read trl_keep which contain index of trials kept after pre-proc
    load(fullfile(subjects(subj_id).results, '01_eye_tracker', 'alltrials.mat') );
    assert(height(Behav)== size(data_eye_aligned,1),'Trial numbers don''t match!')


    allBehav{thesub} = Behav;
    alleye{thesub} = data_eye_aligned;


    stimuli = {'Body', 'Line'};

end


%% plot grand average of x y location 

data_eye_aligned = [];


for thestim = 1:2
    for thesub = 1:Nsub
        trial_stim = strcmp(allBehav{thesub}.stimulus,'Bar') +1; % 1 for line, 2 for body

        % use only the same facing direction the same stimulus
        trials = find(trial_stim==thestim); % 1 for line, 2 for body

        thedata= alleye{thesub}(trials, :,:);

        coordi = squeeze(mean(thedata,1, "omitmissing"));
        deg = atand(coordi(:,1)./coordi(:,2));
        avg_eye(thestim,:,:,thesub) = [coordi, deg];
    end
end

grand_ave = mean(avg_eye, 4);

% another way to compute angle is to use the grand average
for thestim = 1:2
    grand_ave(thestim,:,5) = atand(grand_ave(thestim,:,1)./grand_ave(thestim,:,2));
end

%% plot

MarkerSize = 5;
Col.pink = [255 153 153] /255;
Col.green = [91 153 153] /255;
Col.black = [0 0 0];
Col.grey = [.4 .4 .4];


% Each participant
figure

dim2draw = [1 2 5];
labels = {'horizontal', 'vertical','degree'};
for i = [3]
    for thesub = 1:Nsub
        thedim = dim2draw(i);
        subplot(1,Nsub,thesub)
        plot(t_points, avg_eye(1,:,thedim,thesub), 'Color',Col.pink); hold on;
        plot(t_points, avg_eye(2,:,thedim,thesub),'Color',Col.green); hold on;

        fill([-1.5 -1 -1 -1.5], [-30,-30, 30, 30],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
        fill([0, 2.4, 2.4, 0],  [-30,-30, 30, 30],'blue','LineStyle', 'none', 'FaceAlpha',0.1)

        legend(stimuli)
        title([labels{i} '-sub ' num2str(thesub)])
        xline(0); yline(0);
        set(gca, 'YDir', 'reverse') % to make the position consistent with positions on the screen
    end
end


% Group average
figure
dim2draw = [1 2 5];
labels = {'horizontal', 'vertical','degree'};
for i = 1:length(dim2draw)
    thedim = dim2draw(i);
    subplot(1,length(dim2draw),i)
    plot(t_points, grand_ave(1,:,thedim),'Color',Col.pink); hold on;
    plot(t_points, grand_ave(2,:,thedim),'Color',Col.green); hold on;

    %legend(stimuli)
    title([labels{i}])
    xline(0);
    if thedim==5
        yline(45);
        fill([-1.5 -1 -1 -1.5], 100*[-1,-1, 1, 1],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
        fill([0, 2.4, 2.4, 0], 100*[-1,-1, 1, 1],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
        ylabel('Angle (degree)')
        
    else
        yline(0);
        fill([-1.5 -1 -1 -1.5], [-30,-30, 30, 30],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
        fill([0, 2.4, 2.4, 0],  [-30,-30, 30, 30],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
        ylabel('Position (pixel)')
    end
    box off
    set(gca, 'YDir', 'reverse') % to make the position consistent with positions on the screen
end


%% plot grand average by target angle

fig_pos = figure;
fig_psize = figure;
fig_ang = figure;
colors = hot(7);

for thesub = 1:Nsub

    angles = sort(unique(allBehav{thesub}.angleTarget));
    angles
    figure(); % Activate first figure

    for thestim = 1:2
        trial_stim = strcmp(allBehav{thesub}.stimulus,'Bar') +1; % 1 for line, 2 for body

        for theang = 1:length(angles)

            % use only the same facing direction the same stimulus for
            % decoding
            trials = find(trial_stim==thestim & allBehav{thesub}.angleTarget==angles(theang)); % 1 for line, 2 for body

            thedata= alleye{thesub}(trials, :,:);

            avg = squeeze(mean(thedata,1, "omitmissing"));
            deg = atand(avg(:,1)./avg(:,2));

            avg_eye_aligned{thesub}(:,:,thestim,theang)  = [avg,deg];

            % plot

            subplot(1,2,thestim);hold on;
            %plot(t_points, avg_eye_aligned{thestim}(:,1),'Color', colors(theang,:));
            ylim([-40, 40])
            plot(t_points, squeeze(avg_eye_aligned{thesub}(:,2,thestim,theang)),'Color', colors(theang,:));

            fill([-1.5 -1 -1 -1.5], [-100,-100, 100, 100],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
            fill([0, 2.4, 2.4, 0],  [-100,-100, 100, 100],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
            axis ij
            yline(0)
            %legend({'horizontal','vertical'})
            title(['eye position: ' stimuli{thestim} ])

            % figure(fig_ang)
            % subplot(2,1,thestim)
            % plot(t_points, avg_eye_aligned{thestim}(:,5)); hold on;
            % yline(45)
            % fill([-1.5 -1 -1 -1.5], [-100, -100, 100, 100],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
            % fill([0, 2.4, 2.4, 0],  [-100, -100, 100, 100],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
            %
            % title(['eye position: ' stimuli{thestim} ])
            %
            % % plot pupil size
            % figure(fig_psize);
            % plot(t_points, avg_eye{thestim}(:,4)); hold on;
            % fill([-1.5 -1 -1 -1.5], [-2000,-2000, 500, 500],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
            % fill([0, 2.4, 2.4, 0],  [-2000,-2000, 500, 500],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
            % title(['pupil size: ' stimuli{thestim} ])

        end
    end
end


%% average
figure

for thestim = 1:2
    trial_stim = strcmp(allBehav{thesub}.stimulus,'Bar') +1; % 1 for line, 2 for body

    for theang = 1:length(angles)
        thedata = cellfun(@(x) x(:,:,thestim,theang), avg_eye_aligned, 'UniformOutput', false);
        group = reshape(cell2mat(thedata), [], 5, Nsub);
        group_mean(:,:,thestim,theang) = mean(group, 3);

        group_mean(:,5,thestim, theang) = atand(group_mean(:,1,thestim, theang)./group_mean(:,2, thestim, theang));

        subplot(1,2,thestim);hold on;
        %plot(t_points, avg_eye_aligned{thestim}(:,1),'Color', colors(theang,:));
        ylim([-40, 80])
        plot(t_points, squeeze(group_mean(:,5,thestim,theang)),'Color', theang/6*[1 1 1]);

    end
    fill([-1.5 -1 -1 -1.5], [-100,-100, 100, 100],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
    fill([0, 2.4, 2.4, 0],  [-100,-100, 100, 100],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
    axis ij
    yline(0)
    %legend({'horizontal','vertical'})
    title(['eye position: ' stimuli{thestim} ])

end

%% average over time window

winSize = 0.4;
timewindows = [-1.5:winSize:3-winSize; (-1.5+winSize):winSize:3];
nWin = size(timewindows, 2);
data_avg_time = [];
for i = 1:size(timewindows, 2)
    t = t_points>timewindows(1,i) & t_points<timewindows(2,i);

    thedata = mean(group_mean(t,:,:,:), 1, "omitmissing");
    data_avg_time = cat(1,data_avg_time, thedata);

end   

figure
for thestim = 1:2

    for theang = 1:length(angles)

        subplot(1,2,thestim);hold on;
        %plot(t_points, avg_eye_aligned{thestim}(:,1),'Color', colors(theang,:));
        ylim([-40, 80])
        plot(timewindows(1,:), squeeze(data_avg_time(:,5,thestim,theang)),'Color', theang/6*[1 1 1]);

    end
    fill([-1.5 -1 -1 -1.5], [-100,-100, 100, 100],'yellow','LineStyle', 'none', 'FaceAlpha',0.1)
    fill([0, 2.4, 2.4, 0],  [-100,-100, 100, 100],'blue','LineStyle', 'none', 'FaceAlpha',0.1)
    axis ij
    yline(0)
    %legend({'horizontal','vertical'})
    title(['eye position: ' stimuli{thestim} ])

end

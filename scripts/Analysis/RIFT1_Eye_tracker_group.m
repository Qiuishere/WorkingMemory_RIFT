% 29.10.2024 Songyun at DCC
% check:
% 1. whether the duration of oddball is longer than distractor
% 2. whether the oddball is fixated more frequent than distractor
% 3. in the condition when sbj keep fixating after stimuli onset, do they
% fixated at oddball more frequent?
addpath('/project/3018083.01/scripts/subfun')
addpath('/home/predatt/sonbai/Documents/MATLAB/fieldtrp/')
addpath(genpath('/home/predatt/sonbai/toolbox/BEIM_toolbox-main/matlabv/'))

clear all; clc; close all
ft_defaults
ft_warning off

Results_fold      = '/project/3018083.01/results/';
file_path         = '/01_eye_tracker/fixation_each_trial.mat';
save_path         = '/project/3018083.01/group results/01_beh';
cd(Results_fold)
ss              = indir;
N_sbj           = numel(ss);

subjects = datainfo();
sbj2exclude = [11 12];

colors          = linspecer(2);

%% read individual results and run group analysis

% initialize the matrix 
Dur_SizeOdd    = nan(N_sbj, 1);
Dur_ShapeOdd   = nan(N_sbj, 1);
Dur_SizeDistr  = nan(N_sbj, 1);
Dur_ShapeDistr = nan(N_sbj, 1);

% how frequent were fixations on the oddball
P_odd_size     = nan(N_sbj, 1); 
P_odd_shape    = nan(N_sbj, 1); 

% check whether first fix is at the oddball
P_FirstOdd_pre_size      = nan(N_sbj, 1); 
P_FirstOdd_post_size     = nan(N_sbj, 1); 
P_FirstOdd_pre_shape     = nan(N_sbj, 1); 
P_FirstOdd_post_shape    = nan(N_sbj, 1); 

% accumulate duration
Acc_Dur_SizeOdd_post    = nan(N_sbj, 1); 
Acc_Dur_ShapeOdd_post   = nan(N_sbj, 1); 
Acc_Dur_SizeStd_post    = nan(N_sbj, 1); 
Acc_Dur_ShapeStd_post    = nan(N_sbj, 1); 

Acc_Dur_SizeOdd_pre     = nan(N_sbj, 1); 
Acc_Dur_ShapeOdd_pre    = nan(N_sbj, 1); 
Acc_Dur_SizeStd_pre     = nan(N_sbj, 1); 
Acc_Dur_ShapeStd_pre     = nan(N_sbj, 1); 

for i = 1:N_sbj
    % load the eye tracking data
    if exist(fullfile(Results_fold, ss(i).name ,file_path)) && ~ismember(i, sbj2exclude) 
        load(fullfile(Results_fold, ss(i).name,file_path))
        % load which trial is kept after the pre-processing
        load(fullfile(subjects(i).dir, 'preproc-artifacts-rejectvisual-muscle.mat'),...
            'tri_keep');

        fix = fixation_check(fixation_each_trial,tri_keep);
    
        Dur_SizeOdd(i)    = mean(fix.OddDur.size);
        Dur_ShapeOdd(i)   = mean(fix.OddDur.shape);
        Dur_SizeDistr(i)  = mean(fix.DistrDur.size);
        Dur_ShapeDistr(i) = mean(fix.DistrDur.shape);
    
        P_odd_size(i)     = numel(fix.OddDur.size)/(numel(fix.DistrDur.size) + numel(fix.OddDur.size));
        P_odd_shape(i)    = numel(fix.OddDur.shape)/(numel(fix.DistrDur.shape) + numel(fix.OddDur.shape));
     
        P_FirstOdd_pre_size(i)    = fix.first_fix_oddPre.size / fix.N_trial_PreFix.size;
        P_FirstOdd_pre_shape(i)   = fix.first_fix_oddPre.shape / fix.N_trial_PreFix.shape;
        P_FirstOdd_post_size(i)   = fix.first_fix_oddPost.size / fix.N_trial_PostFix.size;
        P_FirstOdd_post_shape(i)  = fix.first_fix_oddPost.shape / fix.N_trial_PostFix.shape;

        Acc_Dur_SizeOdd_post(i) = mean(fix.AccDur_fix_oddPost.size);
        Acc_Dur_ShapeOdd_post(i) = mean(fix.AccDur_fix_oddPost.shape);
        Acc_Dur_SizeStd_post(i) = mean(fix.AccDur_fix_stdPost.size)/3;
        Acc_Dur_ShapeStd_post(i) = mean(fix.AccDur_fix_stdPost.shape)/3;

        Acc_Dur_SizeOdd_pre(i) = mean(fix.AccDur_fix_oddPre.size);
        Acc_Dur_ShapeOdd_pre(i) = mean(fix.AccDur_fix_oddPre.shape);
        Acc_Dur_SizeStd_pre(i) = mean(fix.AccDur_fix_stdPre.size)/3;
        Acc_Dur_ShapeStd_pre(i) = mean(fix.AccDur_fix_stdPre.shape)/3;

    end 
end





%% violin plot to see the duration difference of size oddball
fg1 = figure('Position',[705 208 589 688]);
points     = [Dur_SizeOdd Dur_SizeDistr];
plot_violin(points,'x',1:2,'showdots',1,'color',colors,'rescaling',.25,'connline',[.5 .5 .5],'dotssize',60,'jitfactor',.05,'dotshift',-.02,'allline',10)
xlim([0.5 2.5]);ylim([200 450]);%yticks([0.45:.1:.95])
format_figure(nan,nan,'','Duration(ms)',[],[],20)
set(gca,'xtick',1:2,'xticklabel',{'Oddball','Standard'})
box on
tl = title('Size'); tl.FontWeight = 'normal';

saveas(fg1,fullfile(save_path, 'fix_duration_size.jpg'));

ttest(Dur_SizeOdd, Dur_SizeDistr)
%% violin plot to see the duration difference of shape oddball
fg2 = figure('Position',[705 208 589 688]);
points     = [Dur_ShapeOdd Dur_ShapeDistr];
plot_violin(points,'x',1:2,'showdots',1,'color',colors,'rescaling',.25,'connline',[.5 .5 .5],'dotssize',60,'jitfactor',.05,'dotshift',-.02,'allline',10)
xlim([0.5 2.5]);ylim([200 450]);%yticks([0.45:.1:.95])
format_figure(nan,nan,'','Duration(ms)',[],[],20)
set(gca,'xtick',1:2,'xticklabel',{'Oddball','Standard'})
box on
tl = title('Shape'); tl.FontWeight = 'normal';

saveas(fg2,fullfile(save_path, 'fix_duration_shape.jpg'));

ttest(Dur_ShapeOdd, Dur_ShapeDistr)

%% bar plot of p_odd to compare with 25%
plt                 = [];
plt.x               = 1:2;
plt.color           = colors([2 1],:);
plt.err_width       = 2.5;
plt.err_color       = [.5 .5 .5];
plt.cap_size        = 0;
plt.dotsout         = .2;
plt.sizedot         = 20;
% plt.allline         = true;

fg3 = figure('position',[1024 369 689 497]);
freq_oddball            = [P_odd_size P_odd_shape];
pl1                 = plot_bar(freq_oddball,plt);

ylim([0 0.3])
format_figure(0.25,nan,'oddball type','probability of fixation');
% format_legend([pl1(2).h(1) pl1(2).h(2)],{'size','shape'});
set(gca,'xtick',[1 2],'xticklabel',{'size','shape'})


saveas(fg3,fullfile(save_path, 'oddball_prob.jpg'));

ttest(P_odd_size, P_odd_shape)


%% check whether fixation condition affect first fixation type
freq_FirstFix            = [P_FirstOdd_pre_size P_FirstOdd_post_size];

fg4 = figure('position',[1024 369 689 497]);
pl1                 = plot_bar(freq_FirstFix,plt);

ylim([0 0.35])
format_figure(0.25,nan,'fixation type','probability of first fixation');
% format_legend([pl1(2).h(1) pl1(2).h(2)],{'size','shape'});
set(gca,'xtick',[1 2],'xticklabel',{'pre','post'})
tl = title('Size'); tl.FontWeight = 'normal';


saveas(fg4,fullfile(save_path, 'first_fix_size.jpg'));

[h,p,~,~] = ttest(P_FirstOdd_pre_size - 0.25, 0)
ttest(P_FirstOdd_post_size - 0.25, 0)


%% shape oddball
freq_FirstFix            = [P_FirstOdd_pre_shape P_FirstOdd_post_shape];

fg5 = figure('position',[1024 369 689 497]);
pl1                 = plot_bar(freq_FirstFix,plt);
ylim([0 0.35])
format_figure(0.25,nan,'fixation type','probability of first fixation');
% format_legend([pl1(2).h(1) pl1(2).h(2)],{'size','shape'});
set(gca,'xtick',[1 2],'xticklabel',{'pre','post'})
tl = title('Shape'); tl.FontWeight = 'normal';


saveas(fg5,fullfile(save_path, 'first_fix_shape.jpg'));

ttest(P_FirstOdd_pre_shape - 0.25, 0)
ttest(P_FirstOdd_post_shape - 0.25, 0)

%% accumulated time on oddball / standard
fg6 = figure('Position',[329 182 1214 764]);

subplot(121)
points     = [Acc_Dur_SizeOdd_post Acc_Dur_SizeStd_post];
plot_violin(points,'x',1:2,'showdots',1,'color',colors,'rescaling',.25,'connline',[.5 .5 .5],'dotssize',60,'jitfactor',.05,'dotshift',-.02,'allline',10)
xlim([0.5 2.5]);%ylim([200 450]);%yticks([0.45:.1:.95])
format_figure(nan,nan,'','Accumulated Duration(ms)',[],[],20)
set(gca,'xtick',1:2,'xticklabel',{'Oddball','Standard'})
box on
tl = title('Size'); tl.FontWeight = 'normal';

subplot(122)
points     = [Acc_Dur_ShapeOdd_post Acc_Dur_ShapeStd_post];
plot_violin(points,'x',1:2,'showdots',1,'color',colors,'rescaling',.25,'connline',[.5 .5 .5],'dotssize',60,'jitfactor',.05,'dotshift',-.02,'allline',10)
xlim([0.5 2.5]);%ylim([200 450]);%yticks([0.45:.1:.95])
format_figure(nan,nan,'','Accumulated Duration(ms)',[],[],20)
set(gca,'xtick',1:2,'xticklabel',{'Oddball','Standard'})
box on
tl = title('Shape'); tl.FontWeight = 'normal';
sgtitle('Accumulated Fixation Duration postfix')

saveas(fg6,fullfile(save_path, 'Accumulated Fixation Duration post_fix.jpg'));

%% pre-fixation
fg7 = figure('Position',[329 182 1214 764]);

subplot(121)
points     = [Acc_Dur_SizeOdd_pre Acc_Dur_SizeStd_pre];
plot_violin(points,'x',1:2,'showdots',1,'color',colors,'rescaling',.25,'connline',[.5 .5 .5],'dotssize',60,'jitfactor',.05,'dotshift',-.02,'allline',10)
xlim([0.5 2.5]);%ylim([200 450]);%yticks([0.45:.1:.95])
format_figure(nan,nan,'','Accumulated Duration(ms)',[],[],20)
set(gca,'xtick',1:2,'xticklabel',{'Oddball','Standard'})
box on
tl = title('Size'); tl.FontWeight = 'normal';

subplot(122)
points     = [Acc_Dur_ShapeOdd_pre Acc_Dur_ShapeStd_pre];
plot_violin(points,'x',1:2,'showdots',1,'color',colors,'rescaling',.25,'connline',[.5 .5 .5],'dotssize',60,'jitfactor',.05,'dotshift',-.02,'allline',10)
xlim([0.5 2.5]);%ylim([200 450]);%yticks([0.45:.1:.95])
format_figure(nan,nan,'','Accumulated Duration(ms)',[],[],20)
set(gca,'xtick',1:2,'xticklabel',{'Oddball','Standard'})
box on
tl = title('Shape'); tl.FontWeight = 'normal';
sgtitle('Accumulated Fixation Duration prefix')

saveas(fg7,fullfile(save_path, 'Accumulated Fixation Duration pre_fix.jpg'));
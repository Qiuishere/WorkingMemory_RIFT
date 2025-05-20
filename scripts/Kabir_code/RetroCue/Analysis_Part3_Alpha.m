%% Dissociating External and Internal Attentional Selection | ANALYSIS
%  RIFT Analysis | RetroCue | part 3 of 5
%
%  This file
%  - computes alpha power across left and right-cued trials
%
%  Requirements
%  - list of participants (corresponding to directory names) at the start
%  - directory for preprocessed .mat file and saving in lines 25 and 119-120
%
%  Outputs (if saving is on, saved in saved/"participant")
%  - attend_left_new_3cyc_initbaseline: TFRs for left cued trials | (1 x channels x frequencies x timepoints)
%  - attend_right_new_3cyc_initbaseline: TFRs for right cued trials | (1 x channels x frequencies x timepoints)

%% Settings

clear all; close all; clc;
addpath("C:\Users\Arora003\Documents\[002] Data\InternalAttention_RetroCue_Project_64\Analysis\Helper_Functions");
addpath('C:\Program Files\MATLAB\R2022a\toolbox\fieldtrip-20220707'); ft_defaults;

participants = ["H0", "HE", "LI", "CM", "DB", "B0", "C0", "N0", "O0", ...
                "F0", "NE", "NA", "AL", "SI", "P0", "S0", "K0", "CA", ...
                "SC", "V0", "CR", "MN"];

dir = "D:\Data\RetroCue\";

% FFT info
fminmax = 1:321; padd_it = 1; do_evoked_fft = 1;
num_wavelets = 3; baseline_TFR = 1; bslntype = 'db';
tfr_baseline = [-3.2 -3]; timespan = -3.4:0.01:2;
freq_range = 8:0.2:13.5;

% Other info
correct_trials = 2; % 0 only incorrect, 1 only correct, 2 all
bad_chans = [{'T7'}, {'T8'}, {'Fp1'}, {'Fp2'}, {'TP7'}, {'TP8'}, {'AF7'}, {'AF3'}, {'Fpz'}, {'CP6'}, {'C6'}, {'C5'}]; % exclude frontal
f_s = 2048;

% Remove bad channels
Bchannel = {'EEG'};
for i = 1:length(bad_chans)
    Bchannel{end+1} = ['-' bad_chans{i}];
end

%% Compute Alpha-Power

for sub = 1:24
    tic
    disp("---- Working on participant "+num2str(sub)+" of "+num2str(length(subjects))+" ----");

    subject = subjects(sub);

    % load preprocessed data
    disp("Importing data")
    load(char(dir+subject+"\preprocessing\clean\"+subject+"_preprocessed.mat"));

    numtrls(sub) = size(data_temp_3.trialinfo,1);

    cfg = [];
    cfg.channel = Bchannel;
    data = ft_selectdata(cfg,data_temp_3);

    % correct answers
    data.trialinfo(:,15) = data.trialinfo(:,5)==data.trialinfo(:,13);
    
    % select correct/incorrect trials
    cfg = [];
    switch correct_trials
        case 0
            cfg.trials = find(data.trialinfo(:,15)==0); %incorrect
        case 1
            cfg.trials = find(data.trialinfo(:,15)==1); %correct
        case 2
            %all
    end
    data = ft_selectdata(cfg,data);

    % select attend left/right trials
    cfg = [];
    cfg.trials = logical(data.trialinfo(:,9))==0;
    dataattL = ft_selectdata(cfg,data);
    cfg.trials = logical(data.trialinfo(:,9))==1;
    dataattR = ft_selectdata(cfg,data);

    cfg = []; cfg.resamplefs = 512;
    dataattL = ft_resampledata(cfg, dataattL);
    dataattR = ft_resampledata(cfg, dataattR);
         
    % TFR settings
    cfg = [];
    cfg.method = 'wavelet';
    cfg.keeptrials = 'no';
    cfg.channel = {'EEG'}; %{'CP5','CP3','CP1','P5','P3','P1','PO3','CP6','CP4','CP2','P6','P4','P2','PO4'};
    cfg.foi = freq_range; % Hz
    cfg.toi = timespan; % sec
    cfg.width = num_wavelets;
    cfg.gwidth = 3;

    % do TFR
    disp("Running ft_freqanalysis")
    evalc('dataattL_freq = ft_freqanalysis(cfg, dataattL)');    
    evalc('dataattR_freq = ft_freqanalysis(cfg, dataattR)');

    if baseline_TFR
        cfg = [];
        cfg.baseline     = tfr_baseline;
        cfg.baselinetype = bslntype; %'relative', 'relchange' 'db' 
        [dataattL_freq] = ft_freqbaseline(cfg, dataattL_freq);
        [dataattR_freq] = ft_freqbaseline(cfg, dataattR_freq);
    end

    if 0
       cfg = [];
       ft_singleplotTFR(cfg,dataattL_freq); 
    end
      
    attendLeft(1, :, :, :) = dataattL_freq.powspctrm;
    attendRight(1, :, :, :) = dataattR_freq.powspctrm;

    disp("Saving AttendLeft info"); save([char(dir+"saved/"+subjects(sub)+"/attend_left_new_3cyc_initbaseline")],'attendLeft','-v7.3');    
    disp("Saving AttendRight info"); save([char(dir+"saved/"+subjects(sub)+"/attend_right_new_3cyc_initbaseline")],'attendRight','-v7.3');  

    % Code duration
    run_dur(sub) = toc;
    disp("This subject took "+num2str(floor(round(run_dur(sub))/60))+" min, "+ num2str(round(run_dur(sub)) - (floor(round(run_dur(sub))/60))*60) +" sec.");
    avg_time = mean(run_dur)*(length(subjects)-sub);
    disp("Expected time remaining is "+num2str(floor(round(avg_time)/60))+" min, "+ num2str(round(avg_time) - (floor(round(avg_time)/60))*60) +" sec.");
end

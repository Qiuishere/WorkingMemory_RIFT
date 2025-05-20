%% Dissociating External and Internal Attentional Selection | ANALYSIS
%  RIFT Analysis | RetroCue | part 1 of 5
%
%  This file
%  - converts a preprocessed FieldTrip dataset into a condition-wise matrix
%
%  Requirements
%  - list of participants (corresponding to directory names) at the start
%  - directory for preprocessed .mat file and saving in lines 27 and 87
%
%  Outputs (if saving is on, saved in saved/"participant")
%  - eeg_matrix: Preprocessed data in matrix format | {conditions x participants}(channels x timepoints x trials)

%% Settings

close all; clc; clear all;

% Fieldtrip 
addpath('C:\Program Files\MATLAB\R2022a\toolbox\fieldtrip-20220707'); ft_defaults;
addpath(genpath("C:\Users\Arora003\Documents\[002] Data\InternalAttention_RetroCue_Project_64")) % add helper functions
saving = 1;

participants = ["H0", "HE", "LI", "BE", "CM", "DB", "B0", "C0", "N0", "O0", ...
                "F0", "NE", "NA", "MG", "AL", "SI", "P0", "S0", "K0", "CA", ...
                "SC", "V0", "CR", "MN"];

dir = "D:\Data\RetroCue\";

%% Convert

for p = 1:length(participants)
    tic;

    % Folders
    if saving
        if ~isfolder(dir+"saved\"+participants(p))
            mkdir(dir+"saved\"+participants(p))
        end
    end
    
    % Collection info
    f_s = 2048; % Hz
    
    % Parameters
    baseline = [-3.4 -2.9];
    
    % Condition information
    conds = [1, 0; 1, 1; 0, 0; 0, 1]; % 60|68 AL, 60|68 AR, 68|60 AL, 68|60 AR
    data_eeg = cell(size(conds, 1), 1); % {conditions, 1}(channels, timepoints. trials)
    
    disp("--- Importing data from participant "+num2str(p)+" of "+num2str(length(participants))+" ---");

    % Import data
    temp_p = load(dir+participants(p)+"\preprocessing\clean\"+participants(p)+"_preprocessed.mat");
    temp_p = temp_p.data_temp_3;

    t = temp_p.time{1}; num_trials = length(temp_p.trialinfo);

    % Baseline correct
    disp("Baseline correcting");
    [~, b1] = min(abs(t-baseline(1)));
    [~, b2] = min(abs(t-baseline(2)));  
    for ii = 1:length(temp_p.trialinfo)
        temp_p.trial{ii} = ft_preproc_baselinecorrect(temp_p.trial{ii}, b1, b2);
    end
    data = temp_p;

    % Split data into four conditions
    for cond = 1:size(conds, 1)
        % Prepare cfg
        cfg = []; cfg.channel = 'EEG'; 
%         cfg.trials = find(all([(data.trialinfo(:, 8)' == conds(cond, 1));  ... % attend/flicker side
%                                (data.trialinfo(:, 9)' == conds(cond, 2)); ... % attend/flicker side
%                                (data.trialinfo(:, 5)' == data.trialinfo(:, 13)')])); % correct trials

        cfg.trials = find(all([(data.trialinfo(:, 8)' == conds(cond, 1));  ... % attend/flicker side
                               (data.trialinfo(:, 9)' == conds(cond, 2))])); % attend/flicker side      
        % Select trials
        evalc('temp_cond = ft_selectdata(cfg, data)');

        for tr = 1:length(temp_cond.trialinfo)
            data_eeg{cond, 1}(:, :, tr) = temp_cond.trial{tr}(1:64, :);
        end
        disp("Done with condition "+num2str(cond)+" of "+num2str(size(conds, 1)))
    end
    
    if saving; disp("Saving big EEG matrix"); save([char(dir+"saved\"+participants(p)+"\eeg_matrix_all")],'data_eeg','-v7.3'); end % ~2 GB per participant

    % Code duration
    run_dur(p) = toc;
    disp("This subject took "+num2str(floor(round(run_dur(p))/60))+" min, "+ num2str(round(run_dur(p)) - (floor(round(run_dur(p))/60))*60) +" sec.");
    avg_time = mean(run_dur)*(length(participants)-p);
    disp("Expected time remaining is "+num2str(floor(round(avg_time)/60))+" min, "+ num2str(round(avg_time) - (floor(round(avg_time)/60))*60) +" sec.");
end
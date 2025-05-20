%% Dissociating External and Internal Attentional Selection | ANALYSIS
%  RIFT Analysis | PreCue | part 2 of 5
%
%  This file
%  - computes coherence at several frequencies including the stimuli frequencies
%  - plots coherence topoplots across attentional conditions / time intervals
%
%  Requirements
%  - list of participants (corresponding to directory names) at the start
%  - directory for eeg_matrix and saving in lines 24, 74, and 89
%
%  Outputs (if saving is on, saved in saved/"participant")
%  - spectrogram_coherence_data_extended: Coherence Traces across different trial sets | (frequencies x 1 x conditions x channels x timepoints)

%% Settings
clear all; close all; clc
addpath(genpath("C:\Users\Arora003\Documents\[002] Data\InternalAttention_RetroCue_Project_64")) % add helper functions
addpath('C:\Program Files\MATLAB\R2022a\toolbox\fieldtrip-20220707'); ft_defaults;

participants = ["LH", "ME", "VE", "EA", "MA", "JU", "SA", "UR", "NE", "PL"...
                "MO", "PH", "DE", "GA", "CA", "EU", "IO", "TI", "RH", "IA"...
                "TT", "OB", "UM", "AR"];

dir = "D:\Data\PreCue\";

% Parameters
f_s = 2048; % Hz
dt = 1/f_s; t_f = 2+2.5; % sec, window of RetroCue data is [-2, 2.5] with cue onset at 0
freq = 56.8:0.8:67.2; % Hz
saving = 1;
time = (0:dt:t_f-dt)'; % sec
delta = 1.9; % Hz, bandpass filter width

%% Coherence

for p = 1:length(participants)

    % Import data
    tic; disp("*** "+participants(p)+" ***")
    load(char(dir+"saved/"+participants(p)+"/eeg_matrix_all"));
    
    % Reshape matrices
    disp("%%% Reshaping data matrix for coherence function %%%");
    data_eeg_ordered = cell(4, 1);

    for cond = 1:size(data_eeg, 1)
        data_eeg_ordered{cond, 1} = permute(data_eeg{cond, 1}, [1, 3, 2]); % (channels x trials x timepoints)
    end
    
    % Computing coherence
    disp("%%% Computing coherence %%%")
    
    coh_60_l = []; coh_68_l = []; coh_60_r = []; coh_68_r = [];
    coherence_data = zeros(length(freq), size(data_eeg, 2), 4, size(data_eeg_ordered{1, 1}, 1), size(data_eeg_ordered{1, 1}, 3)); % frequencies, participants, conditions, channels, timepoints
    
    for i = 1:length(freq)
        signal = sin(2*pi*freq(i)*time + 0.5*pi);
        for ch = 1:size(data_eeg_ordered{1, 1}, 1)
            temp = ft_preproc_bandpassfilter(squeeze(data_eeg_ordered{1, 1}(ch, :, :)), 2048, [freq(i)-delta freq(i)+delta], 4, 'but', 'twopass', [], [], 'hamming');
            coherence_data(i, 1, 1, ch, :) = coherence(temp', signal);
    
            temp = ft_preproc_bandpassfilter(squeeze(data_eeg_ordered{2, 1}(ch, :, :)), 2048, [freq(i)-delta freq(i)+delta], 4, 'but', 'twopass', [], [], 'hamming');
            coherence_data(i, 1, 2, ch, :) = coherence(temp', signal);        
    
            temp = ft_preproc_bandpassfilter(squeeze(data_eeg_ordered{3, 1}(ch, :, :)), 2048, [freq(i)-delta freq(i)+delta], 4, 'but', 'twopass', [], [], 'hamming');
            coherence_data(i, 1, 3, ch, :) = coherence(temp', signal);
    
            temp = ft_preproc_bandpassfilter(squeeze(data_eeg_ordered{4, 1}(ch, :, :)), 2048, [freq(i)-delta freq(i)+delta], 4, 'but', 'twopass', [], [], 'hamming');
            coherence_data(i, 1, 4, ch, :) = coherence(temp', signal);             
        end
        disp("Done with "+num2str(freq(i))+"Hz")
    end
    
    if saving; disp("Saving coherence data"); save([char(dir+"saved/"+participants(p)+"/spectrogram_coherence_data_extended_all_bp1")],'coherence_data','-v7.3'); end

    % Code duration
    run_dur(p) = toc;
    disp("This subject took "+num2str(floor(round(run_dur(p))/60))+" min, "+ num2str(round(run_dur(p)) - (floor(round(run_dur(p))/60))*60) +" sec.");
    avg_time = mean(run_dur)*(length(participants)-p);
    disp("Expected time remaining is "+num2str(floor(round(avg_time)/60))+" min, "+ num2str(round(avg_time) - (floor(round(avg_time)/60))*60) +" sec.");
end

%% Import Coherence

% Import Coherence Data
coherence_data_imported = []; 
for p = 1:length(participants)
    disp("Importing from participant "+num2str(p)+" of "+num2str(length(participants)));
    load([char(dir+"saved\"+participants(p)+"\spectrogram_coherence_data_extended_all.mat")]);
    coherence_data_imported(:, p, :, :, :) = coherence_data(:, 1, :, :, :);
                % ^This is frequencies x participants x conditions x channels x timepoints
end

load('C:\Users\Arora003\Documents\[002] Data\InternalAttention_PreCue_Project\Analysis\saved\t_short.mat'); t = tm2; t = t - 2;
load('C:\Users\Arora003\Documents\[002] Data\InternalAttention_RetroCue_Project_64\Analysis\saved\skymap.mat')
freq_labels = 56.8:0.8:67.2; freq_idxs = [5, 10]; % 60 at index 5, 64 at index 10

%% Coherence spectrogram, averaged everything 

close all
spectrogram_avgd_data = squeeze(mean(coherence_data_imported, [2, 3, 4]));
colormap(skymap)

imagesc('XData', t, 'YData', freq_labels, 'CData', squeeze(spectrogram_avgd_data(:, :)));
xlim([-0.5, 2.5]); ylim([freq_labels(1), freq_labels(end)]); xline([-1, 0, 1], 'w--'); yline([60, 64], 'k-', 'LineWidth', 2)
xlabel("Time w.r.t. stimuli onset"); ylabel("Frequency (Hz)"); colorbar; 

%% Coherence topoplot, averaged everything + condition-wise
close all
topoplot_avgd_data = squeeze(mean(coherence_data_imported, [2, 3, 5]));
topoplot_conditions_data = squeeze(mean(coherence_data_imported, [2, 5]));

% Overall
figure;
for f = 1:2    
    sp = subplot(1, 2, f);
    [cfg, topotemplate] = plot_topotemplate(topoplot_avgd_data(freq_idxs(f), :), ''); 
    cfg.colormap = 'skymap'; cfg.figure = sp; cfg.zlim = [0 0.03];
    ft_topoplotER(cfg,topotemplate);
end

% Condition-wise
for f = 1:2  
    figure;
    for cond = 1:4
        sp = subplot(2, 2, cond);
        [cfg, topotemplate] = plot_topotemplate(topoplot_conditions_data(freq_idxs(f), cond, :), ''); 
        cfg.colormap = 'skymap'; cfg.figure = sp; cfg.zlim = [0 0.03];
        ft_topoplotER(cfg,topotemplate);
    end
end

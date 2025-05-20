%% Dissociating External and Internal Attentional Selection | ANALYSIS
%  RIFT Analysis | RetroCue | part 4 of 5
%
%  This file
%  - prepares various data types for an LME
%  - runs LMEs to explore predictors of the RIFT response and Alpha lateralization
%
%  Requirements
%  - list of participants (corresponding to directory names) at the start
%  - directory for importing in line 16
%
%  Outputs (if saving is on, saved in saved/"participant")
%  - 

%% Settings
dir = "D:\Data\RetroCue\";

% Participant Info
participants = ["H0", "HE", "LI", "BE", "CM", "DB", "B0", "C0", "N0", "O0", ...
                "F0", "NE", "NA", "MG", "AL", "SI", "P0", "S0", "K0", "CA", ...
                "SC", "V0", "CR", "MN"];

% Collection info
f_s_eye = 250; % Hz
f_s_eeg = 2048; % Hz

% Screen parameters
dist_from_screen = 72; screen_dim = [48 27.2]; % cm
screen_pxs = [1920 1080];
screen_dva = rad2deg(2*atan((0.5*screen_dim)/dist_from_screen));
dva_per_px = screen_dva(1)/screen_pxs(1);

%% Import (1/4) - Task data

load('C:\Users\Arora003\Documents\[002] Data\InternalAttention_RetroCue_Project_64\Analysis\saved\t.mat'); % timeline, -3.4 to 2s
attend_sides_mat = []; flicker_sides_mat = [];

for p = 1:length(participants)
    % Import good trialnumbers
    load(char("D:\Data\RetroCue\"+participants(p)+"\preprocessing\clean\"+participants(p)+"_badtrial_idxs"));
    trials_logical = ones(1, 490); trials_logical(1:10) = 0; trials_logical(removed_trials+10) = 0; % remove 10 practice trials + bad EEG trials

    exp_info = readtable(char("D:\Data\ToUpload\Dissociating_External_and_Internal - Eyetracking Data\RetroCue\"+participants(p)+"\trial_info.csv"));
    att_sides_temp = exp_info.attend_sides; fli_sides_temp = exp_info.flicker_sides;

    % Select preprocessed trials
    attend_sides_mat{p} = exp_info.attend_sides(find(trials_logical));
    flicker_sides_mat{p} = exp_info.flicker_sides(find(trials_logical));
end

%% Import (2/4) - Timepoint-by-timepoint RIFT hilbert magnitude

magnitude_60 = []; magnitude_64 = [];

for p = 1:length(participants)
    temp = load(char(dir+"saved/"+participants(p)+"/hilbert_mags_continuous"));

    magnitude_60{p} = squeeze(mean(temp.magnitudes(1, :, :, :), 2));
    magnitude_64{p} = squeeze(mean(temp.magnitudes(2, :, :, :), 2));

    disp("Done with participant "+num2str(p))
end

% Coherence traces attended and unattended, 60 and 64
conds = [1, 0; 1, 1; 0, 0; 0, 1]; % 60 < 64. 60 > 64, 64 < 60, 64 > 60;
magnitude_traces = zeros(2, length(participants), 4, size(magnitude_64{1}, 2)); % freqs x participants x conds x timepoints
attended_magnitude = []; unattended_magnitude = []; rift_interval = [0.23, 0.72];

for p = 1:length(participants)
    for cond = 1:4
        fli_idxs = flicker_sides_mat{p} == conds(cond, 1);
        att_idxs = attend_sides_mat{p} == conds(cond, 2);

        magnitude_traces(1, p, cond, :) = squeeze(mean(magnitude_60{p}((fli_idxs & att_idxs), :), 1));
        magnitude_traces(2, p, cond, :) = squeeze(mean(magnitude_64{p}((fli_idxs & att_idxs), :), 1));
    end

    % One value per trial - use significant interval from group data
    attended_magnitude{p} = zeros(1, size(magnitude_64{p}, 1));
    att_60_idxs = xor(flicker_sides_mat{p}, attend_sides_mat{p}); att_64_idxs = ~att_60_idxs;
    t_range = dsearchn(t', rift_interval(1)):dsearchn(t', rift_interval(2));
    
    attended_magnitude{p}(att_60_idxs) = squeeze(mean(magnitude_60{p}(att_60_idxs, t_range), 2, "omitnan"));
    attended_magnitude{p}(att_64_idxs) = squeeze(mean(magnitude_64{p}(att_64_idxs, t_range), 2, "omitnan"));

    unattended_magnitude{p}(att_64_idxs) = squeeze(mean(magnitude_60{p}(att_64_idxs, t_range), 2, "omitnan"));
    unattended_magnitude{p}(att_60_idxs) = squeeze(mean(magnitude_64{p}(att_60_idxs, t_range), 2, "omitnan"));    
end

%% Import (3/4) - Gaze Position data

t_eye = (-1000:4:1600)/1000;
gazebias = []; x_pos_mat = []; gaze_interval = [0.39, 1]; % Interval significant at group level
[participants_alph, participants_alph_idxs] = sort(participants); % Gaze data is stored with alphabetized participant list

for p = 1:length(participants)
    % Import good trialnumbers
    load(char("D:\Data\RetroCue\"+participants(p)+"\preprocessing\clean\"+participants(p)+"_badtrial_idxs"));
    trials_logical = ones(1, 490); trials_logical(1:10) = 0; trials_logical(removed_trials+10) = 0; % remove 10 practice trials + bad EEG trials

    % Import gaze data
    p_fr = find(participants(p) == participants_alph); % Find alphabetized participant index
    load(char("C:\Users\Arora003\Documents\[002] Data\Re_ gaze bias re-analysis\saved_data\trialwise_gazepos_pp"+num2str(floor(p_fr/10))+num2str(mod(p_fr, 10))));
    temp = x_pos_bl(trials_logical == 1, :); temp = temp.*repmat(((attend_sides_mat{p}-0.5)*2), 1, length(t_eye)); % Flip direction of attend right trials
    x_pos_mat{p} = temp;

    % Select good trials for time-averaged metric
    gazebias{p} = squeeze(mean(x_pos_bl(trials_logical == 1, dsearchn(t_eye', gaze_interval(1)):dsearchn(t_eye', gaze_interval(2))), 2, "omitnan"));
    gazebias{p} = gazebias{p}.*((attend_sides_mat{p}-0.5)*2);
end

mean_gazes = zeros(1, length(participants));
for p = 1:length(participants)
    mean_gazes(1, p) = mean(gazebias{p}, "omitnan");
end

%% Import (4/4) - Alpha Power data, convert into trialwise lateralizations

load("C:\Users\Arora003\Documents\[002] Data\InternalAttention_RetroCue_Project\Analysis\ReAnalyzed\freq_range_modified_new.mat"); 
timespan = -3.4:0.01:2; 

alpha_traces = [];
for p = 1:length(participants)
    disp("Importing from participant "+num2str(p)+" of "+num2str(length(participants)));
    load([char(dir+"saved/"+participants(p)+"/attend_left_new_3cyc_initbaseline_trialwise")]);
    load([char(dir+"saved/"+participants(p)+"/attend_right_new_3cyc_initbaseline_trialwise")]);
    alpha_traces_l{p} = squeeze(mean(attendLeft, 4));
    alpha_traces_r{p} = squeeze(mean(attendRight, 4));
end

load("C:\Users\Arora003\Documents\[002] Data\InternalAttention_RetroCue_Project\Analysis\ReAnalyzed\sig_neg_chans_alpha.mat");
load("C:\Users\Arora003\Documents\[002] Data\InternalAttention_RetroCue_Project\Analysis\ReAnalyzed\sig_pos_chans_alpha.mat");

alpha_trialwise_lat = [];
baseline_t = dsearchn(timespan', -0.6):dsearchn(timespan', -0.1);
alpha_t = dsearchn(timespan', 0.31):dsearchn(timespan', 1.05);

for p = 1:length(participants)
    att_l_trials = logical(attend_sides_mat{p} == 0);
    att_r_trials = logical(attend_sides_mat{p} == 1);

    alpha_trialwise_lat{p} = zeros(size(alpha_traces_l{p}, 2)+size(alpha_traces_r{p}, 2), 1);

    alpha_trialwise_lat{p}(att_l_trials) = squeeze(mean(alpha_traces_l{p}(:, sig_pos_chans, alpha_t), [2, 3])) - squeeze(mean(alpha_traces_l{p}(:, sig_neg_chans, alpha_t), [2, 3]));
    alpha_trialwise_lat{p}(att_r_trials) = squeeze(mean(alpha_traces_r{p}(:, sig_neg_chans, alpha_t), [2, 3])) - squeeze(mean(alpha_traces_r{p}(:, sig_pos_chans, alpha_t), [2, 3]));    
end

% Mean lateralization per participant
for p = 1:length(participants)
    mean_lateralization(p) = mean(alpha_trialwise_lat{p}, 1);
end

%% GLME - RIFT Magnitude

clc; 
glme_riftmagnitude = []; x_1_all = []; x_2_all = []; x_3_all = []; x_4_all = []; x_5_all = []; y_all = [];

for p = 1:length(participants)
    % y - Magnitudes from coherence
    y = [attended_magnitude{p}'; unattended_magnitude{p}'];

    % (1/5) - Cued or uncued?
    x_1 = [ones(1, length(attended_magnitude{p}))'; zeros(1, length(attended_magnitude{p}))'];

    % (2/5) - Is the magnitude measured from 60 (0) or 64 (1)?
    att_60_idxs = xor(flicker_sides_mat{p}, attend_sides_mat{p}); att_64_idxs = ~att_60_idxs;
    x_2 = [att_64_idxs; ~att_64_idxs];

    % (3/5) - Towardness
    x_3 = [gazebias{p}; -gazebias{p}];

    % (4/5) - Alpha Lateralization
    x_4 = [alpha_trialwise_lat{p}; -alpha_trialwise_lat{p}];

    % (5/5) - Participant
    x_5 = p*ones(1, 2*length(attended_magnitude{p}))';

    % Prepare matrices
    x_1_all = [x_1_all; x_1]; x_2_all = [x_2_all; x_2]; x_3_all = [x_3_all; x_3]; x_4_all = [x_4_all; x_4]; x_5_all = [x_5_all; x_5]; y_all = [y_all; y];
end

% Convert to tables
table_for_glme = array2table([logical(x_1_all), logical(x_2_all), x_3_all, x_4_all, x_5_all, y_all], 'VariableNames', {'CuedUncued', 'Frequency', 'Towardness', 'AlphaLat', 'Participant', 'Magnitude'});
table_for_glme.CuedUncued = logical(table_for_glme.CuedUncued); table_for_glme.Frequency = logical(table_for_glme.Frequency);
glme_riftmagnitude = fitglme(table_for_glme, 'Magnitude ~ 1 + CuedUncued + Frequency + Towardness + AlphaLat + (1 + CuedUncued + Frequency + Towardness + AlphaLat|Participant)', 'CheckHessian', true);
disp(glme_riftmagnitude)

%% GLME - Alpha Lateralization

clc; 
glme_alphalat = []; x_1_all = []; x_2_all = []; x_3_all = []; x_4_all = []; x_5_all = []; y_all = [];

for p = 1:length(participants)
    % y - Alpha
    y = [alpha_trialwise_lat{p}];

    % (1/3) - Towardness
    x_1 = [gazebias{p}];

    % (2/3) - RIFT Magnitude
    x_2 = [attended_magnitude{p}'];

    % (5/5) - Participant
    x_3 = p*ones(1, length(attended_magnitude{p}))';

    % Prepare matrices
    x_1_all = [x_1_all; x_1]; x_2_all = [x_2_all; x_2]; x_3_all = [x_3_all; x_3]; y_all = [y_all; y];
end

% Convert to tables
table_for_glme = array2table([x_1_all, x_2_all, x_3_all, y_all], 'VariableNames', {'Towardness', 'Magnitude', 'Participant', 'AlphaLat'});
glme_alphalat = fitglme(table_for_glme, 'AlphaLat ~ 1 + Magnitude + Towardness + (1 + Magnitude + Towardness|Participant)', 'CheckHessian', true);
disp(glme_alphalat)

%% Figure

% Extract coefficient info
coeffs_nointercept = glme_riftmagnitude.Coefficients([2:end], :);

close all; figure('units', 'normalized', 'outerposition', [0 0.1 0.5 0.7]);
ax_curr = axes('FontName', 'Nexa-light', 'Position', [0.25 0.5 0.7 0.4]); axes(ax_curr);

plot_glm(coeffs_nointercept, [-0.065, 0.04], 0, {'Cued or Uncued'; 'Frequency (60/64)'; 'Gaze Bias'; 'Alpha Lateralization'});

% Extract coefficient info
coeffs_nointercept = glme_alphalat.Coefficients([2:end], :);

%close all; figure('units', 'normalized', 'outerposition', [0 0.1 0.5 0.5*0.7]);
ax_curr = axes('FontName', 'Nexa-light', 'Position', [0.25 0.2 0.7 0.2]); axes(ax_curr);

plot_glm(coeffs_nointercept, [-0.6, 0.602], 1, {'Gaze Bias'; 'RIFT Response'});

















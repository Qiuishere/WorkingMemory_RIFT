%% Dissociating External and Internal Attentional Selection | ANALYSIS
%  RIFT Analysis | RetroCue | part 5 of 5
%
%  This file
%  - imports and preprocesses gaze position data (blink_correction, saccade removal, baselining)
%
%  Requirements
%  - list of participants (corresponding to directory names) at the start
%  - directory for eeg_matrix and saving in lines 24, 
%
%  Outputs (if saving is on, saved in saved/"participant")
%  - 

%% Settings

close all; clear all; clc;
addpath(genpath("C:\Users\Arora003\Documents\[002] Data\InternalAttention_RetroCue_Project_64")) % add data and helper functions

% Participant Info
participants = ["H0", "HE", "LI", "BE", "CM", "DB", "B0", "C0", "N0", "O0", ...
                "F0", "NE", "NA", "MG", "AL", "SI", "P0", "S0", "K0", "CA", ...
                "SC", "V0", "CR", "MN"];

dir = "D:\Data\RetroCue\";

% Collection info
f_s = 250; % Hz
trial_latency = [-1, +4.55]; % sec, has to include 0

% Screen info
dist_from_screen = 72; screen_dim = [48 27.2]; % cm
screen_pxs = [1920 1080];
screen_dva = rad2deg(2*atan((0.5*screen_dim)/dist_from_screen));
dva_per_px = screen_dva(1)/screen_pxs(1);

%% Import eyetracking data
disp("(1/7) Importing data from edf files");

for p = 1:length(participants)
    disp("Participant "+num2str(p)+" of "+num2str(length(participants)));

    % Import data
    filename = dir+participants(p)+"\\Eyetracking\\"+participants(p)+".edf";
    edf_data{p} = Edf2Mat(char(filename));
end

% Convert into a big matrix
disp("(2/7) Converting edf into a matrix");

x_data_blinkcorrected = []; y_data_blinkcorrected = [];

% Trial number info
for p = 1:length(participants)
    triggers{p} = find(strncmpi(edf_data{1, p}.Events.Messages.info, "Trial", 5)); % ignore triggers during calibration
    num_trials(p) = length(triggers{p});
    sampling_freqs{p} = edf_data{1, p}.RawEdf.RECORDINGS(1).sample_rate;

    for eye = 1:2
        if (participants(p) == "CR") || (participants(p) == "MN") % only left eye
            x_data_blinkcorrected{1, p}(eye, :) = draft_blink_edge_correction(edf_data{1, p}.Samples.pupilSize(:, 1)', edf_data{1, p}.Samples.posX(:, 1)', f_s);
            y_data_blinkcorrected{1, p}(eye, :) = draft_blink_edge_correction(edf_data{1, p}.Samples.pupilSize(:, 1)', edf_data{1, p}.Samples.posY(:, 1)', f_s);
        else
            x_data_blinkcorrected{1, p}(eye, :) = draft_blink_edge_correction(edf_data{1, p}.Samples.pupilSize(:, eye)', edf_data{1, p}.Samples.posX(:, eye)', f_s);
            y_data_blinkcorrected{1, p}(eye, :) = draft_blink_edge_correction(edf_data{1, p}.Samples.pupilSize(:, eye)', edf_data{1, p}.Samples.posY(:, eye)', f_s);
        end
    end
    disp("Done blink correction for participant "+num2str(p)+" of "+num2str(length(participants)));
end

% Are all trial numbers same?
if nnz(diff(num_trials))
    disp("Warning: Unequal numbers of trials across participants");
end

%% Prepare variables
x_pos_unmatched = zeros(num_trials(1), floor(sum(abs(trial_latency))*f_s), 2, length(participants)); % both eyes
y_pos_unmatched = x_pos_unmatched; % trials x timepoints x eyes x participants

% Extract position traces and store
for p = 1:length(participants)
    f_s_temp = sampling_freqs{p};

    for tr = 1:num_trials(1)
        trig = triggers{p}(tr);
    
        % Start time as per trigger
        [~, trig_time_idx(tr)] = min(abs(edf_data{1, p}.Samples.time - edf_data{1, p}.Events.Messages.time(1, trig)));
    
        % Start and end times
        start_t = trig_time_idx(tr) + floor(trial_latency(1)*f_s_temp);
        end_t = start_t + floor(sum(abs(trial_latency))*f_s_temp) - 1;

        % Extract position with either original trace (250Hz) or downsampled trace (500Hz)
        x_pos_temp = []; y_pos_zeros = [];
        if f_s_temp == 250
            x_pos_temp = x_data_blinkcorrected{p}(:, start_t:end_t)';
            y_pos_temp = y_data_blinkcorrected{p}(:, start_t:end_t)';
        elseif f_s_temp == 500
            for eye = 1:2
                x_pos_temp(:, eye) = downsample(x_data_blinkcorrected{p}(eye, start_t:end_t-1), 2);
                y_pos_temp(:, eye) = downsample(y_data_blinkcorrected{p}(eye, start_t:end_t-1), 2);
            end
        end        
    
        % Store position
        x_pos_unmatched(tr, :, :, p) = x_pos_temp;
        y_pos_unmatched(tr, :, :, p) = y_pos_temp;     

    end
end

%% Mark trials with saccades
disp("(3/7) Marking trials with saccades");

% Parameters
screen_pxs = [1920 1080]; num_trials = 490;
saccade_threshold = 2; %dva
x_deviation = []; y_deviation = []; eccentricities_px = []; saccade_thresh_dur = 0.05;

% Identify timepoints exceeding threshold
for p = 1:length(participants)
    for eye = 1:2
        x_deviation(:, :, eye, p) = x_pos_unmatched(:, [floor(1*f_s):floor(2*f_s) floor(3.4*f_s):floor(5*f_s)], eye, p) - mean(x_pos_unmatched(:, :, eye, p), "all", "omitnan");
        y_deviation(:, :, eye, p) = y_pos_unmatched(:, [floor(1*f_s):floor(2*f_s) floor(3.4*f_s):floor(5*f_s)], eye, p) - mean(y_pos_unmatched(:, :, eye, p), "all", "omitnan");        

        eccentricities_px(:, :, eye, p) = sqrt(x_deviation(:, :, eye, p).^2 + y_deviation(:, :, eye, p).^2);
    end
end
points_beyond_threshold = eccentricities_px > (saccade_threshold/dva_per_px);

% Average eyes
eyes_averaged = mean(points_beyond_threshold, 3);
saccade_trials = zeros(num_trials(1), length(participants));

for p = 1:length(participants)
    for tr = 1:num_trials(1)
        saccade_trials(tr, p) = (nnz(eyes_averaged(tr, :, p)) > saccade_thresh_dur*f_s);        
    end
    disp("--- Participant "+num2str(p)+" made saccades in "+num2str(nnz(saccade_trials(:, p)/num_trials(1)))+" trials out of "+num_trials(1))
end

%% Cue-matching of stimulus onset-locked data
disp("(4/7) Cue-matching the data");

% Portion to select
ccl = [-3, 1.5]; % cue-centered latency, sec, must include 0
num_trials = 490*ones(1, 24);

% Variables for cue-matched data
x_pos_uncorrected = zeros(num_trials(1), floor(sum(abs(ccl))*f_s) + 1, 2, length(participants));
y_pos_uncorrected = x_pos_uncorrected; % trials x timepoints x eyes x participants

for p = 1:length(participants)
    % Load collection info
    load(dir+participants(p)+"\\Experimental\\data.mat");

    % Cue delays per trial
    cue_delays = [data.cue_delay(1,1:10) reshape(data.cue_delay(2:end,:)', [], 1)'];
    post_trigger_time_sec = -trial_latency(1) + 2.5 + cue_delays; % 2.5 sec between stimuli onset and base maintenance end
    post_trigger_time_eytckr = floor(post_trigger_time_sec*f_s); % Convert to samples  

    % Direction of attention
    temp = [data.attend_sides(1,1:10) reshape(data.attend_sides(2:end,:)', [], 1)'];
    attend_sides(p,:) = temp(1:num_trials(1));

    % Extract post-cue portion
    for tr = 1:num_trials
        x_pos_uncorrected(tr, :, :, p) = x_pos_unmatched(tr, post_trigger_time_eytckr(tr) + (ccl(1)*f_s):post_trigger_time_eytckr(tr)+(ccl(2)*f_s), :, p);
        y_pos_uncorrected(tr, :, :, p) = y_pos_unmatched(tr, post_trigger_time_eytckr(tr) + (ccl(1)*f_s):post_trigger_time_eytckr(tr)+(ccl(2)*f_s), :, p);
    end
end

%% Mark trials with NaNs/blinks
disp("(3/7) Marking trials with NaNs");

% Parameters
screen_pxs = [1920 1080];
saccade_threshold = 2.5; %dva
nan_free_period = [0, 0.4] + (-ccl(1)); % cue centered;

% Identify timepoints exceeding threshold
for p = 1:length(participants)
    for tr = 1:size(x_pos_uncorrected, 1)
        nan_or_not = anynan(x_pos_uncorrected(tr, floor(nan_free_period(1)*f_s):floor(nan_free_period(2)*f_s), :, p));    
        nan_trials(tr, p) = nan_or_not;
    end
    disp("--- Participant "+num2str(p)+" had NaNs in "+num2str(nnz(nan_trials(:, p)/num_trials(p)))+" trials out of "+num_trials(p))

    nan_trials_p = nan_trials(:, p);
end

%% NaN Corrections

min_nan_length = 0.07; % sec, patches < this will be interpolated

% Prepare variables
x_pos_unbaselined = x_pos_uncorrected; y_pos_unbaselined = y_pos_uncorrected;
disp("(5/7) Interpolating short periods of NaNs")

% Interpolate patches that are small enough
fix_counter = 0;
for p = 1:length(participants)
    counter_sh = 0; counter_lo = 0;
    for tr = 1:num_trials
        nans_logical = isnan(squeeze(x_pos_uncorrected(tr, :, :, 1)));

        for eye = 1:2
            state_changes = diff([0 nans_logical(:, eye)' 0]); % when this is 1/-1, a NaN block has started/ended 
    
            % Find NaN start/ends
            idxs_NaN_starts{tr} = find(state_changes==1); 
            idxs_NaN_ends{tr} = find(state_changes==-1);

            % If NaNs are present,
            if ~isempty(idxs_NaN_starts{tr})
                % then for each NaN block,
                for i = 1:length(idxs_NaN_starts{tr})
                    % which is within acceptable range,
                    block_length = idxs_NaN_ends{tr}(i)-idxs_NaN_starts{tr}(i);

                    if block_length < min_nan_length*f_s
                        % interpolate the block.
                        fix_counter = fix_counter + 1;

                        if idxs_NaN_starts{tr}(i) == 1
                            % if trial start is already a NaN, interpolate with block end point
                            x_pos_unbaselined(tr, idxs_NaN_starts{tr}(i):idxs_NaN_ends{tr}(i)-1, eye, p) = ...
                                linspace(x_pos_uncorrected(tr, idxs_NaN_ends{tr}(i), eye, p), x_pos_uncorrected(tr, idxs_NaN_ends{tr}(i), eye, p), block_length);
                            y_pos_unbaselined(tr, idxs_NaN_starts{tr}(i):idxs_NaN_ends{tr}(i)-1, eye, p) = ...
                                linspace(y_pos_uncorrected(tr, idxs_NaN_ends{tr}(i), eye, p), y_pos_uncorrected(tr, idxs_NaN_ends{tr}(i), eye, p), block_length);

                        elseif idxs_NaN_ends{tr}(i) == length(state_changes)
                            % if trial end is a NaN, interpolate with block start point
                            x_pos_unbaselined(tr, idxs_NaN_starts{tr}(i):idxs_NaN_ends{tr}(i)-1, eye, p) = ...
                                linspace(x_pos_uncorrected(tr, idxs_NaN_starts{tr}(i)-1, eye, p), x_pos_uncorrected(tr, idxs_NaN_starts{tr}(i)-1, eye, p), block_length);
                            y_pos_unbaselined(tr, idxs_NaN_starts{tr}(i):idxs_NaN_ends{tr}(i)-1, eye, p) = ...
                                linspace(y_pos_uncorrected(tr, idxs_NaN_starts{tr}(i)-1, eye, p), y_pos_uncorrected(tr, idxs_NaN_starts{tr}(i)-1, eye, p), block_length);

                        else
                            % if both start and end point are available
                            x_pos_unbaselined(tr, idxs_NaN_starts{tr}(i)-1:idxs_NaN_ends{tr}(i), eye, p) = ...
                                linspace(x_pos_uncorrected(tr, idxs_NaN_starts{tr}(i)-1, eye, p), x_pos_uncorrected(tr, idxs_NaN_ends{tr}(i), eye, p), block_length+2);
                            y_pos_unbaselined(tr, idxs_NaN_starts{tr}(i)-1:idxs_NaN_ends{tr}(i), eye, p) = ...
                                linspace(y_pos_uncorrected(tr, idxs_NaN_starts{tr}(i)-1, eye, p), y_pos_uncorrected(tr, idxs_NaN_ends{tr}(i), eye, p), block_length+2);

                        end
                    end
                end
            end           
        end
    end
end


%% Baseline Correction
disp("(6/7) Baseline correcting and converting to dva")

% Baseline period
baseline_int = [-1 -0.2] - ccl(1); % sec, wrt. cue onset, + cue latency(1)

% Prepare variables
x_pos_px = x_pos_unbaselined; y_pos_px = y_pos_unbaselined;

x_pos_px = x_pos_unbaselined - repmat(mean(x_pos_unbaselined(:, floor(baseline_int(1)*f_s):floor(baseline_int(2)*f_s), :, :), 2, "omitnan"), 1, size(x_pos_unbaselined, 2), 1, 1);
y_pos_px = y_pos_unbaselined - repmat(mean(y_pos_unbaselined(:, floor(baseline_int(1)*f_s):floor(baseline_int(2)*f_s), :, :), 2, "omitnan"), 1, size(y_pos_unbaselined, 2), 1, 1);

% DVA conversion
x_pos_dva = x_pos_px.*dva_per_px;
y_pos_dva = y_pos_px.*dva_per_px;


% Note trials with outlier baselines
baseline_outlier_trials = [];
for p = 1:length(participants)
    baselines = mean(x_pos_unbaselined(:, floor(baseline_int(1)*f_s):floor(baseline_int(2)*f_s), :, :), 2, "omitnan");
    baselines_avgd = squeeze(mean(baselines(:, :, :, p), 3, "omitnan"));
    baselines_y = mean(y_pos_unbaselined(:, floor(baseline_int(1)*f_s):floor(baseline_int(2)*f_s), :, :), 2, "omitnan");
    baselines_avgd_y = squeeze(mean(baselines_y(:, :, :, p), 3, "omitnan"));    

    % Find stdev
    threshold = 2*std(baselines_avgd, "omitnan"); avg_base = mean(baselines_avgd, "omitnan");
    threshold_y = 2*std(baselines_avgd_y, "omitnan"); avg_base_y = mean(baselines_avgd_y, "omitnan");

    % Trial idxs
    baseline_outlier = ((baselines_avgd) > (avg_base + threshold)) | ((baselines_avgd) < (avg_base - threshold));
    baseline_outlier_y = ((baselines_avgd_y) > (avg_base_y + threshold_y)) | ((baselines_avgd_y) < (avg_base_y - threshold_y));

    baseline_outlier_trials(:, p) = baseline_outlier | baseline_outlier_y;
end


%% Convert to DVA
disp("(7/7) Averaging trials without saccades")

% prepare variables
x_pos_dva_trlavgd = zeros([2, size(x_pos_dva, [2 3 4])]); y_pos_dva_trlavgd = x_pos_dva_trlavgd;
x_pos_dva_p_avgd = zeros([2, size(x_pos_dva, [2 3])]); y_pos_dva_p_avgd = x_pos_dva_p_avgd;
x_pos_dva_no_sac = []; y_pos_dva_no_sac = [];

avoid_trials = saccade_trials | baseline_outlier_trials;

for p = 1:length(participants)
    disp("--- Participant "+num2str(p)+", omitting "+num2str(nnz(avoid_trials(:, p)/num_trials(p)))+" trials out of "+num_trials(p))

    for att = 0:1 % left or right
        % Select trials
        accepted_trials = all([attend_sides(p, :)==att; ~(avoid_trials(:, p)')]);
        x_pos_dva_no_sac{p} = x_pos_dva(accepted_trials, :, :, p); y_pos_dva_no_sac{p} = y_pos_dva(accepted_trials, :, :, p);

        % Average
        x_pos_dva_trlavgd(att+1, :, :, p) = squeeze(mean(x_pos_dva_no_sac{p}, 1, "omitnan"));         
        y_pos_dva_trlavgd(att+1, :, :, p) = squeeze(mean(y_pos_dva_no_sac{p}, 1, "omitnan")); 
    end
end

eyes_to_use = zeros(1, length(participants));
x_pos_dva_e_avgd = []; y_pos_dva_e_avgd = [];
% Average eyes
for p = 1:length(participants)
    switch eyes_to_use(p)
        case 0
            x_pos_dva_e_avgd(:, :, p) = squeeze(mean(x_pos_dva_trlavgd(:, :, :, p), 3)); y_pos_dva_e_avgd(:, :, p) = squeeze(mean(y_pos_dva_trlavgd(:, :, :, p), 3));
        case 1
            x_pos_dva_e_avgd(:, :, p) = squeeze(x_pos_dva_trlavgd(:, :, 1, p)); y_pos_dva_e_avgd(:, :, p) = squeeze(y_pos_dva_trlavgd(:, :, 1, p));
        case 2
            x_pos_dva_e_avgd(:, :, p) = squeeze(x_pos_dva_trlavgd(:, :, 2, p)); y_pos_dva_e_avgd(:, :, p) = squeeze(y_pos_dva_trlavgd(:, :, 2, p));
    end
end

% Filter
filt_duration = 0.025; % s

win = floor(filt_duration*f_s); b = (1/win)*ones(1,win); a = 1;
x_pos_dva_e_avgd_filt = filter(b, a, x_pos_dva_e_avgd, [], 2);
y_pos_dva_e_avgd_filt = filter(b, a, y_pos_dva_e_avgd, [], 2);

%% Participant Average - x

% exclude 8, 14 based on saccades
pts = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]; %1:length(participants);

x_pos_dva_p_avgd = squeeze(mean(x_pos_dva_e_avgd_filt(:, :, pts), 3)); y_pos_dva_p_avgd = squeeze(mean(y_pos_dva_e_avgd_filt(:, :, pts), 3)); 
x_pos_dva_p_std = squeeze(std(x_pos_dva_e_avgd_filt(:, :, pts), [], 3)); y_pos_dva_p_std = squeeze(std(y_pos_dva_e_avgd_filt(:, :, pts), [], 3)); 

% Difference traces
x_pos_diff_trace = squeeze(x_pos_dva_e_avgd_filt(1, :, pts) - x_pos_dva_e_avgd_filt(2, :, pts));
y_pos_diff_trace = squeeze(y_pos_dva_e_avgd_filt(1, :, pts) - y_pos_dva_e_avgd_filt(2, :, pts));

% CIs for time series and differences
ci = [];
for att = 1:2
    ci(att, :, :) = bootci(5000, @mean, squeeze(x_pos_dva_e_avgd_filt(att, :, pts))')';
end

ci_diff = bootci(5000, @mean, squeeze(x_pos_diff_trace'));

% Cluster test
sig_clusters = [];
[~, ~, ~, ~, sig_clusters] = OneDimClusterPermutation(-x_pos_diff_trace', 0.01/f_s, 2000);

sig_traces = [];
[~, ~, ~, ~, sig_traces(1, :)] = OneDimClusterPermutation(squeeze(x_pos_dva_e_avgd_filt(1, :, pts))', 0.01/f_s, 2000);
[~, ~, ~, ~, sig_traces(2, :)] = OneDimClusterPermutation(squeeze(-x_pos_dva_e_avgd_filt(2, :, pts))', 0.01/f_s, 2000);

%% Plot 1/3 - x-direction movements
close all;

% Prepare plotting info
t = linspace(ccl(1), ccl(2), size(x_pos_dva_p_avgd, 2));
 
down = squeeze(ci(:, :, 1)); up = squeeze(ci(:, :, 2));

% Plot - eye-averaged position traces post-retrocue across x
figure('units','normalized','outerposition',[0.1 0.1 0.4 0.8]); 
y_val = 0.15;

% Traces
subplot(5, 1, [1, 2, 3]); hold on;
fill([t, fliplr(t)], [-down(1, :), fliplr(-up(1, :))], 'r', 'FaceAlpha', 0.1, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(t, -squeeze(x_pos_dva_p_avgd(1, :)), 'r')

fill([t, fliplr(t)], [-down(2, :), fliplr(-up(2, :))], 'b', 'FaceAlpha', 0.1, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(t, -squeeze(x_pos_dva_p_avgd(2, :)), 'b')

ylabel("Deviation (dva)"); xticks([]); xticklabels([]); legend("Left Cued", "Right Cued"); xlim([-0.5 ccl(2)]); 
ylim([-y_val y_val+0.001]); xline(0, 'k--', 'HandleVisibility', 'off');  yline(0, 'k-', 'HandleVisibility', 'off', 'LineWidth', 2);
plot(t(find(squeeze(sig_traces(1, :)))), -0.05*ones(size(t(find(squeeze(sig_traces(1, :)))))), 'r.', 'linewidth', 2, 'HandleVisibility', 'off');
plot(t(find(squeeze(sig_traces(2, :)))), -0.06*ones(size(t(find(squeeze(sig_traces(2, :)))))), 'b.', 'linewidth', 2, 'HandleVisibility', 'off');

% Difference
subplot(5, 1, [4, 5]); hold on;
fill([t, fliplr(t)], [ci_diff(1, :), fliplr(ci_diff(2, :))], 'k', 'FaceAlpha', 0.1, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(t, squeeze(mean(x_pos_diff_trace, 2)), 'k');
plot(t(find(sig_clusters(:))), 0.02*ones(size(t(find(sig_clusters(:))))), 'g.', 'linewidth', 2, 'HandleVisibility', 'off');

xlabel("Time (s) w.r.t. retro-cue onset"); ylabel("R-L Gaze Bias (dva)"); xlim([-0.5 ccl(2)]); 
y_val = 0.12; ylim([-y_val y_val]); yticks([-0.05, 0, 0.05]);
xline(0, 'k--', 'HandleVisibility', 'off'); yline(0, 'k-', 'HandleVisibility', 'off', 'LineWidth', 2);

sgtitle("Post-cue eye movements, x-direction")

if 0 
    save FigsDataPython/Eye_Ret/t.mat t
    save FigsDataPython/Eye_Ret/ci.mat ci
    save FigsDataPython/Eye_Ret/x_pos_dva_p_avgd.mat x_pos_dva_p_avgd
    save FigsDataPython/Eye_Ret/sig_traces.mat sig_traces
    save FigsDataPython/Eye_Ret/ci_diff.mat ci_diff
    save FigsDataPython/Eye_Ret/x_pos_diff_trace.mat x_pos_diff_trace
    save FigsDataPython/Eye_Ret/sig_clusters.mat sig_clusters
end


%% Participant Average - y

y_pos_dva_p_avgd = squeeze(mean(y_pos_dva_e_avgd_filt(:, :, pts), 3));

% Difference traces
y_pos_diff_trace = squeeze(y_pos_dva_e_avgd_filt(1, :, pts) - y_pos_dva_e_avgd_filt(2, :, pts));
y_pos_diff_trace = squeeze(y_pos_dva_e_avgd_filt(1, :, pts) - y_pos_dva_e_avgd_filt(2, :, pts));

% CIs for time series and differences
ci_y = [];
for att = 1:2
    ci_y(att, :, :) = bootci(5000, @mean, squeeze(y_pos_dva_e_avgd_filt(att, :, pts))')';
end

ci_diff_y = bootci(5000, @mean, squeeze(y_pos_diff_trace'));

% Cluster test
sig_clusters_y = [];
[~, ~, ~, ~, sig_clusters_y] = OneDimClusterPermutation(y_pos_diff_trace', 0.01/f_s, 2000);

sig_traces_y = [];
[~, ~, ~, ~, sig_traces_y(1, :)] = OneDimClusterPermutation(squeeze(y_pos_dva_e_avgd_filt(1, :, pts))', 0.01/f_s, 2000);
[~, ~, ~, ~, sig_traces_y(2, :)] = OneDimClusterPermutation(squeeze(y_pos_dva_e_avgd_filt(2, :, pts))', 0.01/f_s, 2000);

%% Plot 2/3 - y-direction movements

close all

% Prepare plotting info
t = linspace(ccl(1), ccl(2), size(y_pos_dva_p_avgd, 2)); 
down = squeeze(ci_y(:, :, 1)); up = squeeze(ci_y(:, :, 2));

% Plot - eye-averaged position traces post-retrocue across x
figure('units','normalized','outerposition',[0.1 0.1 0.4 0.8]); 
y_val = 0.15;

% Traces
subplot(5, 1, [1, 2, 3]); hold on;
fill([t, fliplr(t)], -[down(1, :), fliplr(up(1, :))], 'r', 'FaceAlpha', 0.1, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(t, -squeeze(y_pos_dva_p_avgd(1, :)), 'r')

fill([t, fliplr(t)], -[down(2, :), fliplr(up(2, :))], 'b', 'FaceAlpha', 0.1, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(t, -squeeze(y_pos_dva_p_avgd(2, :)), 'b')

ylabel("Deviation (dva)"); xticks([]); xticklabels([]); xlim([-0.5 ccl(2)]); 
ylim([-y_val y_val+0.001]); xline(0, 'k--', 'HandleVisibility', 'off');  yline(0, 'k-', 'HandleVisibility', 'off', 'LineWidth', 2);
plot(t(find(squeeze(sig_traces_y(1, :)))), 0.05*ones(size(t(find(squeeze(sig_traces_y(1, :)))))), 'r.', 'linewidth', 2, 'HandleVisibility', 'off');
plot(t(find(squeeze(sig_traces_y(2, :)))), 0.06*ones(size(t(find(squeeze(sig_traces_y(2, :)))))), 'b.', 'linewidth', 2, 'HandleVisibility', 'off');

% Difference
subplot(5, 1, [4, 5]); hold on;
fill([t, fliplr(t)], [ci_diff_y(1, :), fliplr(ci_diff_y(2, :))], 'k', 'FaceAlpha', 0.1, 'HandleVisibility', 'off', 'EdgeColor', 'none');
plot(t, squeeze(mean(y_pos_diff_trace, 2)), 'k');
plot(t(find(sig_clusters_y(:))), -0.05*ones(size(t(find(sig_clusters_y(:))))), 'g.', 'linewidth', 2, 'HandleVisibility', 'off');

xlabel("Time (s) w.r.t. retro-cue onset"); ylabel("Difference (dva)"); xlim([-0.5 ccl(2)]); 
y_val = 0.1; ylim([-y_val y_val]); yticks([-0.05, 0, 0.05]);
xline(0, 'k--', 'HandleVisibility', 'off'); yline(0, 'k-', 'HandleVisibility', 'off', 'LineWidth', 2);

sgtitle("Post-cue eye movements, y-direction")


%% Plot 3/3 - number of trials spread
close all

figure('units','normalized','outerposition',[0.5 0 0.12 0.75])

trial_inclusion_rates = 1 - sum(avoid_trials, 1)/size(avoid_trials, 1);
x = ones(1, length(participants)) + 0.02*rand(1, length(participants));
c = linspace(1, 7, length(x));

scatter(x, trial_inclusion_rates, [], c, "filled"); 

xlim([0.98 1.04]); ylim([0 1]); xticks([]); ylabel("Proportion of trials included")
function [interpolated_data, blink_timepoints] = interpolate_blinks(pupil_data, data_to_correct, f_s)
% % This function uses pupil velocity data to identify periods corresponding to blinks, and interpolates
%   a given feature of the eyetracking data with a pchip interpolation over those periods.  
% % Inputs:
%   pupil_data - 1xn array, pupil timeseries
%   data_to_correct - 1xn array, must be the same dimension as pupil_data (could also be pupil_data itself, or xposition/yposition/etc)
%   f_s - int, sampling frequency in Hz
% % Outputs:
%   interpolated_data - 1xn array, corrected for blinks
%   blink_timpoints - n_blink x 2 array, identified onset and offset times for each blink
% % Notes:
%   - Authored by Kabir Arora, October 2023
%   - Based on Hershman, R., Henik, A. & Cohen, N. A novel blink detection method based on pupillometry noise. Behav Res 50, 107–114 (2018). https://doi.org/10.3758/s13428-017-1008-1
%   - Threshold choices are arbitrary and may depend on use purposes, amend below:

% Parameters
filter_width_time  = 10; % in ms, window for smoothing
max_velocity = 10; max_acceleration = 2; % in units/sec or /sec^2, used to differentiate "normal" data from blink-time data

%% Identify blank periods

% Case 1: both starting and ending with a blink
if (pupil_data(1) == 0) && (pupil_data(end) == 0)
    nan_period_endpoints = zeros(length(find(diff(pupil_data==0) == 1))+1, 2);

    nan_period_endpoints(1, 1) = 1;

    nan_period_endpoints(2:end, 1) = find(diff(pupil_data==0) == 1);
    nan_period_endpoints(1:end-1, 2) = find(diff(pupil_data==0) == -1);

    nan_period_endpoints(end, 2) = length(pupil_data);

% Case 2: only starting with a blink
elseif pupil_data(1) == 0
    nan_period_endpoints = zeros(length(find(diff(pupil_data==0) == 1))+1, 2);

    nan_period_endpoints(1, 1) = 1;

    nan_period_endpoints(2:end, 1) = find(diff(pupil_data==0) == 1);
    nan_period_endpoints(:, 2) = find(diff(pupil_data==0) == -1);

% Case 3: only ending with a blink        
elseif pupil_data(end) == 0
    nan_period_endpoints = zeros(length(find(diff(pupil_data==0) == 1)), 2);

    nan_period_endpoints(:, 1) = find(diff(pupil_data==0) == 1);
    nan_period_endpoints(1:end-1, 2) = find(diff(pupil_data==0) == -1);

    nan_period_endpoints(end, 2) = length(pupil_data);

% Case 4: neither end is a blink
else
    nan_period_endpoints = zeros(length(find(diff(pupil_data==0) == 1)), 2);

    nan_period_endpoints(:, 1) = find(diff(pupil_data==0) == 1);
    nan_period_endpoints(:, 2) = find(diff(pupil_data==0) == -1);

end

% first measured timepoint instead of last NaN timepoints
nan_period_endpoints(:, 2) = nan_period_endpoints(:, 2) + 1;

%% Smooth data 

sampling_interval = round(1000/f_s);
filter_width = ceil(filter_width_time/sampling_interval); % amount of samples to smooth 
smooth_data    = smooth(pupil_data, filter_width);    

smooth_data(smooth_data==0) = NaN;
diff_smooth_data = diff(smooth_data); % functionally the derivative

%% Improve estimate of blink periods

for blink = 3:size(nan_period_endpoints, 1)-1

    % Starting point
    onset_candidate = nan_period_endpoints(blink, 1);
    offset_candidate = nan_period_endpoints(blink, 2);

    % Relevant portions
    pre_derivative = diff_smooth_data(1:onset_candidate);
    post_derivative = diff_smooth_data(offset_candidate+1:end);

    % Second derivative
    pre_diff_2 = diff(pre_derivative);
    post_diff_2 = diff(post_derivative);

    % Find two closest points on the derivative with
    %   - magnitude <10 units
    %   - difference <2 units
    
    final_onset = 0; idx = 0;
    while ~final_onset
        if (pre_derivative(end-idx) < max_velocity) && (pre_derivative(end-idx-1) < max_velocity) && (pre_diff_2(end-idx) < max_acceleration)
            final_onset = 1; onset_actual = idx;
        end
        idx = idx + 1;
    end

    final_offset = 0; idx = 1;
    while (~final_offset) && (idx < length(post_derivative))
        if (post_derivative(idx) < max_velocity) && (post_derivative(idx+1) < max_velocity) && (post_diff_2(idx) < max_acceleration)
            final_offset = 1; offset_actual = idx;
        end
        idx = idx + 1;
    end    

    actual_endpoints(blink, 1) = onset_candidate - idx;
    actual_endpoints(blink, 2) = offset_candidate + idx - 1;
end


%% Use blink periods to interpolate required data

% % Option 1/2, interpolate across blink
% blink_period = false(1, length(pupil_data));
% 
% for blink = 2:size(nan_period_endpoints, 1)-1
%     blink_period(actual_endpoints(blink, 1):actual_endpoints(blink, 2)) = true;
% end
% 
% data_corrected = data_to_correct;
% data_corrected(blink_period) = NaN;
% 
% % Interpolate
% interpolated_data = fillmissing(data_corrected, 'pchip');

% Option 2/2, replace blink edges with NaNs as well
data_corrected = data_to_correct;
for blink = 2:size(nan_period_endpoints, 1)-1
    data_corrected(max(actual_endpoints(blink, 1), 1):max(actual_endpoints(blink, 2), 1)) = NaN;
end
interpolated_data = data_corrected;

%% Test it out over two examples
% close all
% 
% n_int = 2; % two intervals
% interval{1} = 1000:1300; interval{2} = 7075:7200;
% blink_idxs = [6 20];
% 
% for int = 1:n_int
%     figure;
% 
%     subplot(1, 2, 1); hold on;
%     plot(interval{int}, interpolated_data(interval{int}), 'DisplayName', "Raw Data");
%     scatter(nan_period_endpoints(blink_idxs(int), :), interpolated_data(nan_period_endpoints(blink_idxs(int), :)), 400, 'r.', 'DisplayName', "NaN Period")
%     scatter(actual_endpoints(blink_idxs(int), :), interpolated_data(actual_endpoints(blink_idxs(int), :)), 400, 'b.', 'DisplayName', "Blink Period");
%     ylabel("Pupil Size"); title("Blink Period as per Algorithm"); legend(); xlim([interval{int}(1), interval{int}(end)]);
% 
%     subplot(1, 2, 2); hold on;
%     plot(interval{int}, interpolated_data(interval{int}), 'DisplayName', "Interpolated Data");    
%     plot(interval{int}, interpolated_data(interval{int}), 'DisplayName', "Raw Data");
%     title("Interpolation"); legend(); xlim([interval{int}(1), interval{int}(end)]);
% end
end

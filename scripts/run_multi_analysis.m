clear all;clc
sbj_group = [2:10 13:16]; % modify make_trial for sbj 5
for isbj = 1:numel(sbj_group)
    selected_sbj = sbj_group(isbj);
    RIFT02_spectral(selected_sbj)
    % RIFT01_Eye_tracker(selected_sbj)
end

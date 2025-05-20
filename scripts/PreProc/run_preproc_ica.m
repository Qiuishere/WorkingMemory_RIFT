function run_preproc_ica(subj_id)

addpath('/project/3018085.01/scripts/subfun')
subjects = Start_up;


fprintf('*** SUBJECT %02d : ica ***\n', subj_id);

try
% load data
load(fullfile(subjects(subj_id).dir, 'preproc-data-artreject-400hz.mat'),...
  'data');


%% run ICA
cfg = [];
cfg.method = 'runica';
cfg.demean = 'no';
cfg.channel = {'MEG'}; % only do ICA on MEG channels, not the refchans

comp = ft_componentanalysis(cfg, data);

unmixing = comp.unmixing;
topolabel = comp.topolabel;

save(fullfile(subjects(subj_id).dir, 'preproc-ica-weights.mat'),...
  'unmixing', 'topolabel');

catch ME
    save(fullfile(subjects(subj_id).dir, 'errorMessage.mat'))
end
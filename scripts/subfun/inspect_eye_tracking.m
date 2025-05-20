% Fieldtrip tutorial https://www.fieldtriptoolbox.org/getting_started/eyelink/
addpath('/home/predatt/sonbai/RIFT/subfun/')
addpath('/home/predatt/sonbai/Documents/MATLAB/fieldtrp/')
addpath('/home/predatt/sonbai/RIFT/Analysis Script')
ft_defaults
ft_warning off

filename_eye = '/home/predatt/sonbai/RIFT/behavioural/02/edf/s11.asc';
cfg = [];
cfg.dataset          = filename_eye;
cfg.montage.tra      = eye(4);
cfg.montage.labelorg = {'1', '2', '3', '4'};
cfg.montage.labelnew = {'EYE_TIMESTAMP', 'EYE_HORIZONTAL', 'EYE_VERTICAL', 'EYE_DIAMETER'};
data_eye_general = ft_preprocessing(cfg);

event_eye = ft_read_event(filename_eye);

% remember to remove the diameter channel and change the scale, the
% platform at EYE_horizontal/vertical should overlap with FIX
cfg = [];
cfg.viewmode       = 'vertical';
cfg.preproc.demean = 'yes';
cfg.event          = event_eye;
ft_databrowser(cfg, data_eye_general);


%% read the EEG data, this could also come from an EEG dataset
cd('/home/predatt/sonbai/RIFT/PreProc')
subjects = datainfo();

subj_id = 2;
fprintf('*** SUBJECT %02d : epoching ***\n', subj_id);

cfg = [];
cfg.dataset = subjects(subj_id).rawmeg;
cfg.trialdef.eventtype      = 'UPPT001';
cfg.channel = {'UADC*'}; % read only the eye-tracking moch channels
cfg.trialdef.eventvalue     = 11;
cfg.trialdef.prestim        = 5.5;
cfg.trialdef.poststim       = .5;
cfg.continuous  = 'yes';
cfg = ft_definetrial(cfg);
data_meg = ft_preprocessing(cfg);

% read in eye tracking data
cfg = [];
cfg.dataset = '/home/predatt/sonbai/RIFT/behavioural/02/edf/s11.asc';
cfg.trialdef.eventtype      = 'INPUT';
cfg.trialdef.eventvalue = 11; % stimuli offset
cfg.trialdef.prestim = 5.5; % 5.5;
cfg.trialdef.poststim = .5; %.5;
cfg = ft_definetrial(cfg);
data_eye_general = ft_preprocessing(cfg);



% plot the trigger from MEG and eye tracker
uadc005 = find(strcmp(data_meg.label, 'UADC005'));
uadc006 = find(strcmp(data_meg.label, 'UADC006'));


itrial = 110;
figure
subplot(3,1,1)
plot(data_eye_general.time{itrial}, data_eye_general.trial{itrial}(2,:))
grid on

subplot(3,1,2)
plot(data_meg.time{itrial}, data_meg.trial{itrial}(uadc005,:))
grid on

%% resample MEG to 1000Hz and combine with eye-tracking
cfg = [];
cfg.time = data_meg.time;
data_eye_resampled = ft_resampledata(cfg, data_eye);

cfg = [];
data_combined = ft_appenddata(cfg, data_meg, data_eye_resampled);



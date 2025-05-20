%% NOTE: this is supposed to be run as a script, not a function
addpath('/home/predatt/qiuhan/Documents/Toolboxs/fieldtrip-20250114')
addpath('/project/3018085.01/scripts/subfun')
ft_defaults
ft_warning off


subjects = datainfo();

subj_id = 6;

fprintf('*** SUBJECT %02d : loading data... ***\n', subj_id);

load(fullfile(subjects(subj_id).dir, 'preproc-data-artreject-400hz.mat'),...
  'data');
load(fullfile(subjects(subj_id).dir, 'preproc-ica-weights.mat'),...
  'unmixing', 'topolabel');


cfg = [];

cfg.method = 'predefined mixing matrix';
cfg.demean = 'no';
cfg.channel = {'MEG'};
cfg.topolabel = topolabel;
cfg.unmixing = unmixing;
comp = ft_componentanalysis(cfg, data);

cfg = [];
cfg.viewmode = 'component';
cfg.layout = 'CTF151.lay';
ft_databrowser(cfg, comp);





%% Alternatively, plot the components for visual inspection
cfg           = [];
cfg.layout    = 'CTF151.lay';
cfg.component = [46, 69];
cfg.marker    = 'off';
ft_topoplotIC(cfg, comp)



% plot ECGs together to see if they synchronize
figure
cfg = [];
cfg.channel = {'component010','component050'};       % specify the component(s) that should be plotted
cfg.layout    = 'CTF151.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_databrowser(cfg, comp)


%% write down and save
fprintf('*** SUBJECT %02d : save the identified components!!! ***\n', subj_id);

badcomps = [1, 10, 16, 50] ;
badcomps_reasons = { 'heart beat','eye', 'heart beat','eye',  }; %{'heart beat','heart beat','eye'};


assert(numel(badcomps) == numel(badcomps_reasons));

save(fullfile(subjects(subj_id).dir, 'preproc-ica-badcomps.mat'),...
  'badcomps', 'badcomps_reasons');

fprintf('saved.\n');

% just save it, and reject the components when needed using
% load_clean_data()
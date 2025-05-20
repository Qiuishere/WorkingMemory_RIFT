function [cfg, topotemplate] = plot_topotemplate(data,temp,chans)
load('topotemplate.mat')
topotemplate.powspctrm = data;

cfg = [];
if strcmp('custom',temp)
    cfg.layout = 'customlay.mat';
else
    cfg.layout = 'biosemi64.lay';
end
if exist('chans','var')
    cfg.channel = chans;
    topotemplate.label = chans;
end
cfg.comment = 'no';
cfg.interactive = 'no';
end
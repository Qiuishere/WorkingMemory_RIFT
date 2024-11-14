addpath('.\subfun')
init 

RIFT_order = readmatrix('tag_position.txt');
% for testing
% Screen('Preference', 'SkipSyncTests', 1)
% Screen('Preference','VisualDebugLevel', 0)
%%=========================================================================
% get main parameters and subject code
commandwindow;
% 
prm            = RIFT_prm('practice'); 
prm.exp.ntrials  = 4;  % run only the first 8 trials
prm.exp.blocks   = 2;
prm.exp.design(prm.exp.ntrials + 1 :end,:)  = [];
prm.exp.design{3:end,5}  = 2; % BLOCK   
prm.exp.design{3:end,9}  = 1; % FIXMODE


% change the RIFT position based on sbj id
id             = subject_id(prm, []); 

%% balance the RIFT freq across 
if str2double(id.number) <= 40
    order = RIFT_order(str2double(id.number),:);
    temp = prm.stim.tag_frex;
    prm.stim.tag_frex(order) = temp;
end
%%=========================================================================
% run the practice
RIFT_main(prm,id)
% %%=========================================================================


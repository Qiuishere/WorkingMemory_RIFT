function SbjResp = MemoryTask(data,block)
%%
% This function should show several iamge where subject can give an 
% response to say whether they observe the same object in the past blocks
% To do:
% show the stimuli
% report whether yes or no
%%
default('location',[]);
prm                = data.prm;
sub_result         = tbl_subset(data.results,'block',block);
ntrial             = size(sub_result,1);
AllShownImg        = [];
for i = 1:ntrial
    AllShownImg        = [AllShownImg sub_result(i,:).stim_idx{1}];
end
ShownImgIdx     = unique(AllShownImg);  
UnShownImgIdx      = setdiff(1:numel(prm.stim.AllStimPath),ShownImgIdx);

if numel(UnShownImgIdx) < 2
    prm.task.nShown = 4;
    prm.task.nUnShown = 0;
end


assert(numel(ShownImgIdx)+numel(UnShownImgIdx) == numel(prm.stim.AllStimPath), ...
    'Number of images incorrect');

SelectedImgIdx   = [ShownImgIdx(randperm(numel(ShownImgIdx),prm.task.nShown)) ...
                 UnShownImgIdx(randperm(numel(UnShownImgIdx),prm.task.nUnShown))];
GroudTruth    = [ones(1,prm.task.nShown) zeros(1,prm.task.nShown)];
ShuffleIdx    = randperm(4); % shuffle the order of stimuli

SelectedImg   = prm.stim.AllStimPath(SelectedImgIdx);
SelectedImg   = SelectedImg(ShuffleIdx);
GroudTruth    = GroudTruth(ShuffleIdx);
Rt            = zeros(1,numel(GroudTruth));
SubjectAns    = zeros(1,numel(GroudTruth));

%
for iStim = 1:numel(SelectedImg)  
   [img, ~, alpha]  = imread(SelectedImg(iStim)); 
    sti(:, :, 1:3)   = repmat(img,[1 1 3]);
    sti(:, :, 4)     = alpha;
    Texture = Screen('MakeTexture', prm.monitor.window, sti);
    Screen('DrawTexture', prm.monitor.window, Texture, [], []);
    DrawFormattedText(prm.monitor.window, 'Press LEFT if you have seen this object, otherwise press RIGHT', ...
               'center',3*prm.monitor.height/4, prm.monitor.white, [], [], [], 1.5);
    
%     fprintf('clearing responses\n');
    prm.trigger.btsi.clearResponses();
    Screen('Flip',prm.monitor.window);
    prm.trigger.btsi.sendTrigger(prm.trigger.MemorySti); % onset of memory stimuli
    % get subject response (Inf = no timeout; true = return on button press)
    %fprintf('wait for response\n');
    [this_resp,this_rt] = prm.trigger.btsi.getResponse(Inf, true);

    prm.trigger.btsi.sendTrigger(prm.trigger.response);
    Rt(iStim) = this_rt;
    SubjectAns(iStim) = this_resp;

    Screen('Flip',prm.monitor.window);

    WaitSecs(1); 
end

SbjResp.SelectedImg   = SelectedImg';
SbjResp.GroudTruth    = GroudTruth;
SbjResp.rt            = Rt;
SbjResp.SubjectAns    = SubjectAns;
Screen('Flip',prm.monitor.window);   
% Screen('Close', response_texture) 

end
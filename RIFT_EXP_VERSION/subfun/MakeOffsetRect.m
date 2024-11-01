function rect = MakeOffsetRect(ptb, stim, dx, dy, quadrant)
% modified form Eelke's original function: 
% Step 1: position of each stimuli is first assigned according to
% StiMatrix [sub_w sub_h] and quardrant. Now one of its corner is at the
% center of the screen.
% Step 2: move the stimuli based on its size [stim_w, stim_h] so the
% position of center are in line with PTB format
% Step 3: add an additional offset [dx dy] if more than one stimuli are to
% be draw within the quardrant
% ptb: a struct of monitor information from InitPsychtoolbox.m
% stim: para of stimuli
% dx, dy: x and y offset from the center of the screen
% quardrant: [1 4]
% Nstimli：how many stimuli will be draw in each quardrant, a matrix of 
% [Nrow * Ncolumn]
%                                                   Songyun  15.01.2024 DCC
stim_w = size(stim.grating, 1);
stim_h = size(stim.grating, 2);


sub_w = ptb.width/4;
sub_h = ptb.height/4;

center = [ptb.width/2, ptb.height/2, ptb.width/2, ptb.height/2];
if nargin < 4
    % origin is center of actual screenbuffer
    rect = [-stim_w/2+dx, ptb.height/2-stim_h/2+dy, ptb.width/2+stim_w/2+dx, +stim_h/2+dy];
   
    
elseif quadrant == 1
    % top left quadrant
    rect = center + [-sub_w-stim_w/2+dx, -sub_h-stim_h/2+dy, -sub_w+stim_w/2+dx, -sub_h+stim_h/2+dy];
    
    % check that we're not attempting to draw outside the selected quadrant
%     assert(rect(1) >= 0 && rect(2) >= 0 && rect(3) <= ptb.width/2 && rect(4) <= ptb.height/2);
    
elseif quadrant == 2
    % top right quadrant
    rect = center + [sub_w-stim_w/2+dx, -sub_h-stim_h/2+dy, sub_w+stim_w/2+dx, -sub_h+stim_h/2+dy];
    
    % check that we're not attempting to draw outside the selected quadrant
%     assert(rect(1) >= ptb.width/2 && rect(2) >= 0 && rect(3) <= ptb.width && rect(4) <= ptb.height/2);
    
elseif quadrant == 3
    % bottom left quadrant
    rect = center + [-sub_w-stim_w/2+dx, sub_h-stim_h/2+dy, -sub_w+stim_w/2+dx, sub_h+stim_h/2+dy];
    
    % check that we're not attempting to draw outside the selected quadrant
%     assert(rect(1) >= 0 && rect(2) >= ptb.height/2 && rect(3) <=  ptb.width/2 && rect(4) <= ptb.height);
    
elseif quadrant == 4
    % bottom right quadrant
    rect = center + [sub_w-stim_w/2+dx, sub_h-stim_h/2+dy, sub_w+stim_w/2+dx, sub_h+stim_h/2+dy];
    
    % check that we're not attempting to draw outside the selected quadrant
%     assert(rect(1) >= ptb.width/2 && rect(2) >= ptb.height/2 && rect(3) <=  ptb.width && rect(4) <= ptb.height);
    
else
    error('invalid quadrant specification');
end

end
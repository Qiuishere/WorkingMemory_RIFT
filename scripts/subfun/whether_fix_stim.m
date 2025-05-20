% Songyun 31.05.2024 at DCCN
% check whehther the fixation is within the region of stimuli
% fix_x, fix_y: represent fixation coordinate
% Rect: 4 by 4 matrix contaning [quardrant, location of stim]
% Each row is [X_leftbottom, Y_leftbottom, width, height]
%
% CAUTION: rectangle set the zero point at the left-lower corner, 
% while PTB set it at left upper corner

function [pos_idx] = whether_fix_stim(fix_x, fix_y, Rect)

    pos_idx = nan;
    for iquadrant = 1:4  
        Inrange_X = Rect(iquadrant,1) < fix_x && fix_x < (Rect(iquadrant,1) + Rect(iquadrant,3));
        Inrange_y = Rect(iquadrant,2) < fix_y && fix_y < (Rect(iquadrant,2) + Rect(iquadrant,4));
        if Inrange_X && Inrange_y
            pos_idx = iquadrant;
        end   
    end
end
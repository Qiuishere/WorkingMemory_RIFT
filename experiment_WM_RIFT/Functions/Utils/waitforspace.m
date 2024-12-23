function [] = waitforspace(prm)

if IsWin
    spacebar = 32;  esckey  =   27;
elseif IsOSX
    spacebar = 44;  esckey  =   41;
end

while 1
    
    if IsOSX
        [a,~,c] = KbCheck(-1);
    else
        [a,~,c] = KbCheck;
    end
    if a && ismember(spacebar,find(c))
        break
    elseif a && ismember(esckey,find(c))
        sca; ShowCursor;
        error('[!!!] Program aborted by user');
    end
    
    
    % for buttonbox, also check for space
    if prm.exp.RealRun ==1
        
        [keydown,~] = prm.trigger.btsi.waitResponse(0.001, true);
        
        if keydown==prm.exp.key.space
            break
        end
    end
    
end
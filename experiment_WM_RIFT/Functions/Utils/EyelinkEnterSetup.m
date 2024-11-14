function EyelinkEnterSetup(el)
% EyelinkEnterSetup(el)

fprintf('\n\nCamera setup started (%.4f)...\n', GetSecs());

elsetupresult = true;
while elsetupresult
    % dummy connection indicated by 'el.dummyconnected == -1'
    if Eyelink('IsConnected') > 0
        elsetupresult = EyelinkDoTrackerSetup(el);
        if elsetupresult
            fprintf('EyelinkDoTrackerSetup failed. Trying again...\n')
            fprintf('There are now 3 seconds time to press escape to quit the experiment until setup is tried again...\n');
            [~, keycode] = KbWait([], 2, GetSecs() + 3);
            if any(find(keycode) == KbName('Escape'))
                error('Aborted EL setup. Experiment terminated.');
            end
        end
    elseif Eyelink('IsConnected') == 0 % ---------- no tracker connection ----------
        error('Tracker not connected!');
    elseif Eyelink('IsConnected') == -1 % --------- no setup, since DUMMYMODE --------
        fprintf('\n**********************\nNo Eyelink Setup, since running in DUMMYMODE = TRUE!\n**********************\n')
        elsetupresult = 0;
    end % no ELSE, since no other possibility
end
fprintf('Camera setup done (%.4f).\n', GetSecs());

return
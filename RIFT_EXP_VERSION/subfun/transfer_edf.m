function transfer_edf(prm)
%function transfer_edf(G)
%transfers the edf file
%try
    fprintf('Receiving data file ''%s''\n', prm.exp.eyeFile);
    status = Eyelink('ReceiveFile');
    fprintf('ReceiveFile status %d\n', status);
    WaitSecs(1.0); % Give Tracker time to execute all commands
    % move EDF file to subfolder / clean-up
    fprintf('Moving edf file to subfolder ''%s''...\n', prm.exp.eyeDir);
    %% FIGURE OUT IF THIS IS CORRECT
    [mvEDFsuccess, mvEDFmessage, mvEDFmessID] = movefile(prm.exp.eyeFile,prm.exp.eyeDir); % check if the eyeDir here makes sense
    if ~mvEDFsuccess
        warning(mvEDFmessID, mvEDFmessage);
    end
    Eyelink('Command', 'clear_screen 0'); % Clear trial image on Host PC at the end of the experiment
    WaitSecs(0.1); % Allow some time for screen drawing
%catch
%    fprintf('Failed receiving the datafile.\n');
%end
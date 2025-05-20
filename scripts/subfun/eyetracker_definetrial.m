function [trl, event, fix] = eyetracker_definetrial(cfg)
    % read the header information and the events from the data
    hdr   = ft_read_header(cfg.dataset);
    hdr.FirstTimeStamp = cfg.alltimestamp(1);
    All_Sample = 1:numel(cfg.alltimestamp);
    event = ft_read_event(cfg.dataset);




    % search for "trigger" events
    value_trg  = [event(find(strcmp(cfg.trialdef.eventtype, {event.type}))).value]';
    timestamp_trg = [event(find(strcmp(cfg.trialdef.eventtype, {event.type}))).timestamp]';
    sample_trg = [event(find(strcmp(cfg.trialdef.eventtype, {event.type}))).sample]';

  
    % since the eye-tracking data is not continuous, the sample should be
    % reset
    for i = 1:numel(value_trg)
        selcted_timestamp = timestamp_trg(i);
        try
            sample_trg(i) = All_Sample(cfg.alltimestamp==selcted_timestamp);
        catch 
            sample_trg(i) = nan;
        end
        
    end


    % figure
    % subplot(211)
    % plot(sample/hdr.Fs, value, '.');
    % subplot(212)
    % plot(timestamp/hdr.Fs, value, '.');


    
    % determine the number of samples before and after the trigger
    pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
    posttrig =  round(cfg.trialdef.poststim * hdr.Fs);
    
    % look for trigger
    trl = [];
    for j = 1:(length(value_trg))
        trg1 = value_trg(j);
        if trg1==cfg.trialdef.eventvalue
          trlbegin = sample_trg(j) + pretrig;
          trlend   = sample_trg(j) + posttrig - 1;
          offset   = pretrig;
          newtrl   = [trlbegin trlend offset];
          trl      = [trl; newtrl];
        end
    end
end

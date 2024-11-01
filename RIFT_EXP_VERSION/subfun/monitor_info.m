function monitor = monitor_info(screen_number,view_distance,screen_size)

default('view_distance',60)
default('screen_number',max(Screen('Screens')))
default('screen_size',[])
% initialize
PsychDefaultSetup(2);
% screen 
monitor            = Screen('resolution',  screen_number);
% if ismac
%     hz             = fix(FrameRate);
%     monitor        = Screen('resolution',screen_number,[],[],hz,32);
% end
% if ismac
%     monitor.hz = fix(FrameRate);
%     temp            = Screen('resolution',  screen_number);
%     set(0,'units','pixels')  
%     %Obtains this pixel information
%     Pix_SS = get(0,'screensize');
%     monitor.screen_pix = Pix_SS;
%     monitor.width = Pix_SS(3);
%     monitor.height = Pix_SS(4);
%     monitor.pixelSize = temp.pixelSize;
%     screen_number = max(Screen('Screens'));
% end
monitor.view_dist  = view_distance;
monitor.dt         = 1/monitor.hz;
monitor.frame_ms   = @(ms) round(ms/(1000/monitor.hz));
[width, height]    = Screen('DisplaySize', screen_number);

monitor.screen_cm  = .1*[width, height];
monitor.screen_pix = [0 0 monitor.width monitor.height];
[x, y]             = RectCenter(monitor.screen_pix);
monitor.center     = [x, y];
monitor.screen_num = screen_number;
% white, black and gray     
monitor.black      = BlackIndex(screen_number);
monitor.white      = WhiteIndex(screen_number);
monitor.grey       = monitor.white / 2;
% pixels, cm and degrees
monitor.pix_in_cm  = monitor.screen_cm(1)/(monitor.width/2);
monitor.deg_to_cm  = 2*view_distance*tan(pi*1/(2*180));
monitor.pix_to_deg = monitor.screen_cm(1)/(monitor.width/2)/ ...
                     view_distance*180/pi;
% monitor.deg_to_pix = fix(monitor.deg_to_cm/monitor.pix_in_cm);
monitor.deg_to_pix = fix(1/monitor.pix_to_deg);
% additional
info.struct        = Screen('Version');
info.comp          = Screen('Computer');
info.arch          = computer;
monitor.info       = info;
% =========================================================================
% update fields when a different screen size is required
if ~isempty(screen_size)
    original           = monitor;
    % overwrite
    [monitor.width, monitor.height] = deal(screen_size(3),screen_size(4));
    monitor.original   = original;
    monitor.rate       = 1/monitor.hz;
    monitor.frame_ms   = @(ms) round(ms/(1000/monitor.hz));
    ratio_to_original  = screen_size(3:4)./original.screen_pix(3:4);
    monitor.screen_cm  = original.screen_cm.*ratio_to_original;    
    monitor.screen_pix = [0 0 monitor.width monitor.height];
    [x, y]             = RectCenter(monitor.screen_pix);
    monitor.center     = [x, y];
    monitor.screen_num = monitor.original.screen_num;
    % white, black and gray     
    monitor.black      = BlackIndex(monitor.screen_num);
    monitor.white      = WhiteIndex(monitor.screen_num);
    monitor.grey       = monitor.white / 2;
    % pixels, cm and degrees
    monitor.pix_in_cm  = monitor.screen_cm(1)/(monitor.width/2);
    monitor.deg_to_cm  = 2*monitor.original.view_dist*tan(pi*1/(2*180)); % only determined by viewing distance
    monitor.deg_to_pix = fix(monitor.deg_to_cm/monitor.pix_in_cm);
end
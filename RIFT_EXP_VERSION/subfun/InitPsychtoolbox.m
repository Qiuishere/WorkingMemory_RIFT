function ptb = InitPsychtoolbox(ptb, is_live)
% InitPsychtoolbox performs some basic setup calls for PsychToolbox, as
% well gets the indices of the screen's white colour, screen dimensions,
% etc.

% initialize a structure in which to store useful information

ptb.is_live = is_live;

% skip the ptb synchronization tests when not running actual experiment
if ~is_live
    Screen('Preference', 'SkipSyncTests', 1);
    PsychDebugWindowConfiguration([], 1); % change second argument to reduce screen opacity
    commandwindow();
else
    Screen('Preference', 'SkipSyncTests', 0);
    PsychDebugWindowConfiguration([], 1);
    commandwindow();
%     HideCursor();
end
HideCursor();
% some default setup
AssertOpenGL();

% ensure consistent mapping of keyboard buttons to character labels
KbName('UnifyKeyNames');

% get available screens and draw onto the one with highest ID (usually an
% extra screen on which Matlab is not running)
screens = Screen('Screens');

ptb.screen_id = max(screens);

% get pointer to window on which to draw
[ptb.window, ptb.win_rect] = Screen('OpenWindow', ptb.screen_id,...
    127);%, [0 0 1920 1080]);%, [], [], [], [], 4);

% query the frame duration
ptb.ifi = Screen('GetFlipInterval', ptb.window);

% set the text font and size
Screen('TextFont', ptb.window, 'Calibri');
Screen('TextSize', ptb.window, 18);
Screen('TextStyle', ptb.window, 0);

% get the size of the on screen window
[ptb.width, ptb.height] = Screen('WindowSize', ptb.window);

% set the blend function for the screen
Screen('BlendFunction', ptb.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% initialize the audio device
% InitializePsychSound();
% ptb.snd = PsychPortAudio('Open', [], [], 0, 48000);

% set script to high priority (1 on Windows)
Priority(1);

% initialize GetSecs (to ensure the mexfile is loaded)
GetSecs();
WaitSecs(0.001);

% clear the screen
Screen('Flip', ptb.window);

end
% GENERAL SETTINGS for equipments: keyboard, monitor, log file

%% Clean up


AssertOpenGL;

%% Set random seed

[seed, whichGen] = ClockRandSeed;


%% Screen preferences

% Suppress non-critical tests and warnings      (TEST ELSEWHERE!)
Screen('Preference', 'VisualDebuglevel',    3);
Screen('Preference', 'SuppressAllWarnings', 1);
Screen('Preference',  'SkipSyncTests',      2);


%% Keyboard

KbName('UnifyKeyNames')
% participant keys
if RealRun
    Keys.name = {'Index', 'Middle'};
    Keys.number = [97, 98];
else
    Keys.name = {'UpArrow', 'DownArrow'};
    Keys.number = [KbName('UpArrow') KbName('DownArrow')];
end



RespKeys = Keys.number;
% restrict keys that are read out by KbCheck (speeds up KbCheck)
enabledKeys = [RespKeys KbName('ESCAPE') KbName('SPACE')];
RestrictKeysForKbCheck(enabledKeys);

%% Open screen & settings==================================================
sca;
AssertOpenGL;

text.Size = 24;
text.Color = prm.monitor.white;
text.Font = 'Arial';


prm.w = SetScreen('Window', prm.monitor.screenId,'OpenGL', 1, 'BGColor', 255* prm.monitor.grey, 'Blending', 'Transparent',...
    'Debug', 0,'FontSize',text.Size,'FontColor', text.Color, 'FontType',text.Font);

topPriorityLevel = MaxPriority(prm.w.Number);
Priority(topPriorityLevel);

if RealRun
    HideCursor(prm.monitor.screenId);
end


%% Make directories for participant (if needed)============================
if RealRun==1
    RunDir = fullfile(prm.exp.datedir, sprintf('Sub%02d', SubNo), RunType);
else
    RunDir = fullfile(prm.exp.datedir, sprintf('Sub%02d', SubNo), RunType);
end

if ~exist(RunDir, 'dir')
    mkdir(RunDir);
end


%% Set diary and record the current script

LogDir = fullfile(RunDir, 'Logs');
if ~exist(LogDir, 'dir')
    mkdir(LogDir);
end

RunStart = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
diary(fullfile(LogDir, sprintf('Sub%02d_%s_%g_%s.txt', SubNo, RunType, RunStart)));


ScriptName = {'Adjustment_withMask_RIFT','Settings_General_MEG', 'RUN'};
myCode = cell(length(ScriptName),1);
for i = 1:length(ScriptName)
    fid = fopen([ScriptName{i} '.m'],'r'); % 'r': open file for reading
    myCode{i} = fscanf(fid,'%300c');
    fclose('all');
end
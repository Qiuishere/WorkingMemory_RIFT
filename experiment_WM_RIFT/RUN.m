%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% MEG study: probing working memory biases with RIFT                      %
%                                                                         %
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                                                                         %
%    Qiu Han                          %
%    Partly based on Songyun Bai's code                                   %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
sca; clear all; close all; clc; commandwindow;
addpath(genpath('Functions'));
   
SubNo = 4;
RunType = 'pra';

global Environment
Environment = 1; % 1 - office, 2 - MEG lab

global RealRun
RealRun = 0;

if Environment==2 && RealRun ==1
    useEye = 1;
else useEye = 0;
end

prm = prm_RIFT(RunType, RealRun, Environment, useEye);


%% define experiment Runs

% experiment information
Runs = {'pra', ...
        'test'};

% get the current runno from data folder, as the default input
[RunNo, NOccur] = EstiRun(SubNo, Runs, prm.exp.datedir); % NOccur returns the index 
% of each run inside their own category (e.g. the second training run)

% enquire about current run
prompt      =   {'Enter next run number:',...
    'Enter Participant number:',...
    };
title       =   sprintf('Participant %g', SubNo);
defInput    =   {num2str( RunNo), num2str(SubNo), num2str(1)};
Input       =   inputdlg( prompt, title, [1 60], defInput);
RunNo       =   str2double(Input(1));

NextRun     = Runs{RunNo};

RunInfo     = [RunNo, numel(Runs), NOccur(RunNo)]; % which run, how many runs in total, which occurrence of run type
    

if  strcmp(prm.exp.RunType , 'test') && RealRun
    subInfo = GetParticipantInfo('age','gender','hand');
end
    

%% setup propixx mode

% eye-link calibaration
if prm.exp.eyelink_live
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram', 0); % 2 for 480, 5 for 1440 Hz, 0 for normal
    Datapixx('RegWrRd');
    
    prm = eyelink_init(prm,id); % initialize the eye tracker parameters
    EyelinkEnterSetup(prm.eyelink);% enter the set up for the eyelink
    Screen('FillRect',prm.monitor.window,127); % make it grey again
end


if prm.exp.is_live % switch to 1440Hz mode
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram', 5); % 2 for 480, 5 for 1440 Hz, 0 for normal
    Datapixx('RegWr');

    Datapixx('SetPropixxLedCurrent', 0, 8);
    Datapixx('SetPropixxLedCurrent', 1, 10);
    Datapixx('SetPropixxLedCurrent', 2, 9);

    Datapixx('RegWrRd');
end

 
% start the eye tracker
if prm.exp.eyelink_live
    Eyelink('StartRecording');
    WaitSecs(0.050);
    Eyelink('Message', 'STARTEXP');
end

%% Run experiment

switch NextRun

    case 'pra'
        Adjustment_withMask_Rift(SubNo, prm);
    case 'test'

end

clearvars; sca; ShowCursor;

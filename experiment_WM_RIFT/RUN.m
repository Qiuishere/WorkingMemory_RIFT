%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% MEG study: probing working memory biases with RIFT                      %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%    Qiu Han                          %
%    Partly based on Songyun Bai's code                                   %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
sca; clear; close all; clc; commandwindow;
addpath(genpath('Functions'));
   
SubNo = 4;

Environment = 1; % 1 - office, 2 - MEG lab. this controls which screen/projector to use

RealRun = 0;  % this controls if use MEG (i.e. send triggers) and store date in Data_Lab

if Environment==2 && RealRun ==1
    useEye = 1; % this controls if use eye-Tracker
else
    useEye = 0;
end


%% define experiment Runs

% experiment information
Runs = {'pra', ...
        'test'};

% enquire about current run
prompt      =   {'Enter next run number (1: practice; 2: test):',...
    'Enter Participant number:',...
    'RealRun?'...
    'Use eye-tracker?'};
title       =   sprintf('Participant %g', SubNo);
defInput    =   {num2str( 1), num2str(SubNo), num2str(RealRun), num2str(useEye)};
Input       =   inputdlg( prompt, title, [1 60], defInput);


%update all entries based on input
RunType     = Runs{str2double(Input{1})};
SubNo       = str2double(Input{2});
RealRun     = str2double(Input{3});
useEye      = str2double(Input{4});


% initiate corresponding parameters
prm = prm_RIFT(SubNo, RunType, RealRun, Environment, useEye);


% gather participant's demographic data
if  strcmp(prm.exp.RunType , 'test') && RealRun
    prm.subInfo = GetParticipantInfo('age','gender','hand');
end
    

%% setup propixx mode

% eye-link calibaration
if useEye
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram', 0); % 2 for 480, 5 for 1440 Hz, 0 for normal
    Datapixx('RegWrRd');
    
    prm = eyelink_init(prm,id); % initialize the eye tracker parameters
    EyelinkEnterSetup(prm.eyelink);% enter the set up for the eyelink
    Screen('FillRect',prm.monitor.window,127); % make it grey again
end


 
% start the eye tracker
if useEye
    Eyelink('StartRecording');
    WaitSecs(0.050);
    Eyelink('Message', 'STARTEXP');
end

%% Run experiment


[prm, T] = Adjustment_withMask_RIFT(prm);

clearvars; sca; ShowCursor;

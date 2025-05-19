Scripts for both sitmuli presentation and MEG data analysis

# Experiment: ./experiment_WM_RIFT
there are a few versions of scripts:

**Adjustment_retroCue_flip**: a classical retro-cue paradigm, behavior only. This paradigm was not used in the final MEG experiment

**Adjustment_withMask_behavior**: Another version of behavior task, without retro-cue

**Adjustment_withMask_RIFT**: the version that was actually used in the final MEG experiment, including interface with MEG and eye tracker. It should be called as a function by RUN.m In this version, a stimulus is drawn at four quadrants of the screen simultaneously, and this process is repeated for three color channels to achieve 1440hz refresh rate

**Adjustment_withMask_RIFT_across_12**: A new way of achieving 1440 hz: 4 quadrants * 3 color channels are drawn all at the same time. This should be the fasted way in terms of code efficiency. But since in reality it doesn't differ from Adjustment_withMask_RIFT much in term of timing, we used the previous one just to ensure the color was presented right.


**UpDownTask***: Behavior scripts for an Up-down task. Was not adopted for MEG experiment


# MEG and eye tracking analysis: ./scripts

## Preprocessing: ./PreProc
in order:

**run_preproc_beforeica**: noise rejection: filtering, motion (manual), muscle noise (manual)

**run_preproc_ica**: simply run ica and save the weights

**run_preproc_afterica**: manually identify noisy ICs and save their lables. This step doesn't save a cleaned version of data after IC rejection. Instead, the saved ICs will be rejected when loading data for further analysis using subfun/load_clean_data.m

## frequency analysis: ./Analysis

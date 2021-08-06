# MonkeysChokeUnderPressure_behavior
Data and analysis scripts for the paper "Monkeys exhibit a paradoxical decrease in performance in high-stakes scenarios", Smoulder, Pavlovsky, Marino et al., 2021. Correspondence on the manuscript should go to Steven Chase (schase@andrew.cmu.edu) or Aaron Batista (apb10@pitt.edu). Questions about data/preprocessing/reproducing figures should go to me (asmoulder3195@gmail.com). Happy to chat about anything regarding the paper or data!

Currently (8/6/21) uploaded are all pre-processed data (https://doi.org/10.6084/m9.figshare.13547435) and analysis + figure generation scripts, run originally in MATLAB r2018a. These are in the "AnalysisAndFigures" directory. 

To reproduce the figures/tables of the manuscript:
- Download the data files from figshare to some folder
- Download the MATLAB files from this repo
- In MATLAB, add the folder with the data files to your path
  - I.e.  addpath(*INSERT DIRECTORY PATH HERE*) 
- Run the scripts for each of the figures/tables!
  
A description of the data files is below. If any of the functions/scripts are unclear or don't work, please let me know!

  
  
Data format:
  
- For speed+accuracy task behavior files (including Rare-Large/Common-Jackpot experiments):
  - Each file is one structure (behavior) that has a length equal to the number of trials
  - The fields are as follows:
    - direction: label for reach target direction. For Monkeys N/F, 1 = 0 degrees from the horizontal, 2 = 180 degrees. For E, 1-8 correspond to 0:45:315 degrees
    - reward: label for reward size. 1 = Small, 2 = Medium, 3 = Large, 4 = Jackpot, 5 = Rare-Large/Common-Jackpot. For corresponding reward values, see Table S1
    - trialStatusLabels: label for trial status. See below for full explanation for each trial's classification. 1 = success, all other values are failures.
    - reactionTimes: the time between go cue and exiting the start target (ms)
    - peakSpeed: max hand speed (m/s)
    - avgReachSpeed: average hand speed from reaction time to end of homing (m/s)
    - centerExitSpeed: hand speed at the time the start target was exited (m/s)
    - t1_trialStartTimes: time when animal's hand entered center target (for all t#_xxxx fields, 0 is when the center target shows up) (ms)
    - t2_targetOnsetTimes: time when the target shows up onscreen. Should be 200ms+/-screen lag from t1 (ms)
    - t3_goCueTimes: time of go cue (ms)
    - t4_reactionTimes: time when 20% of peak speed was achieved (ms)
    - t4b_reactionTimes_20mm: time when hand speed reaches >= 0.02 m/s (ms)
    - t4c_reactionTimes_exit: time when hand exited the center target
    - t5_peakSpeedTimes: time of peak speed (ms)
    - t6_reachEndTimes: time of hand speed going below 10% of peak speed (ms)
    - t6b_reachEndTimes_20mm: time of hand speed going below 0.02 m/s (ms)
    - t7_postReachStateTimes: Either time of reach target entry or time of failure, depending on if the trial ever went into the reach target (ms)
    - t8_trialEndTimes: time of reward or failure
    - kinematics_updated: structure containing kinematics for the trial
      - position: ntime x 2 matrix with XY position of cursor on screen (mm)
      - velocity: ntime x 2 matrix with cursor velocity (mm/ms)
      - acceleration: ntime x 2 matrix with cursor acceleration (mm/ms^2)
      - time: ntime x 1 vector of time values (ms)
      - rotatedPosition: ntime x 2 matrix of cursor position that has been rotated such that the reach target location is at [X Y] = [85 0]mm (mm)
      - distanceFromCenter: ntime x 1 vector of cursor distance from the center of the center target (mm)
      - distanceFromEndTarget: ntime x 1 vector of cursor distance from the center of the end target (mm)
      - speed: ntime x 1 vector of cursor speeds
      - ballisticPredEndpoint: 1 x 2 XY position of the ballistic endpoint prediction (see text for how this was calculated) (mm)
      - rotatedBallisticPredEndpoint: 1 x 2 ballistic endpoint prediction, rotated such that the end target location for the trial is at [X Y] = [85 0]mm (mm)
    - stateTrans: ignore; was used in preprocessing to get out event timepoints
    - distInTarget: how much distance was covered by the hand while in the reach target (mm)
    - timeInTarget: how long the hand was in the reach target (ms)
    - trialIndex: trial number for a given day
    - day: session number
    - reach_timing: structure with reach timing information:
      - reactionTime: time hand exited center target (same as reactionTimes in the full structure) (ms)
      - earlyTime: time from exiting the center target to covering 1/3 of the distance to the reach target (ms)
      - midreachTime: time between 1/3 to 2/3 of distance to the target (ms)
      - homingTime: time from 2/3 of distance to the target to within 1mm of entry (ms)
  
- For precision task ("tubes") experiment data:
    - Each file is one structure (behavior) that has a length equal to the number of trials
    - Fields follow the same conventions as those in the speed+accuracy task
  
- For choice task experiments:
    - Most comparisons are for one reward versus another in a nreward x nrewards matrix. 1-4 = SMLJ, 5 = CJ/RL
    - Each file has one structure, choiceBehavior, that contains the following field:
    - comboFracCorrect: nreward x nrewards upper-triangular matrix with the fraction of the time the correct (higher) reward was selected for the given comparison
      - the indices of these square matrices indicate the comparison. For instance, position (1,3) compares Small versus Large
    - comboN: nreward x nrewards upper-triangular matrix with the number of times each combination of reward conditions was presented
    - comboNCorrect: nrewards x nrewards upper-triangular matrix with the number of times the correct (higher) reward was selected
    - decision0States: ignore; from preprocessing of data to identify if one target vs. other was selected
    - decision1States: ignore; from preprocessing of data to identify if one target vs. other was selected
    - decisionChoices: ntrials x 1 vector of labels for which decision was made (0 or 1)
    - dirAngles: angles for target presentation corresponding to the directionOptions
    - directionOptions: ntrials x 2 matrix of labels with target locations for choice 0 (first col) or 1 (2nd col)
    - directions: unique direction option values
    - rewardOptions: ntrials x 2 matrix of reward labels for choice 0 (first col) or 1 (2nd col)
    - sameRewDirCounts: ndirections x 1 vector showing how many times a given direction was selected when the same reward was presented at both targets
      - this is useful for evaluating if the animals have bias towards one target or another. Monkey N is a fantastic example of this...
    - sameRewDirFracs: ndirections x 1 vector showing fraction of times a given direction was selected when the same reward was presented at both targets
    - wrongDecSelectedDirs: ndirections x 1 vector showing the total number of times a wrong decision was made towards the given target
  
  Other:
    - SAVESTATE_forBREPExample2.mat has preprocessed information from a single trial that's used for the schematics in Fig. 3A
  
  
Trial Status Labels are as follows (% of all trials for all subjects in normal speed+accuracy task):
  - 1 = success (69.85%)
  - -11 = false start, hand began reach towards target before go cue (1.85%)
  - -12 = delay drift, cursor drifted slowly out of the center target before the go cue. Often (but not always) "cheating" towards the target (4.25%)
  - -22 = overshoot (1), trials where the cursor goes past the center of (but never enters) the reach target (0.97%)
  - -23 = undershoot (1), hand speed dropped below 10% of peak speed outside of the center target, ran out of time before corrective movement achieved target (2.77%)
  - -24 = undershoot (2), based on time peak speed, subject was too slow to have any chance of making it to the target in time (4.33%)
  - -25 = undershoot (3), based on time of peak speed, subject could have made it to target with his fastest-observed homing time, but did not make it in time (4.09%)
  - -31 = overshoot (2), trials where the cursor "scuffs" the target, grazing the inside for a short period of time (1.21%)
  - -32 = overshoot (3), cursor "blows through" the target, coming out the far side at >10% of peak speed (5.00%)
  - -34 = target hold drift (1), cursor lands in target but drifts out before reward (2.77%)
  - -35 = target hold drift (2), cursor lands in target near edge and "jitters" out before reward (0.59%)
  
The following are trial status labels that are not analyzed in this manuscript and are removed before analyses
  - 0/-111/-14 = error, trial that did not occur due to system glitch/lag or had cursor drops; should be removed before analysis (1.00%)
  - -13/-20 = quitout, Large hand motion away from the reach target, indicating quitting (0.51%)
  - -21 = no attempt, cursor basically never leaves center target (reaction time > 500ms) (0.65%)
  - -33 = early return, cursor begins returning towards center (speed > 0.1 m/s) and exits reach target before reward (0.05%)

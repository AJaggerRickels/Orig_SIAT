README for the Implicit Association Task

*****************************************************************SUMMARY**************************************************************************


*****************************************************************TRIAL STRUCTURE******************************************************************

      [         instructions         ] ->
      ^                              ^
INSTRUCT_ONSET                   TASK_ONSET


        ~3s                   2s window of collecting response  
[fixation shown] -> [target word] + [left word shown} + [right word shown] -> 
^                   ^              
TRIAL_ONSET        TARGET_ONSET      RESPONSE


            10s
      [block break] -> 
      ^     
BREAK 

*****************************************************************INPUT DETAILS********************************************************************

EACH LINE CODES: one trial
COLUMN 1: Trial Type (block_type)
COLUMN 2: Target Word
COLUMN 3: Duration 
COLUMN 4: Left Word
COLUMN 5: Right Word
COLUMN 6: Correct Category
TRIAL ORDER IS: randomized within block

RANDOMIZATION SCHEME:

Trial_Type randomized order.
Fixation durations are chosen randomply from a list of 14 durations. Average 3 seconds

*****************************************************************OUTPUT DETAILS*******************************************************************

INSTRUCT_ONSET (1)
response_time: not used
response: not used
result: not used

TASK_ONSET (2)
response_time: time between INSTRUCT_ONSET and TASK_ONSET
response: not used
result: not used

TRIAL_ONSET (3)
response_time: not used
response: not used
result: not used

FIXATION_ONSET (4)
response_time: not used
response: not used
result: duration of the cue

TARGET_ONSET (5)
response_time: not used
response: not used
result: the target word


RESPONSE (7)
response_time: time since TARGET_ONSET
response: -1 for left, 1 for right
result: 0 for incorrect, 1 for correct


BREAK (9)
response_time: not used
response: not used
result: not used

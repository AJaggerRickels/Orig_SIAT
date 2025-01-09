
import StimToolLib, os, random, operator
from psychopy import visual, core, event, data, gui, sound

# Suicide Implicit Association Task (IAT)

class GlobalVars:
    #This class will contain all module specific global variables
    #This way, all of the variables that need to be accessed in several functions don't need to be passed in as parameters
    #It also avoids needing to declare every global variable global at the beginning of every function that wants to use it
    def __init__(self):
        self.win = None #the window where everything is drawn
        self.clock = None #global clock used for timing
        self.x = None #X fixation stimulus
        self.output = None #The output file
        self.msg = None
        self.ideal_trial_start = None #ideal time the current trial started
        self.this_trial_output = '' #will contain the text output to print for the current trial
        self.trial = None #trial number
        self.trial_type = None #current trial type
        self.offset = 0.008 #8ms offset--request window flip this soon before it needs to happen->should get precise timing this way
        self.break_instructions = ['''You may now take a short break.''']
        self.line_location_range = 24
        #self.line_location_range = 0.7 #amount the lines (vertical and horizontal) can move left/right and up/down


event_types = {
    'INSTRUCT_ONSET':1,
    'TASK_ONSET':2,
    'TRIAL_ONSET':3,
    'FIXATION_ONSET':4, 
    'TARGET_ONSET':5,
    'RESPONSE':6,
    'BREAK' : 7,
    'TASK_END':StimToolLib.TASK_END
    }


def show_fixation(trial_start, duration):
    """
    Show Fixation 
    """

    StimToolLib.mark_event(g.output, g.trial, g.trial_type, event_types['FIXATION_ONSET'], trial_start, 'NA', 'NA', duration)#, g.session_params['signal_parallel'], g.session_params['parallel_port_address'])
    g.fixation.draw()
    g.win.flip()
    StimToolLib.just_wait(g.clock, trial_start + duration - g.offset) # wait



def addWordLine(word):
    """
    Add new line to the words that includ an 'or'
    return Word with new Line
    """
    newword = ''
    words = word.split(' or ')
    
    # Check if there's an or in the word
    if len(words) < 2:
        newword = word
    else:
        newword = words[0] + '\n' + 'or\n' + words[1]
    return newword


def do_one_trial(trial):
    """
    Function for displaying 1 trial
    """ 
    
    g.trial_type = trial[0]
    target_word = trial[1]
    duration = int(trial[2])
    left_word = addWordLine(trial[3])
    right_word = addWordLine(trial[4])
    correct_word = trial[5]

    g.win.flip()
    trial_start = g.clock.getTime()
    #mark trial start
    StimToolLib.mark_event(g.output, g.trial, g.trial_type, event_types['TRIAL_ONSET'], trial_start, 'NA', 'NA', 'NA')#, g.session_params['signal_parallel'], g.session_params['parallel_port_address'])

    show_fixation(trial_start, duration)

    
    left_stim = visual.TextStim(g.win, text = left_word, units = 'norm', pos = (-0.55,0.4), height = .2)
    right_stim = visual.TextStim(g.win, text = right_word, units = 'norm', pos = (0.55, 0.4), height = .2)

    left_stim.draw()
    right_stim.draw()

    target_stim = visual.TextStim(g.win, text = target_word, units = 'norm', pos = (0,-0.3), height = .2)
    target_stim.draw()

    g.win.flip()
    target_start = g.clock.getTime()
    
    StimToolLib.mark_event(g.output, g.trial, g.trial_type, event_types['TARGET_ONSET'], g.clock.getTime() , '', '', target_word)#, g.session_params['signal_parallel'], g.session_params['parallel_port_address']) 
    
    event.clearEvents() # Clear Events in the event buffer
    notResponsed = True
    while g.clock.getTime() < target_start + 2: #show target for 2.5s
        StimToolLib.check_for_esc()
        k = event.getKeys([g.session_params[g.run_params['select_1']], g.session_params[g.run_params['select_2']]]) 
        if k and notResponsed:
            time_of_response = g.clock.getTime()
            if k[0] == g.session_params[g.run_params['select_1']]:
                resp = -1

                # Show feedback. change color of the text
                left_stim.setColor('yellow')
                left_stim.draw()

                right_stim.setColor('white')
                right_stim.draw()   

                target_stim.draw()
            else:
                resp = 1
                # Show feedback
                
                left_stim.setColor('white')
                left_stim.draw()

                right_stim.setColor('yellow')
                right_stim.draw()

                target_stim.draw()

            g.win.flip()
            if correct_word == '1': #Left is correct category
                if k[0] == g.session_params[g.run_params['select_1']]:
                    correct = 1
                else:
                    correct = 0
            else: #Right is correct category
                if k[0] == g.session_params['down']:
                    correct = 1
                else:
                    correct = 0
            notResponsed = False # Set Flag to only let them respond once
            StimToolLib.mark_event(g.output, g.trial, g.trial_type, event_types['RESPONSE'], time_of_response, time_of_response - target_start, resp, correct)#, g.session_params['signal_parallel'], g.session_params['parallel_port_address'])
            
    
    g.trial = g.trial + 1 # increment trial number 

def run(session_params, run_params):
    global g
    g = GlobalVars()
    g.session_params = session_params
    g.run_params = StimToolLib.get_var_dict_from_file(os.path.dirname(__file__) + '/IAT.Default.params', {})
    g.run_params.update(run_params)
    
    print os.path.exists(os.path.dirname(__file__) + '/IAT.Default.params')
    try:
        run_try()
        g.status = 0
    except StimToolLib.QuitException as q:
        g.status = -1
    StimToolLib.task_end(g)
    return g.status
        
def sort_and_shuffle_blocks(trial_types, targets, durations, left_words, right_words, correct_words):
    # Shuffle the images, left word and right word and correct responses
    # per Block
    g.block_1 = []
    g.block_2 = []
    g.block_3 = []
    g.block_4 = []
    for ttype,target,duration,left_w,right_w,correct_w in zip(trial_types, targets, durations, left_words, right_words, correct_words):
        trial = [ttype, target, duration, left_w, right_w, correct_w]
        if '1' in ttype:
            g.block_1.append(trial)
        if '2' in ttype:
            g.block_2.append(trial)
        if '3' in ttype:
            g.block_3.append(trial)
        if '4' in ttype:
            g.block_4.append(trial)
        
    #shuffle all images lists
    random.shuffle(g.block_1)
    random.shuffle(g.block_2)
    random.shuffle(g.block_3)
    random.shuffle(g.block_4)
    
def showBreakBlock():
    g.win.flip() # clear window
    now  = g.clock.getTime()

    # Show the Text
    text_stim = visual.TextStim(g.win, text = 'This is a new block.\nCategories will now change', units = 'norm', pos = (0,0), height = .1)
    text_stim.draw()

    g.win.flip() # Flip screen

    StimToolLib.just_wait(g.clock, now + 5 - g.offset) # wait 5 seconds
    StimToolLib.mark_event(g.output, 'NA', 'NA', event_types['BREAK'], now, 'NA', 'NA', 'NA')#, g.session_params['signal_parallel'], g.session_params['parallel_port_address'])

    # Show Fixation
    g.fixation.draw()
    g.win.flip()
    now = g.clock.getTime()
    StimToolLib.just_wait(g.clock, now + 5 - g.offset) # wait 5 seconds
    StimToolLib.mark_event(g.output, 'NA', 'NA', event_types['FIXATION_ONSET'], now, 'NA', 'NA', 'NA')#, g.session_params['signal_parallel'], g.session_params['parallel_port_address'])



def run_try():  
#def run_try(SID, raID, scan, resk, run_num='Practice'):
    schedules = [f for f in os.listdir(os.path.dirname(__file__)) if f.endswith('.schedule')]
    if not g.session_params['auto_advance']:
        myDlg = gui.Dlg(title="IAT")
        myDlg.addField('Run Number', choices=schedules, initial=str(g.run_params['run']))
        myDlg.show()  # show dialog and wait for OK or Cancel
        if myDlg.OK:  # then the user pressed OK
            thisInfo = myDlg.data
        else:
            print 'QUIT!'
            return -1#the user hit cancel so exit 
        g.run_params['run'] = thisInfo[0]
    StimToolLib.general_setup(g)

    g.fixation = visual.TextStim(g.win, text = '+', units = 'norm', pos = (0,0), height = .3) # Save the fixation stmiulus
    g.select = visual.TextStim(g.win, text = '----', units = 'norm', pos = (0,0), height = 0.2) # Save select symbol


    schedule_file = os.path.join(os.path.dirname(__file__), g.run_params['run'])
    #param_file = os.path.join(os.path.dirname(__file__),'T1000_DP_Schedule' + str(g.run_params['run']) + '.csv')

    # No images used for schedul so just read schedul as usual

    schedule = open(schedule_file, 'r')
    trial_types = []
    targets = []
    durations = []
    left_words = []
    right_words = []
    correct_words = []

    print schedule
    for idx,line in enumerate(schedule):
    
        if idx == 0:
            continue
        row = line.replace('\n','').split(',')
        trial_types.append(row[0])
        targets.append(row[1])
        durations.append(row[2])
        left_words.append(row[3].replace('"',''))
        right_words.append(row[4].replace('"',''))
        correct_words.append(row[5])

    
    # Shuffle by blocks
    sort_and_shuffle_blocks(trial_types, targets, durations, left_words, right_words, correct_words)


    #load instruction slides
    #slides = []
    #for i in range(1,6):
    #    slides.append(visual.ImageStim(g.win, image=os.path.join(os.path.dirname(__file__),  'media/Instructions/slide' + str(i) + '.bmp'), pos=[0,0], units='pix'))
    
    
    
    #g.x = visual.TextStim(g.win, text="+", units='pix', height=25, color=[-1,-1,-1], pos=[0,0], bold=True)
    g.clock = core.Clock()
    start_time = data.getDateStr()
    param_file = g.run_params['run'][0:-9] + '.params' #every .schedule file can (probably should) have a .params file associated with it to specify running parameters (including part of the output filename)
    StimToolLib.get_var_dict_from_file(os.path.join(os.path.dirname(__file__), param_file), g.run_params)
    g.prefix = StimToolLib.generate_prefix(g)
    fileName = os.path.join(g.prefix + '.csv')
    #g.prefix = 'DP-' + g.session_params['SID'] + '-Admin_' + g.session_params['raID'] + '-run_' + str(g.run_params['run']) + '-' + start_time
    #fileName = os.path.join(os.path.dirname(__file__), 'data/' + g.prefix +  '.csv')
    g.output = open(fileName, 'w')
    
    sorted_events = sorted(event_types.iteritems(), key=operator.itemgetter(1))
    g.output.write('Administrator:,' + g.session_params['admin_id'] + ',Original File Name:,' + fileName + ',Time:,' + start_time + ',Parameter File:,' +  schedule_file + ',Event Codes:,' + str(sorted_events) + ',Trial Types are coded as follows: 8 bits representing [valence neut/neg/pos] [target_orientation H/V] [target_side left/right] [duration .5/1] [valenced_image left/right] [cue_orientation H/V] [cue_side left/right]\n')
    g.output.write('trial_number,trial_type,event_code,absolute_time,response_time,response,result\n')
    StimToolLib.task_start(StimToolLib.IAT_CODE, g)
    instruct_start_time = g.clock.getTime()
    StimToolLib.mark_event(g.output, 'NA', 'NA', event_types['INSTRUCT_ONSET'], instruct_start_time, 'NA', 'NA', 'NA')#, g.session_params['signal_parallel'], g.session_params['parallel_port_address'])
    #StimToolLib.show_title(g.win, g.title)
    #g.output.write('Trial_Type:-1:negative;0:neutral;1:positive,Image,ITI_Onset,ITI_startle,Stimulus_onset,Stimulus_startle,Valence_Rating,Valence_rating_time,Arousal_rating,Arousal_rating_time\n')
    
    # Run 
    #StimToolLib.run_instructions(os.path.join(os.path.dirname(__file__), g.run_params['instruction_schedule']), g)
    StimToolLib.run_instructions_keyselect(os.path.join(os.path.dirname(__file__), g.run_params['instruction_schedule']), g)

    g.trial = 0
    if g.session_params['scan']:
        StimToolLib.wait_scan_start(g.win)
        g.clock.reset()
        start_time = g.clock.getTime()
    else:
        StimToolLib.wait_scan_start(g.win)
        g.clock.reset()
        start_time = g.clock.getTime()
    g.win.flip()
    instruct_end_time = g.clock.getTime()
    StimToolLib.mark_event(g.output, 'NA', 'NA', event_types['TASK_ONSET'], instruct_end_time, instruct_end_time - instruct_start_time, 'NA', 'NA')#, g.session_params['signal_parallel'], g.session_params['parallel_port_address'])
    
    
    g.trial = 0
    ## Run Trials
    # Block 1
    for trial in g.block_1:
        do_one_trial(trial)
    
    showBreakBlock()

    for trial in g.block_2:
        do_one_trial(trial)
    
    showBreakBlock()

    for trial in g.block_3:
        do_one_trial(trial)
    
    showBreakBlock()

    for trial in g.block_4:
        do_one_trial(trial)
    
    
  
 



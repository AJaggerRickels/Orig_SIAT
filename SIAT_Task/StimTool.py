from psychopy import prefs
prefs.general['audioLib'] = [u'pyo', u'pygame']
#prefs.general['audioLib'] = [u'pygame']
#prefs.general[u'audioDriver'] = [u'ASIO4ALL', u'ASIO', u'Audigy']
from pygame import pypm #need to load this before loading psychopy.sound when using pyo 
from psychopy import core, visual, gui, data, misc, event, logging #, microphone VMeter sound
import time, numpy, random, datetime, sys, os, copy, pyo, csv, ast, threading, pygame.pypm, CloseAudio
import StimToolLib
import IAT.IAT


mod_mapping = {
    'Implicit Association Task' : IAT.IAT
    }

def run_task_until_success(task, session_params):
    
    task_and_args = task.split(':') #task names (followed by args) should end with a :
    task = task_and_args[0]
    #if task == 'Stop Camera': #turn off the camera early--typically just before the end so it no longer ties up the output log file
    #    CD.stop_capture()
    #    core.wait(0.2)
    #    return False
    if len(task_and_args) > 1:
        run_params = StimToolLib.convert_run_args_to_dict(task_and_args[1])
    else:
        run_params = {}
    success = False
    just_escaped = False
    while True: #repeat the task until the run is completed successfully
        if session_params['auto_advance'] and not just_escaped: #auto advance to next task--don't ask which task to go to
            status = mod_mapping[task].run(session_params, run_params)
            if status == 0: 
                return False#task completed successfully, will continue to the next one
            else:
                just_escaped = True
            #maybe put a mechanism here to convert from auto_advance to not--so that the experimenter can switch order mid session
        else: #don't auto advance--open a dialogue and ask for the next task to run
            myDlg = gui.Dlg(title="StimTool: File List Mode, reading task order from the specified file")
            myDlg.addField('Program', choices=modules, initial=task)
            myDlg.show()  # show dialog and wait for OK or Cancel
            if myDlg.OK:  # then the user pressed OK
                thisInfo = myDlg.data
                if thisInfo[0] == task: #if the user picked the defualt next task, continue reading from the file    
                    switched = False
                elif thisInfo[0] == 'SKIP': #user chose to skip this task--go on to the next one but don't change to free mode
                    return False
                elif thisInfo[0] == 'SKIP TO': #user picked skip to--return 2 to signal to prompt the user with a list to skip to
                    return 2
                else:
                    switched = True #otherwise set it to free mode
                    session_params['auto_advance'] = False
                if thisInfo[0] == task: #the user selected the default: load any parameters that might be in the file
                    these_params = run_params
                else:
                    these_params = {}
                status = mod_mapping[thisInfo[0]].run(session_params, these_params)
                
                if status != -1: #-1 is returned when a task fails (e.g. user hits escape to quit)
                    return switched
            else:
                print 'QUITTING!'
                return True #say switched task->will start free mode, one more escape to quit
        

def reset_flags():
    flag_file = open('run_flags.txt','w')
    flags={'AAC_reset_total': True}
    for f in flags:
        flag_file.write(f+' '+str(flags[f])+'\n')
    flag_file.close()

if __name__ == '__main__':
    #CD.start_camera()
    #StimToolLib.open_and_close_vmeter() #this hack seems to fix a problem with the vMeter not responding the first time it's used after logging in...
    reset_flags()

    modules = mod_mapping.keys()
    modules.append('SKIP')
    modules.append('SKIP TO')
    modules.sort() #sorted alphabetically so tasks are easier to find
    
    #Enter ID and select tasks dialog
    myDlg = gui.Dlg(title="StimTool")
    myDlg.addField('Subject ID:', StimToolLib.get_var_from_file('Default.params', 'last_subject')) #subject ID
    myDlg.addField('Administrator ID:', StimToolLib.get_var_from_file('Default.params', 'last_admin')) #name of the person administering the session
    task_lists = [f for f in os.listdir('.') if f.endswith('.TL')] #get a list of all ".TL" files (tasklists), which have a list of all tasks to be run in a session
    task_lists.insert(0, 'free')
    myDlg.addField('Experiment Order:', choices=task_lists) #text file with a list of which tasks to run.  Some tasks may be repeated.
    myDlg.show()  # show dialog and wait for OK or Cancel
    if myDlg.OK:  # then the user pressed OK
        thisInfo = myDlg.data
    else:
        print 'QUIT!'
        core.quit()#the user hit cancel so exit 
    sid = thisInfo[0].replace('\'', '.')
    raID = thisInfo[1].replace('\'', '.') #prevent problems with storing/reading strings that have single quotes in them
    try:
        StimToolLib.write_var_to_file('Default.params', 'last_subject', sid)
        StimToolLib.write_var_to_file('Default.params', 'last_admin', raID)
    except:
        print "Could not store last subject and administrator--most likely because the user doesn't have write permission to the StimTool directory."
    session_params = StimToolLib.get_var_dict_from_file('Default.params', {'SID':sid, 'raID':raID}) #initialize default session parameters--SID and raID always come from the dialogue box
    #session_params['vMeter'] = StimToolLib.try_to_open_vMeter()
    session_params['admin_id'] = raID
    
    
    if thisInfo[2] == 'free':
        free = True
    else:
        #read the input file to get the experiment order
        #this file should contain one line /per task, and each line must match an entry in the mod_mapping dictionary    
        with open(thisInfo[task_list_idx]) as f:
            order = f.read().splitlines()
            session_params['task_list']=thisInfo[task_list_idx]
        free = False
        param_file = thisInfo[task_list_idx][0:-3] + '.params' #every .TL file can have a .params file associated with it to specify running parameters
        session_params = StimToolLib.get_var_dict_from_file(param_file, session_params) #overwrite any parameters defined in the session specific file (e.g. BehavioralSession.params)
    #redirect stdout and stderr to a file
    if session_params['redirect_output']:
        StimToolLib.redirect_output(session_params)
    #if session_params['signal_parallel']: #if using the parallel port, make sure it's working
    #    StimToolLib.verify_parallel(session_params)
    
    if session_params['scan']:
        StimToolLib.get_exam_number(session_params) #if it's a scanning session, save the exam number (e.g. S2352)

#    if session_params['record_video']: #start the camera
#        #StimToolLib.initialize_mic()
#        if not CD.start_and_test_camera():
#            StimToolLib.error_popup('ERROR: CAMERA IS NOT WORKING PROPERLY.  UNPLUG THE CAMERA FROM THE USB PORT AND PLUG IT BACK IN, THEN RESTART')

    
    
    
    if not free:
        idx = 0
        while idx < len(order):
            switched_task = run_task_until_success(order[idx], session_params)
            if switched_task == 1: #user switched tasks--go to free mode
                free = True
                break
            elif switched_task == 2: #user picked SKIP TO: so prompt to choose which task to skip to
                myDlg = gui.Dlg(title="StimTool: What task would you like to skip to?")
                myDlg.addField('Program', choices=order[idx:])
                myDlg.show()  # show dialog and wait for OK or Cancel
                if myDlg.OK:  # then the user pressed OK
                    thisInfo = myDlg.data
                    idx = order.index(thisInfo[0]) #this won't work if there are two runs of the same task with the same run---but that shouldn't happen anyway
                else:
                    print 'QUITTING!'
                    break
            else: #advance to the next one
                idx = idx + 1
                    #choose the task/run to skip to here
        #for task in order:
        #    switched_task = run_task_until_success(task, session_params)
        #    if switched_task:
        #        free = True
        #        break
    if free:
        scan = True 
        while True:
            myDlg = gui.Dlg(title="StimTool: Free Mode, hit cancel or escape when done")
            myDlg.addField('Program', choices=modules)
            myDlg.show()  # show dialog and wait for OK or Cancel
            if myDlg.OK:  # then the user pressed OK
                thisInfo = myDlg.data
                status = mod_mapping[thisInfo[0]].run(session_params, {})
            else:
                print 'QUIT!'
                break #core.quit()
    #CD.stop_capture()
    #CloseAudio.run() #seems to close the audio server so the user doesn't have to hit the stop button...?
    #pygame.pypm.Terminate()

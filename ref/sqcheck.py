# filename: sqcheck.py

# to run:
#	python sqcheck.py [user] [sleepInterval_in_seconds]

# Author:
#	Daniel W. Barnette
# 	SNL org 1422
#	dwbarne@sandia.gov  or  dwbarne@gmail.com

# Purpose:
#  checks on current machine's squeue (SLURM command) status for any owner's 
#    jobs, current user if unspecified on command line
# 
#  keeps checking queues until jobs are completed or until terminated by user
#   with CTRL-C
#
#  if all jobs complete, then this script prints out a notice and ends
#
#  default sleep interval between queue checks: 5 seconds
#   -- user-specified sleep interval can be specified on command line

import os, time, sys, platform, subprocess

DEBUG=0

def timeStamp():
# time start date
    timeStartDate = (2012,1,1,0,0,0,6,1,0) # jan 1, 2012
    timeStartDateSinceEpoch_Seconds = time.mktime(timeStartDate)
    timeNow_Seconds = time.time()
    timeNowSinceStartDate_Days = \
        int((timeNow_Seconds - timeStartDateSinceEpoch_Seconds)/3600./24. + 1)
# set date and time, but update before sending to table
    date = time.ctime(time.time())[4:10] + ' ' + time.ctime(time.time())[20:24]
    month = time.ctime(time.time())[4:7]
    dayofweek = time.ctime(time.time())[0:3]
    dayofmonth = time.ctime(time.time())[8:10]
# replace blanks with zeros for printing
    dayofmonth = dayofmonth.replace(' ','0')
    year = time.ctime(time.time())[20:24]
    timeToday = time.ctime(time.time())[11:19]      
# define month dictionary
    dictMonth = {
        'Jan':'01', 'Feb':'02', 'Mar':'03', 'Apr':'04', 'May':'05', 'Jun':'06', 
        'Jul':'07', 'Aug':'08', 'Sep':'09', 'Oct':'10', 'Nov':'11', 'Dec':'12'
        }
# NOTE: for sorting purposes, use date format of yyyy/mm/dd when storing to database
    dateTimeToday = year + '/' + dictMonth[month] + '/' + dayofmonth  + ', ' + timeToday
    
    return dateTimeToday

    
# ==== MAIN ====

print('\n           ***** WELCOME TO SQCHECK.PY *****')

numArgs = len(sys.argv)

hostName = platform.uname()[1]

stringUsage = (
    '\nUSAGE\n' +
    '      $python sqcheck.py [user] [sleepInterval]\n' +
    '   where \'user\' and \'sleepInterval\' are optional.\n\n' +
    '     user: sqcheck will check the status of jobs belonging to \'user\';\n' +
    '             if not specified, will default to username of person\n' +
    '             running sqcheck.\n\n' +
    '     sleepInterval: specifies the number of seconds between checking a job;\n' +
    '             if sqcheck does not find \'user\' jobs in 3 iterations, sqcheck will\n' +
    '             print that no jobs were found and auto exit.\n' +
    '             If not specified, will be, this value is set to\n' +
    '               DEFAULT VALUE: 3 seconds\n' 
    )
if numArgs > 3:
    print(stringUsage)
    sys.exit()
if len(sys.argv) == 2 and sys.argv[1] == 'usage':
    print(stringUsage)
    sys.exit()
        
if DEBUG:
# for 1 command line argument (program name always counts as first argument)
    if numArgs == 2:
        try:
            print('\ntype(int(float(sys.argv[1]))) = %s' % type(int(float(sys.argv[1]))))
        except:
            print('\ntype(sys.argv[1]) = %s' % type(sys.argv[1]))

# for 2 command line arguments
    if numArgs == 3:
        try:
            print('\ntype(int(float(sys.argv[1]))) = %s' % type(int(float(sys.argv[1]))))
        except:
            print('\ntype(sys.argv[1]) = %s' % type(sys.argv[1]))

        try:
            print('\ntype(int(float(sys.argv[2]))) = %s' % type(int(float(sys.argv[2]))))
        except:
            print('\ntype(sys.argv[2]) = %s' % type(sys.argv[2]))
            
# CHECKS
# check if 2 ints or 2 strings have been entered; if so, print 'usage' and quit 
if numArgs == 3:
    try:
        if type(int(float(sys.argv[1]))) == int and type(int(float(sys.argv[2]))) == int:
            print(stringUsage)
            sys.exit()
    except:
        try:
            if type(sys.argv[1]) == str and type(sys.argv[2]) == str:
                print(stringUsage)
                sys.exit()
        except:
            pass
        
# get username whether on Linux or Windows
try:
    userNameTemp = os.environ['USER']
except:
    try:
        userNameTemp = os.environ['USERNAME']
    except:
        userNameTemp = 'UNK'
        stringNoUserName = (
            'No username can be found!\n\n' +
            ' ... exiting\n\n'
            )
        print('\n' + stringNoUserName)
        sys.exit()

if numArgs == 1:
# no values on command line other than program name; set to default values
    userName = userNameTemp
    sleepInterval = 5 
    
elif numArgs == 2:
# either 'user' or 'sleepInteval' has been specified -- figure out which one, 
# ... then use default value for the other
    try:
        value = int(float(sys.argv[1]))
        sleepInterval = value # value is an integer
        userName = userNameTemp
    except:
        userName = sys.argv[1]  # sys.argv[2] is an integer
        sleepInterval = 3
        
else:
    try:
        value = int(float(sys.argv[1]))
        sleepInterval = value
        userName = sys.argv[2]
    except:
        try:
            sleepInterval = int(float(sys.argv[2]))
            userName = sys.argv[1]
        except:
            stringSysArgError = (
                '\nOne of the input command line values must be an integer (sleepInterval),\n'  +
                'and one must be a string (user)\n\n' +
                'For options, type\n\n' +
                '    python sqcheck.py usage\n\n' 
                )
            print(stringSysArgError)
            sys.exit()

stringHeader = (
	'\nRunning sqcheck script\n' +
    '    hostname: %s\n' +
    '    user: %s\n' +
    '    sleep interval: %s sec.'
    ) % (hostName,userName,str(sleepInterval))
print(stringHeader)

timeStart = timeStamp()
print('\nsqcheck date/time start: %s\n' % timeStart)
    
jobRunning = 1
timerCount = 1
totalCount = 0
timerCount_Quit = 3
variablesPerLine = 8
linesEstimatedStartFinish = []

# clear the window
subprocess.call('clear') 

while jobRunning:
    totalCount += 1
# main command for checking queue
# ... open file filled with output from the squeue command
    qC = os.popen('squeue -u %s' % userName)
    
# ... read file as one string
    queueCheck = qC.read()

# split file into separate lines that we can check
    jobList = queueCheck.split('\n') 
    
# remove any blanks at the end of the list
    while not jobList[-1]:
        jobList.pop()
    
    timeCurrent = time.ctime(time.time())[11:19]      
    print('\nqueueCheck #%s on %s at %s' % (totalCount,hostName,timeCurrent))
    print('approx update interval: %s secs\n' % sleepInterval)
    if len(jobList) == 1:
        print(' ..... No jobs are in queue (try %s of %s tries)' % 
          (timerCount,timerCount_Quit))
        if timerCount >= timerCount_Quit :
            print(
                '\n\nNo jobs found in queue on %s in last %s tries.\n\n' +
                'Script sqcheck.py is exiting' 
                ) % (hostName,timerCount_Quit)
            print('\n\nsqcheck date/time start: %s' % timeStart)
            timeFinish = timeStamp()
            print('\nsqcheck date/time finish: %s\n' % timeFinish)
            sys.exit()
        timerCount += 1

    else:

        for line in enumerate(jobList):
            if line[0] == 0:
                print('   %s' % line[1])
            elif line[1] <> '':
                print('%s. %s' % (line[0],line[1]))
# now print status for all jobs
        print('\n-------------------- JOB STATUS ----------------------')
        lineCount = 0
        for line in jobList:
            if line <> '' and lineCount > 0:
                jobID = line.split()[0]
                lESF = os.popen('showstart %s' % jobID)
                linesESF = lESF.read()
                linesEstimatedStartFinish = linesESF.split('\n')
# following if statement is necessary if jobID is found in list but has not yet started 
                if len(linesEstimatedStartFinish) > 1:
                    print('\n' + str(lineCount) + ', ' + linesEstimatedStartFinish[0])
                    print('     ' + linesEstimatedStartFinish[2])
                    print('     ' + linesEstimatedStartFinish[3])
                else:
                    print(str(lineCount) + '. ' + '... job ' + jobID + ' is starting')
            lineCount += 1 

        timerCount = 1
        
# wait before updating job status again
    if timerCount <= timerCount_Quit:
        time.sleep(sleepInterval)

# clear screen and go back up to 'while jobRunning' again 
        subprocess.call('clear') 

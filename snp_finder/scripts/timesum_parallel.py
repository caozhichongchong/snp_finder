import os,glob
import os
allinput_scripts = glob.glob('/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/parallel/*.sh')

def time_sec(realtime):
    realtime = realtime.split(':')
    timeinsec = 0
    timeinsec += float(realtime[-1])
    timeinsec += float(realtime[-2])*60
    if len(realtime) > 2:
        timeinsec += float(realtime[-2]) * 360
    return timeinsec

for input_script in allinput_scripts:
    alltime = []
    errfile = input_script + '.err'
    i = 40
    for lines in open(errfile, 'r'):
        if 'Elapsed' in lines:
            realtime = lines.split('):')[1].split('\n')[0].replace(' ', '')
            alltime.append(
                '%s\t%.2f\n' % (i,time_sec(realtime)))
            i -= 1
    f1 = open('%s.timesum.txt'%(input_script),'w')
    f1.write('thread\ttime\n')
    f1.write(''.join(alltime))
    f1.close()

import glob
import os
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/SNP_model'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/SNP_model_noindel'

def time_sec_old(realtime):
    timeinsec = 0
    if 'h' in realtime:
        timeinsec += float(realtime.split('h')[0])*360
        realtime = realtime.split('h')[1]
    if 'm' in realtime:
        timeinsec += float(realtime.split('m')[0])*60
        realtime = realtime.split('m')[1]
    if 's' in realtime:
        timeinsec += float(realtime.split('s')[0])
    return timeinsec

def time_sec(realtime):
    realtime = realtime.split(':')
    timeinsec = 0
    timeinsec += float(realtime[-1])
    timeinsec += float(realtime[-2])*60
    if len(realtime) > 2:
        timeinsec += float(realtime[-2]) * 360
    return timeinsec

alltime = []
allcov = []
allmem = []
for errfile in glob.glob('%s/*.err'%(input_script)):
    filename = os.path.split(errfile)[-1].split('.SNP.fasta.mapper')[0]
    fasta, SNP = filename.split('.fasta.')
    outfile = errfile.replace('.err','.out')
    time_file = []
    cov_file = []
    mem_file = []
    for lines in open(errfile,'r'):
        if 'Elapsed' in lines:
            realtime = lines.split('):')[1].split('\n')[0].replace(' ','')
            time_file.append(time_sec(realtime))
        if 'Maximum' in lines:
            mem = lines.split('):')[1].split('\n')[0].replace(' ', '')
            mem_file.append(float(mem))
    if False:
        for lines in open(outfile,'r'):
            if 'queries matched' in lines:
                cov_file.append(
                    '%.2f'%(100*float(lines.split(': ')[1].split('/')[0])/float(lines.split('/')[1].split('\n')[0]))
                                )
    if len(time_file)==4:
        if time_file[0] > time_file[1]:
            print('time ', errfile)
        alltime.append('%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\n'%(fasta,SNP,time_file[0],time_file[1],time_file[2],time_file[3]))
        allmem.append(
            '%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\n' % (fasta, SNP, mem_file[0], mem_file[1], mem_file[2], mem_file[3]))
    #allcov.append('%s\t%s\t%s\t%s\n' % (fasta, SNP, cov_file[0], cov_file[1]))


f1 = open('%s/../alltimesum.mo.txt'%(input_script),'w')
f1.write('fasta\tSNP\tmapper\tbowtie2\tminimap2\tbwa\n')
f1.write(''.join(alltime))
f1.close()

f1 = open('%s/../allmemsum.mo.txt'%(input_script),'w')
f1.write('fasta\tSNP\tmapper\tbowtie2\tminimap2\tbwa\n')
f1.write(''.join(allmem))
f1.close()
if False:
    f1 = open('%s/../allcovsum.txt'%(input_script),'w')
    f1.write('fasta\tSNP\tbowtie\tmapper\n')
    f1.write(''.join(allcov))
    f1.close()
import glob
import os
import statistics
from datetime import datetime
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-ref",
                      help="path to all ref",
                      type=str, default='.',
                      metavar='ref/')
# optional input genome
optional.add_argument("-cluster",
                      help="a cluster to run, default is all clusters",
                      type=str, default='',
                      metavar='cluster1')

################################################## Definition ########################################################
args = parser.parse_args()
# set up cutoff
min_maf_for_call = .8 #Remove sample*candidate
min_cov = 6 # at least 6 reads mapped to POS

################################################### Function ########################################################
def load_indel(bowtievcf,allindel):
    ref_indel = []
    Qualified_indel = []
    indel_result = [0,0,0]# TP, same POS diff indel, diff POS
    for lines in open(bowtievcf,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR,POS,USELESS,REF,INDELset = lines_set[:5]
        INDELset = INDELset.split(',')
        Depthset = [float(i) for i in lines_set[9].split(':')[-1].split(',')]
        totaldepth = sum(Depthset)
        if totaldepth >= min_cov:
            for i in range(1,len(Depthset)):
                subdepth = Depthset[i]
                INDEL = INDELset[i-1]
                if subdepth >= min_cov and subdepth/totaldepth >= min_maf_for_call:
                    # qualified indel
                    Qualified_indel, indel_result,ref_indel = compareindelallsimple(allindel, [CHR, POS, REF,INDEL,subdepth,totaldepth], Qualified_indel,
                                                                    indel_result,ref_indel)
    f1 = open(bowtievcf + '.filtered','w')
    f1.write('CHR\tPOS\tREF\tINDEL\tDepth\tTotal_depth\tCompare\n')
    f1.write(''.join(Qualified_indel))
    f1.close()
    return [100 - len(ref_indel)] + indel_result

def compareindelsub(refstring, indelstring):
    matchingpair = 0
    indellen = len(indelstring)
    reflen = len(refstring)
    for i in range(1,indellen):
        if indelstring[(i-1):(i+1)] in refstring:
            matchingpair += 1
    if matchingpair >= min(reflen,indellen) -1 -1:
        return True
    return False

def compareindel(indelset1,indelset2):
    REF1, INDEL1 = indelset1
    REF2, INDEL2 = indelset2
    if '-' not in REF1:
        # deletion
        if len(INDEL2) < len(REF2):
            #if REF2 in REF1 or REF1 in REF2 or (len(REF1)>5 and (REF1[3:] in REF2 or REF1[:-3] in REF2)):
            if REF2 in REF1 or REF1 in REF2 or compareindelsub(REF1, REF2):
                # same deletion or a subset
                return True
        return False
    else:
        # insertion
        if len(INDEL2) > len(REF2):
            if INDEL2 in INDEL1 or INDEL1 in INDEL2 or compareindelsub(INDEL1, INDEL2):
            #if INDEL2 in INDEL1 or INDEL1 in INDEL2 or (len(INDEL1)>5 and (INDEL1[3:] in INDEL2 or INDEL1[:-3] in INDEL2)):
                # same insertion or a subset
                return True
        return False

print(compareindel(['AGCCAAATACA','A'],['AGCCAA','']))
print(compareindel(['----------','ACTGAATAAC'],['AG','AGACTGAATAAC']))
print(compareindel(['----','CCGG'],['','GCCG']))
print(compareindel(['C','CCAAGGCTAAG'],['CAAGGCTC','']))
print(compareindel(['-------------','GGTTCGGACGTGG'],['','GGTTCGGACGTGG']))

def load_refindel(refindel):
    allindel = dict()
    allindel_len = dict()
    for lines in open(refindel, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR, POS, INDEL,REF = lines_set
        POS = int(POS)
        allindel.setdefault(CHR,{})
        if '-' in REF:
            # NEED CHANGE LATER
            #POS = POS - len(REF) - 1
            indel_len =len(INDEL)
        else:
            indel_len = len(INDEL) - len(REF)
        allindel[CHR].setdefault(POS,[REF, INDEL])
        allindel_len[indel_len] = allindel_len.get(indel_len,0) + 1
    return [allindel,allindel_len]

def compareindelallsimple(allindel,indelset,Qualified_indel,indel_result,ref_indel):
    CHR, POS,REF, INDEL,subdepth,totaldepth = indelset
    POS = int(POS)
    Output = False
    if CHR in allindel:
        allindelCHR = allindel[CHR]
        for POSREF in allindelCHR:
            REFREF,INEDLREF = allindelCHR[POSREF]
            if abs(POS - POSREF) <= max(len(REFREF),len(INEDLREF)):
                # indel position within 10 bp
                Output = True
                if compareindel(allindelCHR[POSREF], [REF, INDEL]):
                    indel_result[0] += 1
                    Qualified_indel.append(
                        '%s\t%s\t%s\t%s\t%s\t%s\tTP\n' % (CHR, POS, REF, INDEL, subdepth, totaldepth))
                    if [CHR,POSREF] not in ref_indel:
                        ref_indel.append([CHR,POSREF])
                else:
                    indel_result[1] += 1
                    Qualified_indel.append(
                        '%s\t%s\t%s\t%s\t%s\t%s\tFP_samePOS\n' % (
                            CHR, POS, REF, INDEL, subdepth, totaldepth))
    if not Output:
        indel_result[2] += 1
        Qualified_indel.append(
                '%s\t%s\t%s\t%s\t%s\t%s\tFP_diffPOS\n' % (CHR, POS, REF, INDEL, subdepth, totaldepth))
    return [Qualified_indel,indel_result,ref_indel]

def compareindelall(allindel,indelset,Qualified_indel,indel_result, indel_result_set,ref_indel):
    CHR, POS,REF, INDEL,subdepth,totaldepth = indelset
    POS = int(POS)
    Output = False
    indel_len = len(INDEL) - len(REF)
    if CHR in allindel:
        allindelCHR = allindel[CHR]
        for POSREF in allindelCHR:
            if abs(POS - POSREF) <= 10:
                # indel position within 10 bp
                Output = True
                REF_REF, INDEL_REF = allindelCHR[POSREF]
                if compareindel([REF_REF, INDEL_REF], [REF, INDEL]):
                    if '-' in REF_REF:
                        indel_lenref = len(INDEL_REF)
                    else:
                        indel_lenref = len(INDEL_REF) - len(REF_REF)
                    indel_len = indel_lenref
                    indel_result[0] += 1
                    Qualified_indel.append(
                        '%s\t%s\t%s\t%s\t%s\t%s\tTP\n' % (CHR, POS, REF, INDEL, subdepth, totaldepth))
                    indel_result_set.setdefault(indel_len, [0, 0, 0])
                    indel_result_set[indel_len][0] += 1
                    if [CHR,POSREF] not in ref_indel:
                        ref_indel.append([CHR,POSREF])
                else:
                    indel_result[1] += 1
                    Qualified_indel.append(
                        '%s\t%s\t%s\t%s\t%s\t%s\tFP_samePOS\n' % (
                            CHR, POS, REF, INDEL, subdepth, totaldepth))
                    indel_result_set.setdefault(indel_len, [0, 0, 0])
                    indel_result_set[indel_len][1] += 1
    if not Output:
        indel_result[2] += 1
        Qualified_indel.append(
                '%s\t%s\t%s\t%s\t%s\t%s\tFP_diffPOS\n' % (CHR, POS, REF, INDEL, subdepth, totaldepth))
        indel_result_set.setdefault(indel_len, [0, 0, 0])
        indel_result_set[indel_len][2] += 1
    return [Qualified_indel,indel_result, indel_result_set, ref_indel]

def load_indelmapper(mappervcf, allindel):
    ref_indel = []
    Qualified_indel = []
    indel_result = [0, 0, 0]  # TP, same POS diff indel, diff POS
    indel_result_set = dict()
    insertion = ['', 0,'','',[],[]]#CHR POS REF ALT SUBdepth Totaldepth
    deletion = ['', [0],'','',[],[]]#CHR POS REF ALT SUBdepth Totaldepth
    for lines in open(mappervcf, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR, POS = lines_set[:2]
        if '-' in lines:
            POS = int(POS)
            INDELset = lines_set[3].split(',')
            # using only middle depth
            Depthset = [i for i in lines_set[5].split(';')]
            totaldepth = sum([float(j.split(',')[0]) + float(j.split(',')[1]) for j in Depthset])
            REF = lines_set[2]
            if totaldepth >= min_cov:
                for i in range(1, len(Depthset)):
                    subdepth = sum([float(j) for j in Depthset[i].split(',')])
                    INDEL = INDELset[i-1]
                    # extending indels, loose cutoff
                    if subdepth >= min_cov and subdepth / totaldepth >= min_maf_for_call - 0.1:
                        # qualified indel
                        if POS < 0:
                            #print(insertion, REF, INDEL,subdepth, totaldepth)
                            POS = -POS
                            # insertion
                            if POS != insertion[1]:
                                # diff insertion
                                # sum the previous insertion
                                if len(insertion[3]) > 1:  # insertion >= 2bp
                                    insertion[4] = statistics.mean(insertion[4])  # mean of subdepth
                                    insertion[5] = statistics.mean(insertion[5])  # mean of depth
                                    Qualified_indel, indel_result, indel_result_set,ref_indel = compareindelall(allindel,
                                                                                                      insertion,
                                                                                                      Qualified_indel,
                                                                                                      indel_result,
                                                                                                      indel_result_set,ref_indel)
                                # start of a new insertion, stringent cutoff
                                if subdepth / totaldepth >= min_maf_for_call:
                                    insertion = [CHR, POS, '', INDEL, [subdepth], [totaldepth]]
                                else:
                                    insertion = ['', 0, '', '', [], []]
                            else:
                                insertion[3] += INDEL
                                insertion[4].append(subdepth)
                                insertion[5].append(totaldepth)
                            #print(insertion)
                        elif INDEL == '-':
                            #print(deletion, CHR, POS, subdepth, totaldepth)
                            # deletion
                            if CHR != deletion[0] or abs(int(POS) - int(deletion[1][-1])) > 10:
                                # diff deletion
                                # sum the previous deletion
                                if len(deletion[2]) > 1: # deletion >= 2bp
                                    deletion[4] = statistics.mean(deletion[4]) # mean of subdepth
                                    deletion[5] = statistics.mean(deletion[5]) # mean of depth
                                    deletion[1] = deletion[1][0] # use the first POS
                                    Qualified_indel, indel_result,indel_result_set,ref_indel = compareindelall(allindel, deletion,
                                                                                    Qualified_indel, indel_result, indel_result_set,ref_indel)
                                # start of a new deletion, stringent cutoff
                                if subdepth / totaldepth >= min_maf_for_call:
                                    deletion = [CHR,[POS], REF, '', [subdepth], [totaldepth]]
                                else:
                                    deletion = ['', [0], '', '', [], []]
                            else:
                                # same deletion
                                deletion[1].append(POS)
                                deletion[2] += REF
                                deletion[4].append(subdepth)
                                deletion[5].append(totaldepth)
                            #print(deletion)
    # sum the last deletion
    if len(deletion[2]) > 1:  # deletion >= 2bp
        deletion[4] = statistics.mean(deletion[4]) # mean of subdepth
        deletion[5] = statistics.mean(deletion[5]) # mean of depth
        deletion[1] = deletion[1][0]  # use the first POS
        Qualified_indel, indel_result, indel_result_set,ref_indel = compareindelall(allindel, deletion,
                                                                          Qualified_indel, indel_result,
                                                                          indel_result_set,ref_indel)
    # sum the last insertion
    if len(insertion[3]) > 1:  # insertion >= 2bp
        insertion[4] = statistics.mean(insertion[4])  # mean of subdepth
        insertion[5] = statistics.mean(insertion[5])  # mean of depth
        Qualified_indel, indel_result, indel_result_set,ref_indel = compareindelall(allindel, insertion,
                                                                          Qualified_indel, indel_result,
                                                                          indel_result_set,ref_indel)

    f1 = open(mappervcf + '.indel.vcf.filtered', 'w')
    f1.write('CHR\tPOS\tREF\tINDEL\tDepth\tTotal_depth\tCompare\n')
    f1.write(''.join(Qualified_indel))
    f1.close()
    for indel_len in indel_result_set:
        indel_len_result = indel_result_set[indel_len]
        f2.write('%s\t%s\tmapper\t%s\t%s\t%s\t%s\n' % (indel_len,samplename, indel_len_result[0], max(allindel_len.get(indel_len,0) - indel_len_result[0],0), indel_len_result[1], indel_len_result[2]))
    return [100 - len(ref_indel)] + indel_result

##################################################### Main ##########################################################
allref = glob.glob('%s/%s*.indel.txt'%(args.ref,args.cluster))
allsum = ['Sample\tTool\tFN\tTP\tFP_samePOS\tFP_diffPOS\n']
f2 = open('%s/modelindelsumlen.txt' % (args.i), 'w')
f2.write('Indellen\tSample\tTool\tTP\tFN\tFP_samePOS\tFP_diffPOS\n')
for ref in allref:
    if '.500000.' not in ref:
        print('%s load ref indel %s' % (datetime.now(), ref))
        allindel, allindel_len = load_refindel(ref)
        # print(allindel)
        samplename = os.path.split(ref)[-1].split('.indel.txt')[0]
        if True:
            # bowtie
            bowtievcf = '%s/%s.bowtie.flt.indel.vcf' % (args.i, samplename)
            print('%s process bowtie indel %s' % (datetime.now(), bowtievcf))
            indel_result = load_indel(bowtievcf, allindel)
            allsum.append('%s\tbowtie2\t%s\t%s\t%s\t%s\n' % (samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
            # bwa
            bwavcf = '%s/%s.bwa.flt.indel.vcf' % (args.i, samplename)
            print('%s process bwa indel %s' % (datetime.now(), bwavcf))
            indel_result = load_indel(bwavcf, allindel)
            allsum.append('%s\tbwa\t%s\t%s\t%s\t%s\n' % (samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
            # minimap
            minimapvcf = '%s/%s.minimap.flt.indel.vcf' % (args.i, samplename)
            print('%s process minimap indel %s' % (datetime.now(), minimapvcf))
            indel_result = load_indel(minimapvcf, allindel)
            allsum.append('%s\tminimap2\t%s\t%s\t%s\t%s\n' % (samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))
        # mapper
        mappervcf = '%s/%s.mapper1.vcf.snp' % (args.i, samplename)
        print('%s process mapper indel %s' % (datetime.now(), mappervcf))
        indel_result = load_indelmapper(mappervcf, allindel)
        allsum.append('%s\tmapper\t%s\t%s\t%s\t%s\n' % (samplename, indel_result[0], indel_result[1], indel_result[2], indel_result[3]))

f1 = open('%s/modelindelsum.txt' % (args.i), 'w')
f1.write(''.join(allsum))
f1.close()
f2.close()

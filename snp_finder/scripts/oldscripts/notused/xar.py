# start
# use annotated genes
import glob
allanno = glob.glob('*.sum')

# load annotated genes
annotated = dict()
for lines in open('annotated.genes.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    try:
        geneset = lines_set[-2].split(';')
        for gene in geneset:
            annotated.setdefault(gene.replace(' ',''),lines.split('\n')[0])
    except IndexError:
        pass

# output annotation
confirmed = ['lolD']
for anno in allanno:
    Output = []
    for lines in open(anno,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        geneset = lines_set[6].replace(',',';').split(';')
        funset = lines_set[7].replace(',',';').split(';')
        temp_line = ''
        for gene in geneset:
            gene = gene.replace(' ','')
            if gene in annotated and gene != '':
                newline = annotated[gene]
                funnew = newline.split('\t')[-1]
                if any(fun in funnew for fun in funset) or gene in confirmed:
                    if temp_line == '':
                        temp_line = newline
                    elif temp_line!= newline:
                        if newline.split('\t')[0] != temp_line.split('\t')[0]:
                            Output.append('%s\t%s\n' % (lines.split('\n')[0], temp_line))
                            temp_line = newline
                        elif len(newline) > len(temp_line):
                            temp_line = newline
                else:
                    print(gene,funset)
                    print(funnew + '\n')
        Output.append('%s\t%s\n'%(lines.split('\n')[0],temp_line))
    f1=open(anno + '.new.txt','w')
    f1.write(''.join(Output))
    f1.close()

# add new kegg ko -> rarely used
Output = set()
temp = ['','','']
Not_output = 0
for lines in open('kegg_new','r'):
    lines_set = lines.split('\n')[0].lstrip().split(' ')
    KO = lines_set[0]
    Brite = ' '.join(lines_set)
    if not KO.startswith('K'):
        if Not_output == 0:
            temp[-1] = Brite
        elif Not_output == 1:
            temp[-2] = temp[-1]
            temp[-1] = Brite
        elif Not_output == 2:
            temp[-3] = temp[-2]
            temp[-2] = temp[-1]
            temp[-1] = Brite
        Not_output += 1
    else:
        print(temp)
        Not_output = 0
        Output.add('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(
        KO,' '.join(lines_set[1:]).split(';')[0],
        ' '.join(lines_set[1:]),'',' '.join(temp[-3].split(' ')[1:]),
        ' '.join(temp[-2].split(' ')[1:]),temp[-1].split(' [')[0]))

f1 = open('ko_format_new','w')
f1.write(''.join(list(Output)))
f1.close()
################################################### END ########################################################

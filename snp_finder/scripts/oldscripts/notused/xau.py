# start
# annitator results
annotation = dict()
for lines in open('newannotation','r'):
    lines_set = lines.split('\n')[0].split('\t')
    annotation.setdefault(lines_set[0],lines)

# merge multiple gene names
Output = []
annotation_output = []
for lines in open('all.selected.gene.faa.High_select2.faa.cluster.sum.new.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    fun = lines_set[1].replace(' ','').replace(',',';')
    fun_set = fun.split(';')
    temp_line = ''
    for fun in fun_set:
        if fun in annotation:
            annotation_output.append(fun)
            temp_line += '%s\t%s'%(lines_set[1],annotation[fun])
    if temp_line != '':
        Output.append(temp_line + '\n')

for lines in open('all.trunc.gene.faa.cluster.sum.new.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    fun = lines_set[1].replace(' ','').replace(',',';')
    fun_set = fun.split(';')
    temp_line = ''
    for fun in fun_set:
        if fun in annotation:
            annotation_output.append(fun)
            temp_line += '%s\t%s'%(lines_set[1],annotation[fun])
    if temp_line != '':
        Output.append(temp_line + '\n')

notoutput = [fun for fun in annotation if fun not in annotation_output]
print(notoutput)

f1 = open('newannotation.order.txt','w')
f1.write(''.join(list(Output)))
f1.close()

# merge multiple gene names by uniprot
annotation = dict()
Output = []
for lines in open('newannotation.order.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    geneset = lines_set[0]
    gene = lines_set[2]
    protein = lines_set[3]
    geneprotein = '%s\t%s'%(gene,protein)
    if gene!='None':
        if geneprotein not in annotation:
            annotation.setdefault(geneprotein,lines.split('\n')[0])
        else:
            if annotation[geneprotein].split('\t')[0] != geneset:
                if lines_set[1] != annotation[geneprotein].split('\t')[1]:
                    annotation[geneprotein] = geneset + ';' + annotation[geneprotein] + '\t' + lines.split('\n')[0]
                else:
                    annotation[geneprotein] = geneset + ';' + annotation[geneprotein]
    else:
        Output.append('\t%s\n' % (lines.split('\n')[0]))

Output = []
for geneprotein in annotation:
    Output.append('%s\t%s\n'%(geneprotein,annotation[geneprotein]))

f1 = open('newannotation.order.txt','w')
f1.write(''.join(list(Output)))
f1.close()

# delete annotated genes
Output = []
annotated = dict()
notright = ['lacI','abfA','asa1','recG','tetM','cbpA']
for lines in open('annotated.genes.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    try:
        geneset = lines_set[-2].replace(' ','').replace(',',';').split(';')
        for gene in geneset:
            annotated.setdefault(gene.replace(' ',''),lines.split('\n')[0])
    except IndexError:
        pass

for lines in open('newannotation.order.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    gene = lines_set[2]
    geneset = gene.replace(' ','').replace(',',';').split(';')
    alreadyin = False
    for gene in geneset:
        if gene in annotated and gene not in notright:
            print(gene, annotated[gene])
            print(geneset, lines)
            alreadyin = True
            break
    if not alreadyin:
        Output.append(lines)

f1 = open('newannotation.order.2.txt','w')
f1.write(''.join(list(Output)))
f1.close()

#>> annotate and update annotated.genes.txt
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
confirmed = ['lolD','FAEPRAA2165_RS06535','FAEPRAA2165_RS06280','crtI','cbpA','valS','pepO']
Output = []
for anno in allanno:
    if 'High_select' in anno:
        tag = 'High_select'
    elif 'trunc' in anno:
        tag = 'Trunc'
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
                            Output.append('%s\t%s\t%s\n' % ('\t'.join(lines_set[0:8]), tag, temp_line))
                            temp_line = newline
                        elif len(newline) > len(temp_line):
                            temp_line = newline
                else:
                    print(gene,funset)
                    print(funnew + '\n')
        Output.append('%s\t%s\t%s\n'%('\t'.join(lines_set[0:8]),tag,temp_line))

f1=open(anno + '.all.new.txt','w')
f1.write(''.join(Output))
f1.close()

################################################### END ########################################################

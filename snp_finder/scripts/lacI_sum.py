import os,glob
output_folder = '/scratch/users/anniz44/genomes/HTH/diamond/'
alloutput = glob.glob('%s/*.txt'%(output_folder))
metadata = '/scratch/users/anniz44/genomes/HTH/GMC_BN10_WGS_metadata.all.update.txt'
identity_cutoff = 50

def load_taxameta(metadata):
    taxameta = dict()
    for lines in open(metadata,'r'):
        lines_set = lines.split('\t')
        strainname = lines_set[1]
        species = lines_set[2].split('_')
        if species[0].startswith('Bifido'):
            species = '%s%s'%(species[0][:2],species[1][:2])
            species = species.replace('Biad','BA').replace('Bilo','BL').replace('Bips','BiPs')
            env = '\t'.join(lines_set[12:15])
            taxameta.setdefault(strainname,[species,env + '\tH'])
    return taxameta

def compare_geneanno(geneanno1,geneanno2):
    if geneanno1[1] < geneanno2[1] and geneanno1[2]*0.75 < geneanno2[2]:
        # new gene is better
        return geneanno2
    else:
        return geneanno1

def load_diamondresult(diamondresult,donorspeciessum,donorspeciesanno,taxameta):
    strainname = os.path.split(diamondresult)[-1].split('.txt')[0]
    donor = strainname.split('_')[0]
    if strainname not in taxameta:
        Jaydata = True
        species = strainname.split('_')[1]
        env = 'urban\tindustrialized\twestern\tH'
        if strainname.startswith('P'):
            env = 'urban\tindustrialized\twestern\tP'
    else:
        Jaydata = False
        species,env = taxameta.get(strainname,['',''])
    if donor == 'D131' or (donor == 'D116' and species !='BA') or donor == 'D113' or (donor == 'D114'and species !='BA') :
        species = 'BL'
    donorspecies = '%s_%s'%(donor,species)
    donorspeciesanno.setdefault(donorspecies,'%s'%(env))
    geneanno = dict()
    for lines in open(diamondresult,'r'):
        lines_set = lines.split('\t')
        if not Jaydata:
            genename = lines_set[0].split('prokka.faa_')[1]
        else:
            genename = lines_set[0]
        lacIname = lines_set[1].split('_')[-1]
        if 'unknown' in lines_set[1]:
            lacIname = 'unknown_LacI'
        identity = float(lines_set[2])
        coverage = float(lines_set[3])
        if identity >= identity_cutoff:
            if genename not in geneanno:
                geneanno.setdefault(genename,[lacIname,identity,coverage])
            else:
                geneanno[genename] = compare_geneanno(geneanno[genename],[lacIname,identity,coverage])
    strain_count,geneannosum = donorspeciessum.get(donorspecies,[0,dict()])
    strainsumdict = dict()
    for genename in geneanno:
        lacIname, identity, coverage = geneanno[genename]
        geneannosum[lacIname] = geneannosum.get(lacIname,0) + 1
        strainsumdict.setdefault(lacIname,[0,identity])
        strainsumdict[lacIname][0] += 1
        strainsumdict[lacIname][1] = min(strainsumdict[lacIname][1],identity)
    for lacIname in strainsumdict:
        strainsum.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (strainname,donorspecies, species, lacIname,strainsumdict[lacIname][0],strainsumdict[lacIname][1]))
    strain_count += 1
    donorspeciessum[donorspecies] = [strain_count,geneannosum]
    return [donorspeciessum,donorspeciesanno]

# load taxa meta
taxameta = load_taxameta(metadata)
# process diamond result
strainsum = []
donorspeciessum = dict()
donorspeciesanno = dict()
for diamondresult in alloutput:
    donorspeciessum,donorspeciesanno = load_diamondresult(diamondresult,donorspeciessum,donorspeciesanno,taxameta)
allsum = []
for donorspecies in donorspeciessum:
    env = donorspeciesanno[donorspecies]
    strain_count,geneannosum = donorspeciessum[donorspecies]
    for lacIname in geneannosum:
        allsum.append('%s\t%s\t%s\t%s\t%s\n'%(donorspecies,lacIname,geneannosum[lacIname]/strain_count,strain_count,env))
# output diamond result
f1 = open('%s/../alldiamond.txt'%(output_folder),'w')
f1.write('donor_species\tLacI\tunique_genenum\tstrainnum\turbanism\tlifestyle\tlifestyle2\thealth\n')
f1.write(''.join(allsum))
f1.close()
f1 = open('%s/../alldiamondstrainsum.txt'%(output_folder),'w')
f1.write('strainname\tdonor_species\tspecies\tLacI\tunique_genenum\tminID\n')
f1.write(''.join(strainsum))
f1.close()

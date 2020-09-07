import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import itertools
import random
# set up path
Round = 4
Cluster = True
Tree = True
Paircompare = False
Cov_dis = 20

#input_script_temp = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s_allpair'%(Round)
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s_tree'%(Round)
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
input_script2 = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round%s'%(Round)
genome_root2 = '/scratch/users/anniz44/genomes/donor_species/jay'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round%s/*'%(Round))
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/bwa/0/'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge'%(Round)
vcf_name = '.all.flt.snp.vcf'
ref_filename = '.all.spades%s.fasta'%(Round)
fasta_name = '.fasta.corrected.fasta'
deleting_file = []
fastq_name = '_1.fastq'
Species_replace = dict()
Species_replace.setdefault('BA','Bifidobacterium_adolescentis')
Species_replace.setdefault('BL','Bifidobacterium_longum')
if Round == 4:
    genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*'
    output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge_genome/vcf' % (
        Round)

#input_script_temp = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s_allpair'%(Round)
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s_tree'%(Round)
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
genome_root = '/scratch/users/anniz44/genomes/donor_species/jay/round%s'%(Round)
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round%s/*'%(Round))
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/bwa/0/'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge'%(Round)
vcf_name = '.all.flt.snp.vcf'
ref_filename = '.all.spades%s.fasta'%(Round)
fasta_name = '.fasta.corrected.fasta'
fastq_name = '.sorted.bam'
deleting_file = []
if Round == 4:
    genome_root = '/scratch/users/anniz44/genomes/donor_species/jay/round*'
    output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge_genome/vcf' % (Round)

# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def vcf_to_txt(lines,output_list,cluster_sub=[]):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            if (cluster_sub==[] and i not in deleting_set) or i in cluster_sub:
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputcov(output_name,vcf_file_POS_candidate,cluster_sub=[]):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        if '.fna.flt.snp.vcf' in vcf_file:
            cov_file = vcf_file.split('.fna.flt.snp.vcf')[0] + '.fna.raw.vcf'
        else:
            cov_file = vcf_file.split('.flt.snp.vcf')[0] + '.raw.vcf'
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            cov_file,
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name),
                                       os.path.join(input_script, 'grep.temp.txt')))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()

def outputtree(output_name):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            SNP_alignment_output_parsi.append('S%s    %s\n' % (genomename[-8:], SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    if Tree:
        vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' % (output_name), 'w')
        vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
        vcf_file_filtered.close()

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

def SNP_distance_correct(distanceNEW, distanceREF, SNPREF):
    return int((distanceNEW/distanceREF)*SNPREF)

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    if all_genes == []:
        # database is a gene database
        codon_start = POS - 1 - int((POS - 1) % 3)
        Ref_seq_chr = Ref_seq.get(CHR, 'None')
        return [CHR, POS, codon_start, Ref_seq_chr, Reverse_chr]
    else:
        # database is a contig database
        for a_gene in all_genes:
            POS1, POS2, GENE = a_gene
            if POS >= POS1 and POS <= POS2:
                Ref_seq_chr = Ref_seq.get(GENE, 'None')
                Gene_length = len(Ref_seq_chr)
                if GENE in Reverse:  # reversed
                    POS_gene = Gene_length-(int(POS-POS1))
                    Reverse_chr = 1
                else:
                    POS_gene = int(POS-POS1)+1
                codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
                return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return True
    else:
        return False

def depthcheck(vcf_genome,vcf_fq):
    os.system('cat %s | cut -f 1,2 > %s.temp' %(vcf_genome,vcf_genome))
    os.system('grep -f %s.temp %s --no-group-separator > %s.temp.depth' % (
        vcf_genome,
        vcf_fq,
        vcf_genome))
    Length = dict()
    Depth_set = dict()
    Total = 0
    for lines in open(vcf_genome + '.temp.depth'):
        lines_set = lines.split('\n')[0].split('\t')
        if Total == 0:
            Total = len(lines_set) - 9
        CHR = lines_set[0]
        if CHR not in Length:
            try:
                total_length = CHR.split('size')[1]
            except IndexError:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            total_length = int(total_length)
            Length.setdefault(CHR, total_length)
        total_length = Length[CHR]
        if total_length >= Length_cutoff:
            POS = lines_set[1]
            CHRPOS = '%s\t%s' % (CHR, POS)
            lines_set_sub = lines_set[9:]
            for i in range(0,Total):
                Subdepth_all = lines_set_sub[i]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                if total_sub_depth >= Depth_cutoff:
                    Depth_set.setdefault(CHRPOS, [])
                    Depth_set[CHRPOS].append(i)
    os.system('rm -rf %s.temp*'%(vcf_genome))
    return Depth_set

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff,Depth_set,cluster_sub=[]):
    CHR = lines_set[0]
    POS = int(lines_set[1])
    CHRPOS = '%s\t%s'%(CHR,POS)
    if Depth_set == {} or CHRPOS in Depth_set:
        temp_snp_line = []
        temp_snp_line_frq = []
        temp_snp_line_NS = ['None', 'None', 'None']
        temp_snp_line_AA = ''
        Total_qualify = 0
        Total_qualify_SNP = 0
        Total_qualify_notSNP = 0
        Total_unqualify_alt_freq = 0
        SNP = set()
        SNP_seq = []
        REF = lines_set[3]
        allels_set = [REF]
        Total_subsample = Total
        lines_set_sub = lines_set[9:]
        REF_where=0
        if cluster_sub!= []:
            lines_set_sub = [lines_set[i] for i in cluster_sub]
            Total_subsample = len(cluster_sub)
            if Total_subsample >= 15:
                SNP_presence_cutoff = 0.33  # for a large group of samples
            elif Total_subsample in [3,4]:
                SNP_presence_cutoff = 1  # for a small group of samples
                SNP_presence_sample_cutoff = 2
            elif Total_subsample in [1,2]:
                SNP_presence_cutoff = 1  # for a small group of samples
                SNP_presence_sample_cutoff = 1
                no_SNP_cutoff = 0
        else:
            cluster_sub = list(range(9, len(lines_set)))
        if Total_subsample > 0:
            if '.' not in lines_set[4]:
                allels_set += lines_set[4].split(',')
            Total_alleles = len(allels_set)
            genome_order = 0
            Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
            if Total_subsample > 2:
                REF,REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
            sample_num = 9
            for Subdepth_all in lines_set_sub:
                if sample_num not in deleting_set:
                    genome_order += 1
                    Allels_frq = [0, 0, 0, 0]
                    Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                    total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                    for num_allels in range(0, Total_alleles):
                        allels = allels_set[num_allels]
                        Subdepth_alleles = int(Subdepth[num_allels])
                        if allels in Allels:
                            Allels_frq[Allels[allels]] += Subdepth_alleles
                        else:
                            pass
                    # find major alt and calculate frequency
                    Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                    temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq))
                    SNP_seq.append(REF)  # set as reference
                    if total_sub_depth > 0:
                        MLF = Major_ALT[1] / total_sub_depth
                        if MLF >= Major_alt_freq_cutoff:
                            # major alt frequency cutoff
                            Total_qualify += 1
                            # check for qualified SNP
                            if Major_ALT[0] != REF:
                                if (Depth_set == {} or sample_num - 9 in Depth_set[CHRPOS]): # depth cutoff
                                    # a qualified SNP
                                    temp_snp_line_pass += 'PASS'
                                    Total_qualify_SNP += 1
                                    SNP.add(genome_order)  # only take qualified SNP as valid SNP
                                    SNP_seq[-1] = Major_ALT[0] # only qualified SNP include in alignment
                                else:
                                    print('deleting sample %s of CHRPOS %s' % (sample_num - 9, CHRPOS))
                            else:
                                Total_qualify_notSNP += 1
                        else:
                            # major alt frequency low
                            Total_unqualify_alt_freq += 1
                sample_num += 1
            if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
                    Total_unqualify_alt_freq < Poor_MLF_freq_cutoff and\
                    Total_qualify >= SNP_presence_sample_cutoff and \
                    Total_qualify_SNP >= 1 and Total_qualify_SNP <= Total_qualify - no_SNP_cutoff and\
                    Total_qualify_notSNP >= no_SNP_cutoff:
                # -> qualified SNP, qualified samples cutoff + unqualified samples cutoff, at least 1 qualified SNP, calculate NS
                gene_info = contig_to_gene(CHR, POS)
                if gene_info!= []:
                    Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                    if Ref_seq_chr != 'None':
                        #  observed NS ratio calculated
                        temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                        if codon_start <= POS_gene - 1:
                            Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                            Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_chr = Ref_seq_chr
                            if len(Ref_seq_codon) == 3:
                                Ref_seq_aa = translate(Ref_seq_codon)[0]
                                temp_snp_line_AA += Ref_seq_aa
                                ALT_set = allels_set
                                ALT_set.remove(REF)
                                for ALT in ALT_set:
                                    SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                    SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                    SNP_seq_aa = translate(SNP_seq_codon)[0]
                                    temp_snp_line_AA += SNP_seq_aa
                                    temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                    temp_snp_line_NS[-1]+=temp_NorS
                # output lines and output major alt
                if 'PASS' in temp_snp_line_pass:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'NOPASS'
                if CHR == CHR_old:
                    # same CHR
                    POS_DIS = abs(POS - POS_old)
                    vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
                else:
                    # diff CHR first SNP
                    vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
                POS_old = POS
                CHR_old = CHR
                temp_snp_line.append(CHR)
                temp_snp_line.append(str(POS))
                temp_snp_line.append(REF)
                temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
                vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                    ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
                vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
                vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
                i = 9
                j = 0
                SNP_alignment[reference_name] += REF
                for genomename in SNP_alignment:
                    if genomename != reference_name:
                        if i in cluster_sub:
                            SNP_alignment[genomename] += SNP_seq[j]
                            j += 1
                        i += 1
    else:
        print('deleting CHRPOS %s'%(CHRPOS))
    return [CHR_old,POS_old]

def SNP_check_output(lines_set,CHR_old,POS_old,reference_name):
    CHR = lines_set[0]
    POS = int(lines_set[1])
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    REF_where = 0
    cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        if Total_subsample > 2:
            REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    qualify_loci = 0
                    MLF = Major_ALT[1] / total_sub_depth
                    if Major_ALT[0] != REF and MLF >= Major_alt_freq_cutoff or qualify_loci == 1:
                        # check for qualified SNP
                        SNP.add(genome_order)  # only take qualified SNP as valid SNP
                        SNP_seq[-1] = Major_ALT[0]  # only qualified SNP include in alignment
            sample_num += 1
        gene_info = contig_to_gene(CHR, POS)
        if gene_info != []:
            Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
            if Ref_seq_chr != 'None':
                #  observed NS ratio calculated
                temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                if codon_start <= POS_gene - 1:
                    Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                    Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                    SNP_seq_chr = Ref_seq_chr
                    if len(Ref_seq_codon) == 3:
                        Ref_seq_aa = translate(Ref_seq_codon)[0]
                        temp_snp_line_AA += Ref_seq_aa
                        ALT_set = allels_set
                        ALT_set.remove(REF)
                        for ALT in ALT_set:
                            SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                            SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_aa = translate(SNP_seq_codon)[0]
                            temp_snp_line_AA += SNP_seq_aa
                            temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                            temp_snp_line_NS[-1] += temp_NorS
        temp_snp_line_pass = 'PASS'
        if CHR == CHR_old:
            # same CHR
            POS_DIS = abs(POS - POS_old)
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
        else:
            # diff CHR first SNP
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
        POS_old = POS
        CHR_old = CHR
        temp_snp_line.append(CHR)
        temp_snp_line.append(str(POS))
        temp_snp_line.append(REF)
        temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
        vcf_file_list.append(
            '\t'.join(temp_snp_line) + '\t' + '\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,
                '\t'.join(temp_snp_line_NS), temp_snp_line_AA))
        vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
        vcf_file_list_vcf.append('\t'.join(lines_set[0:9]) + '\t' + '\t'.join(lines_set_sub) + '\n')
        i = 9
        j = 0
        SNP_alignment[reference_name] += REF
        for genomename in SNP_alignment:
            if genomename != reference_name:
                if i in cluster_sub:
                    SNP_alignment[genomename] += SNP_seq[j]
                    j += 1
                i += 1
    return [CHR_old,POS_old]

def dis_pos(current_CHRPOS,pre_CHRPOS):
    POS1 = int(pre_CHRPOS.split('\t')[1])
    POS2 = int(current_CHRPOS.split('\t')[1])
    return POS2 - POS1

def compare_set(Ge_pre,Ge_cu):
    Ge_preset = Ge_pre.replace('\"','').split(';')
    Ge_cuset = Ge_cu.replace('\"', '').split(';')
    if (len(Ge_preset) >= 3 and len(Ge_cuset) >= 2) or (len(Ge_preset) >= 2 and len(Ge_cuset) >= 3):
        return all(elem in Ge_cuset for elem in Ge_preset) or all(elem in Ge_preset for elem in Ge_cuset)
    else:
        return False

def N_ratio(allrec,SNP_N_set):
    N = 0
    for CHRPOS in allrec:
        if SNP_N_set[CHRPOS] == 'N':
            N += 0
    return N / len(allrec)

def cluster_rec(CHR_set,SNP_gen_set,SNP_N_set,CHRPOS_set):
    for CHR in CHR_set:
        potential_rec = dict()
        allCHRPOS = CHR_set[CHR]
        for i in range(1,len(allCHRPOS)):
            current_CHRPOS = allCHRPOS[i]
            for j in reversed(range(0, i)):
                pre_CHRPOS = allCHRPOS[j]
                Ge_pre = SNP_gen_set[pre_CHRPOS]
                Ge_cu = SNP_gen_set[current_CHRPOS]
                Dis = dis_pos(current_CHRPOS,pre_CHRPOS)
                if Dis <= Rec_length_cutoff:
                    if Ge_pre == Ge_cu or compare_set(Ge_pre,Ge_cu):
                        # cluster rec sites
                        potential_rec.setdefault(Ge_pre, set())
                        potential_rec[Ge_pre].add(current_CHRPOS)
                        potential_rec[Ge_pre].add(pre_CHRPOS)
                        potential_rec.setdefault(Ge_cu, set())
                        potential_rec[Ge_cu].add(current_CHRPOS)
                        potential_rec[Ge_cu].add(pre_CHRPOS)
                        potential_rec[Ge_cu].update(list(potential_rec[Ge_pre]))
                        potential_rec[Ge_pre].update(list(potential_rec[Ge_cu]))
                        break
                else:
                    if j == i-1:
                        # output recombination
                        delete_set = []
                        for Ge_pre in potential_rec:
                            allrec = potential_rec[Ge_pre]
                            if len(allrec) >= Rec_SNP_cutoff:
                                CHRPOS_set += allrec
                            delete_set.append(Ge_pre)
                        for Ge_pre in delete_set:
                            potential_rec.pop(Ge_pre, 'None')
                    break
        # output the last recombination
        for Ge_pre in potential_rec:
            allrec = potential_rec[Ge_pre]
            if len(allrec) >= Rec_SNP_cutoff:
                CHRPOS_set += allrec
    return CHRPOS_set

def remove_rec(SNP_file):
    CHRPOS_set = []
    CHR_set = dict()
    SNP_gen_set = dict()
    SNP_N_set = dict()
    # import SNP info
    for lines in open(SNP_file,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR = lines_set[0]
        POS = lines_set[1]
        CHRPOS = '%s\t%s' % (CHR, POS)
        Ge = lines_set[-6]
        NS = lines_set[-2]
        if 'NN' in NS:
            NS = 'N'
        SNP_gen_set.setdefault(CHRPOS,Ge)
        SNP_N_set.setdefault(CHRPOS, NS)
        CHR_set.setdefault(CHR, [])
        CHR_set[CHR].append(CHRPOS)
    CHRPOS_set = cluster_rec(CHR_set,SNP_gen_set,SNP_N_set,CHRPOS_set)
    return CHRPOS_set

def load_ref_vcf(ref_vcf_file):
    ref_chr = dict()
    for files in ref_vcf_file:
        Set_length = False
        for lines in open(files,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, Notused, REF, ALT = lines_set[0:5]
            CHR_POS = '%s__%s'%(CHR, POS)
            ref_chr.setdefault(CHR_POS,[])
            ref_chr[CHR_POS]=[REF,ALT]
    return ref_chr

# set up output
if Tree:
    try:
        os.mkdir(output_dir_merge + '/tree')
    except IOError:
        pass
    try:
        os.mkdir(input_script_sub)
    except IOError:
        pass

# set up cutoff
reference_set = ['reference']
outputname_set = ['filtered']
SNP_presence_cutoff = 0.66  # avg presence in all samples
SNP_presence_sample_cutoff = 3  # num of samples passing the above criteria
Major_alt_freq_cutoff = 0.9 # major alt freq in a genome, do not allow multiple homolougs genes
no_SNP_cutoff = 1
Poor_MLF_freq_cutoff = 1 # no sample should have homologous genes (low major alt freq)
# set up strict cutoff
SNP_presence_cutoff2 = 0.66 # avg coverage in all samples
SNP_presence_sample_cutoff2 = 3  # num of samples passing the above criteria
Poor_MLF_freq_cutoff2 = 1 # no sample should have homologous genes (low major alt freq)
# cluster cutoff Step 3
SNP_total_cutoff_2 = 100
cluster_cutoff = 2
# Depth and recombination cutoff Round4
Depth_cutoff = 10 # covered by 10 reads
Length_cutoff = 2000 # minimum ref contig length
Rec_length_cutoff = 1000 # maximum distance between recombination sites
Rec_SNP_cutoff = 4 # minumum no. of SNPs grouped/clustered as a recombination
end_cutoff = 70 # contig end no SNP calling

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

# run vcf filtering
if Round == 4:
    vcf_name = '.all.flt.snp.vcf.filtered.vcf'
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    outputname_set = ['final']
    ref_filename = '.all.spades*.fasta'
else:
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))

reference_name = reference_set[0]
output_name = outputname_set[0]
for vcf_file in all_vcf_file:
    filesize = 0
    try:
        filesize = int(os.path.getsize(vcf_file + '.%s.vcf'%(output_name)))
    except FileNotFoundError:
        pass
    if filesize == 0:
        SNP_presence_cutoff = SNP_presence_cutoff2  # for group of samples
        SNP_presence_sample_cutoff = SNP_presence_sample_cutoff2
        no_SNP_cutoff = 1
        print(vcf_file)
        Total = 0
        # filter depth
        Depth_set = dict()
        ref_chr = dict()
        if Round == 4:
            vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
            # need change later
            vcf_file_raw = \
            glob.glob(output_dir_merge + '/../../../vcf_round*/merge/' + vcf_ref_file_name + '.all.raw.vcf')[0]
            try:
                # WGS
                # need change later
                vcf_fq = glob.glob(output_dir_merge + '/../../merge/' + vcf_ref_file_name + '*.all.fq.flt.snp.vcf')[0]
                #Depth_set = depthcheck(vcf_file, vcf_fq)
                ref_chr = load_ref_vcf([vcf_fq])
            except IndexError:
                pass
        else:
            vcf_file_raw = vcf_file.replace('.flt.snp.vcf', '.raw.vcf')
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
        database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
        if len(database) > 1:
            print(vcf_file,database)
        database = database[0]
        ref_dir, ref_name = os.path.split(database)
        ref_fna = database.replace('.fasta', '.fna')
        try:
            f1 = open(ref_fna, 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -d %s' % (database, ref_fna))
        Sample_name = []
        deleting_set = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        for lines in open(vcf_file_raw, 'r'):
            if lines.startswith('##bcftoolsCommand=mpileup '):
                # setup samples
                sample_set = lines.split(ref_name + ' ')[1].split('\n')[0].split(' |')[0].split(' ')
                samplenum = 9
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                    Sample_name.append(genomename.replace('.', ''))
                    if genomename in deleting_file:
                        deleting_set.append(samplenum)
                    samplenum += 1
                break
        print('running %s' % (donor_species))
        # load database
        database_file = ref_fna
        Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
        SNP_tree_cmd = []
        SNP_tree_cmd2 = []
        vcf_file_list = []
        vcf_file_list_vcf = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        SNP_alignment = dict()
        SNP_alignment.setdefault(reference_name, '')
        cov_file_list = []
        CHR_old = ''
        POS_old = 0
        for genomename in Sample_name:
            SNP_alignment.setdefault(genomename, '')
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\n')[0].split('\t')
                CHR = lines_set[0]
                POS = int(lines_set[1])
                # a SNP confirmed in WGS mapping
                CHR_POS = '%s__%s' % (CHR, POS)
                if (Round < 4 or CHR_POS in ref_chr) and not contig_end(CHR, POS):
                    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                    if Total == 0:
                        Total = len(lines_set) - 9 - len(deleting_set)
                        if Total >= 15:
                            SNP_presence_cutoff = 0.33  # for a large group of genomes
                        elif Total in [3, 4]:
                            SNP_presence_cutoff = 1  # for a small group of genomes
                            SNP_presence_sample_cutoff = 2
                        elif Total in [1, 2]:
                            SNP_presence_cutoff = 1  # for only 1 or 2 samples, compare to ref
                            SNP_presence_sample_cutoff = 1
                            no_SNP_cutoff = 0
                    if Depth / Total >= SNP_presence_cutoff:
                        # average depth in all samples cutoff
                        if "INDEL" not in lines_set[7] \
                                and (lines_set[6] != 'LowQual'):
                            CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                             CHR_old, POS_old, reference_name,
                                                             SNP_presence_cutoff,
                                                             SNP_presence_sample_cutoff, no_SNP_cutoff, Depth_set)
        outputvcf(output_name)
        outputtree(output_name)

# remove recombination
if Round == 4:
    vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf'
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    print(all_vcf_file)
    outputname_set = ['removerec']
    ref_filename = '.all.spades*.fasta'
    reference_name = reference_set[0]
    output_name = outputname_set[0]
    for vcf_file in all_vcf_file:
        filesize = 0
        try:
            filesize = int(os.path.getsize(vcf_file + '.%s.vcf' % (output_name)))
        except FileNotFoundError:
            pass
        if filesize == 0:
            print(vcf_file)
            Total = 0
            # filter recombination
            vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
            # need change later
            vcf_file_raw = \
                glob.glob(output_dir_merge + '/../../../vcf_round*/merge/' + vcf_ref_file_name + '.all.raw.vcf')[0]
            SNP_file = vcf_file.replace('.final.vcf', '.final.snp.txt')
            CHRPOS_set = remove_rec(SNP_file)
            donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
            database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
            if len(database) > 1:
                print(vcf_file, database)
            database = database[0]
            ref_dir, ref_name = os.path.split(database)
            ref_fna = database.replace('.fasta', '.fna')
            try:
                f1 = open(ref_fna, 'r')
            except FileNotFoundError:
                os.system('prodigal -q -i %s -d %s' % (database, ref_fna))
            Sample_name = []
            deleting_set = []
            Ref_seq = dict()
            Mapping = dict()
            Mapping_loci = dict()
            for lines in open(vcf_file_raw, 'r'):
                if lines.startswith('##bcftoolsCommand=mpileup '):
                    # setup samples
                    sample_set = lines.split(ref_name + ' ')[1].split('\n')[0].split(' |')[0].split(' ')
                    samplenum = 9
                    for samples in sample_set:
                        genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                        Sample_name.append(genomename.replace('.', ''))
                        if genomename in deleting_file:
                            deleting_set.append(samplenum)
                        samplenum += 1
                    break
            print('running %s' % (donor_species))
            # load database
            database_file = ref_fna
            Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
            SNP_tree_cmd = []
            SNP_tree_cmd2 = []
            vcf_file_list = []
            vcf_file_list_vcf = []
            vcf_file_POS = []
            vcf_file_POS_candidate = set()
            SNP_alignment = dict()
            SNP_alignment.setdefault(reference_name, '')
            cov_file_list = []
            CHR_old = ''
            POS_old = 0
            for genomename in Sample_name:
                SNP_alignment.setdefault(genomename, '')
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    CHRPOS = '%s\t%s' % (CHR, POS)
                    if Total == 0:
                        Total = len(lines_set) - 9 - len(deleting_set)
                    if CHRPOS not in CHRPOS_set and not contig_end(CHR, POS):
                        CHR_old, POS_old = SNP_check_output(lines_set,
                                                            CHR_old, POS_old, reference_name)
            outputvcf(output_name)
            outputtree(output_name)

# run parsi tree
if Tree:
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    for vcf_file in all_vcf_file:
        a_parsi_file = vcf_file + '.%s.parsi.fasta'%(output_name)
        if 'RAxML_parsimonyTree' not in a_parsi_file:
            os.system('rm -rf %s %s' % (a_parsi_file + '.out.txt',
                                        a_parsi_file + '.out.tree'))
            SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (a_parsi_file))
            f1 = open(os.path.join(input_script_sub, 'parsi.optionfile.txt'), 'w')
            f1.write(SNP_tree_cmd3)
            f1.close()
            os.system('rm -rf outfile outtree')
            os.system('dnapars < %s/parsi.optionfile.txt > %s/%s.parsi.output\n' % (
                input_script_sub, input_script_sub, os.path.split(a_parsi_file)[-1]))
            os.system('mv outfile %s' % (a_parsi_file + '.out.txt'))
            os.system('mv outtree %s' % (a_parsi_file + '.out.tree'))
    os.system('mv %s/*.parsi* %s/tree' % (
        output_dir_merge, output_dir_merge))

################################################### END ########################################################
################################################### SET PATH ########################################################
# After round 4 filter results of WGS
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import itertools
import random
# set up path
Cov_dis = 20
Round = 4

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s_tree'%(Round)
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round%s'%(Round)
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round%s/*'%(Round))
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/bwa/0/'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge'%(Round)
vcf_name = '.all.flt.snp.vcf'
ref_filename = '.all.spades%s.fasta'%(Round)
fasta_name = '.fasta.corrected.fasta'
deleting_file = []
fastq_name = '_1.fastq'
if Round == 4:
    genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*'
    output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge/' % (
        Round)

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s_tree'%(Round)
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
genome_root = '/scratch/users/anniz44/genomes/donor_species/jay/round%s'%(Round)
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round%s/*'%(Round))
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/bwa/0/'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge/'%(Round)
vcf_name = '.all.flt.snp.vcf'
ref_filename = '.all.spades%s.fasta'%(Round)
fasta_name = '.fasta.corrected.fasta'
fastq_name = '.sorted.bam'
deleting_file = []
if Round == 4:
    genome_root = '/scratch/users/anniz44/genomes/donor_species/jay/round*'
    output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge/' % (Round)

################################################### Function ########################################################
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def vcf_to_txt(lines,output_list,cluster_sub=[]):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            if (cluster_sub==[] and i not in deleting_set) or i in cluster_sub:
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputcov(output_name,vcf_file_POS_candidate,cluster_sub=[]):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            vcf_file.split('.flt.snp.vcf')[0] + '.raw.vcf',
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name)))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()

def outputtree(output_name):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            SNP_alignment_output_parsi.append('S%s    %s\n' % (genomename[-8:], SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' %(output_name), 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

def SNP_distance_correct(distanceNEW, distanceREF, SNPREF):
    return int((distanceNEW/distanceREF)*SNPREF)

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    for a_gene in all_genes:
        POS1, POS2, GENE = a_gene
        if POS >= POS1 and POS <= POS2:
            Ref_seq_chr = Ref_seq.get(GENE, 'None')
            Gene_length = len(Ref_seq_chr)
            if GENE in Reverse:  # reversed
                POS_gene = Gene_length-(int(POS-POS1))
                Reverse_chr = 1
            else:
                POS_gene = int(POS-POS1)+1
            codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
            return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,cluster_sub=[]):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    Total_qualify = 0
    Total_qualify_SNP = 0
    Total_qualify_notSNP = 0
    Total_unqualify_alt_freq = 0
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    if cluster_sub!= []:
        lines_set_sub = [lines_set[i] for i in cluster_sub]
        Total_subsample = len(cluster_sub)
        if Total_subsample >= 15:
            SNP_presence_cutoff = 0.33  # for a large group of samples
        if Total_subsample <= 3:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 2
    else:
        cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
                total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    MLF = Major_ALT[1] / total_sub_depth
                    if total_sub_depth_forward >= Sample_depth_cutoff and \
                            total_sub_depth_reverse >= Sample_depth_cutoff:
                        # forward and reverse cutoff POS detected
                        if MLF >= Major_alt_freq_cutoff:
                        # major alt frequency cutoff
                            Total_qualify += 1
                            # check for qualified SNP
                            if Major_ALT[0] != REF:
                                SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                                # qualified SNP
                                if (total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= SNP_depth_cutoff and \
                                        total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= SNP_depth_cutoff) and \
                                        Major_ALT[1] >= SNP_depth_cutoff:
                                    Total_qualify_SNP += 1
                                    SNP.add(genome_order)  # only take qualified SNP as valid SNP
                                    SNP_seq[-1] = Major_ALT[0]
                            else:
                                Total_qualify_notSNP += 1
                        else:
                            # major alt frequency low
                            Total_unqualify_alt_freq += 1
            sample_num += 1
        if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
                Total_unqualify_alt_freq / Total_subsample <= Poor_MLF_freq_cutoff and\
                Total_qualify >= SNP_presence_sample_cutoff and \
                Total_qualify_SNP >= 1 and Total_qualify_SNP < Total_qualify and\
                Total_qualify_notSNP > 0:
            # -> qualified SNP
            # qualified samples cutoff + unqualified samples cutoff
            # at least 1 qualified SNP
            # calculate NS
            gene_info = contig_to_gene(CHR, POS)
            if gene_info!= []:
                Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                if Ref_seq_chr != 'None':
                    #  observed NS ratio calculated
                    temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                    if codon_start <= POS_gene - 1:
                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_chr = Ref_seq_chr
                        if len(Ref_seq_codon) == 3:
                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                            temp_snp_line_AA += Ref_seq_aa
                            ALT_set = allels_set
                            ALT_set.remove(REF)
                            for ALT in ALT_set:
                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                temp_snp_line_AA += SNP_seq_aa
                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                temp_snp_line_NS[-1]+=temp_NorS
            # output lines and output major alt
            temp_snp_line_pass = 'PASS'
            if CHR == CHR_old:
                # same CHR
                POS_DIS = abs(POS - POS_old)
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
            else:
                # diff CHR first SNP
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
            POS_old = POS
            CHR_old = CHR
            temp_snp_line.append(CHR)
            temp_snp_line.append(str(POS))
            temp_snp_line.append(REF)
            temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
            vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
            vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
            vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
            i = 9
            j = 0
            SNP_alignment[reference_name] += REF
            for genomename in SNP_alignment:
                if genomename != reference_name:
                    if i in cluster_sub:
                        SNP_alignment[genomename] += SNP_seq[j]
                        j += 1
                    i += 1
    return [CHR_old,POS_old]

def SNP_check_all_fq(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    if Major_ALT[0] != REF:
                        SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                        SNP.add(genome_order)
            sample_num += 1
        # calculate NS
        gene_info = contig_to_gene(CHR, POS)
        if gene_info != []:
            Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
            if Ref_seq_chr != 'None':
                #  observed NS ratio calculated
                temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                if codon_start <= POS_gene - 1:
                    Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                    Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                    SNP_seq_chr = Ref_seq_chr
                    if len(Ref_seq_codon) == 3:
                        Ref_seq_aa = translate(Ref_seq_codon)[0]
                        temp_snp_line_AA += Ref_seq_aa
                        ALT_set = allels_set
                        ALT_set.remove(REF)
                        for ALT in ALT_set:
                            SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                            SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_aa = translate(SNP_seq_codon)[0]
                            temp_snp_line_AA += SNP_seq_aa
                            temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                            temp_snp_line_NS[-1] += temp_NorS
        # output lines and output major alt
        temp_snp_line_pass = 'PASS'
        if CHR == CHR_old:
            # same CHR
            POS_DIS = abs(POS - POS_old)
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
        else:
            # diff CHR first SNP
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
        POS_old = POS
        CHR_old = CHR
        temp_snp_line.append(CHR)
        temp_snp_line.append(str(POS))
        temp_snp_line.append(REF)
        temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
        vcf_file_list.append(
            '\t'.join(temp_snp_line) + '\t' + '\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass, '\t'.join(temp_snp_line_NS),
                temp_snp_line_AA))
        vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
        vcf_file_list_vcf.append('\t'.join(lines_set[0:9]) + '\t' + '\t'.join(lines_set_sub) + '\n')
        i = 9
        j = 0
        SNP_alignment[reference_name] += REF
        for genomename in SNP_alignment:
            if genomename != reference_name:
                if i in cluster_sub:
                    SNP_alignment[genomename] += SNP_seq[j]
                    j += 1
                i += 1
    return [CHR_old,POS_old]

def load_ref_vcf(ref_vcf_file):
    ref_chr = dict()
    for files in ref_vcf_file:
        Set_length = False
        for lines in open(files,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, Notused, REF, ALT = lines_set[0:5]
            CHR_POS = '%s__%s'%(CHR, POS)
            ref_chr.setdefault(CHR_POS,[])
            ref_chr[CHR_POS]=[REF,ALT]
    return ref_chr

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return True
    else:
        return False
################################################### Set up ########################################################
# set up steps
SNP_cluster = dict()
cluster_set = set()
reference_set = ['reference']
outputname_set = ['filtered']

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

end_cutoff = 70 # contig end no SNP calling
################################################### Main ########################################################
# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir_merge,'*.flt.snp.vcf'))
reference_name = reference_set[0]
output_name = outputname_set[0]
for vcf_file in all_vcf_file:
    filesize = 0
    try:
        filesize = int(os.path.getsize(vcf_file + '.%s.vcf'%(output_name)))
    except FileNotFoundError:
        pass
    if filesize == 0:
        Total = 0
        donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
        SNP_tree_cmd = []
        SNP_tree_cmd2 = []
        vcf_file_list = []
        vcf_file_list_vcf = []
        Sample_name = []
        deleting_set = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        SNP_alignment = dict()
        SNP_alignment.setdefault(reference_name, '')
        cov_file_list = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        CHR_old = ''
        POS_old = 0
        SNP_cluster_donor_species = dict()
        try:
            # genome mapping file
            # need change later
            vcf_genome = glob.glob(
                output_dir_merge + '/../merge_genome/vcf/' + '_'.join(donor_species.split('_')[:-1]) + '*.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf')[0]
            ref_chr = load_ref_vcf([vcf_genome])
            for cluster_type in cluster_set:
                SNP_cluster_donor_species.setdefault(cluster_type, [])
            for lines in open(os.path.join(input_script_sub_merge, '%s.vcf.sh' % (donor_species)), 'r'):
                if lines.startswith('bcftools mpileup '):
                    # setup samples
                    sample_set = lines.split('.fasta ')[1].split('\n')[0].split(' | ')[0].split(' ')
                    samplenum = 9
                    for samples in sample_set:
                        genomename = os.path.split(samples)[-1].split(fastq_name)[0].split('all')[0].split('.sorted.bam')[0]
                        Sample_name.append(genomename.replace('.', ''))
                        if genomename in deleting_file:
                            deleting_set.append(samplenum)
                        else:
                            SNP_alignment.setdefault(genomename, '')
                            if SNP_cluster != dict() and genomename in SNP_cluster:
                                SNP_cluster_donor_species[SNP_cluster[genomename]].append(samplenum)
                        samplenum += 1
            print('running %s' % (donor_species))
            # load database
            database_file = glob.glob(os.path.join(genome_root,
                                                   '%s/*all.spades*.fna' % (donor_species)))
            if database_file == []:
                database_file = glob.glob(os.path.join(genome_root,
                                                       '%s/*all.spades*.fna' % ('_'.join(donor_species.split('_')[:-1]))))
            if database_file != []:
                Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file[0])
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    # a SNP confirmed in WGS mapping
                    CHR_POS = '%s__%s' % (CHR, POS)
                    if CHR_POS in ref_chr and not contig_end(CHR, POS):
                        Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                        if Total == 0:
                            Total = len(lines_set) - 9 - len(deleting_set)
                        if "INDEL" not in lines_set[7] \
                                and (lines_set[6] != 'LowQual'):
                            CHR_old, POS_old = SNP_check_all_fq(lines_set, '',
                                                                CHR_old, POS_old, reference_name)
            outputvcf(output_name)
            outputtree(output_name)
            try:
                f1 = open(vcf_file + '.%s.cov.txt' % (output_name), 'r')
            except IOError:
                outputcov(output_name, list(vcf_file_POS_candidate), [])
        except IndexError:
            print('missing genome mapping file %s' % (
                    output_dir_merge + '/../merge_genome/vcf/' + '_'.join(
                donor_species.split('_')[:-1]) + '*.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf'))

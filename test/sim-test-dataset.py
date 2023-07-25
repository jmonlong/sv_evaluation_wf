import random as rd
import numpy as np
from Bio.Seq import MutableSeq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

nuc = np.array(["A", "T", "C", "G"])

# simulate some genomic sequence
ref_len = 100000
seqArray = nuc[[int(rd.random()*4) for i in range(ref_len)]]
ref_seq = MutableSeq("".join(seqArray))

# write as fasta
recs = []
recs.append(SeqRecord(ref_seq, id='S', description=""))
SeqIO.write(recs, "ref.fa", "fasta")

# simulate some deletions, insertions
svs = []
for ii in range(100):
    sv = {}
    sv['pos'] = rd.randint(0, ref_len - 1)
    sv['len'] = rd.randint(50, 500)
    if rd.random() > .5:
        sv['type'] = 'DEL'
    else:
        sv['type'] = 'INS'
    svs.append(sv)


# function to write a VCF
def writeSVs(vcf_file, svs, ref_seq):
    vcff = open(vcf_file, 'w')
    vcff.write('##fileformat=VCFv4.1\n')
    vcff.write('##FORMAT=<ID=GT,Number=1,'
               'Type=String,Description="Genotype">\n'
               '##contig=<ID=S,length={}>\n'.format(ref_len))
    vcff.write('##SAMPLE=<ID=sample>\n')
    vcff.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\t'
               'FILTER\tINFO\tFORMAT\tsample\n')
    for var in svs:
        if var['type'] == 'DEL':
            REF = ref_seq[var['pos']:(var['pos']+var['len'])]
            ALT = ref_seq[var['pos']]
        else:
            REF = ref_seq[var['pos']]
            seqArray = nuc[[int(rd.random()*4) for i in range(var['len'])]]
            ins_seq = MutableSeq("".join(seqArray))
            ALT = ref_seq[var['pos']] + ins_seq
        GT = ['0/1', '1/1'][rd.randint(0, 1)]
        vcfrec = 'S\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t{}\n'
        vcff.write(vcfrec.format(var['pos']+1, REF, ALT, GT))
    vcff.close()
    return vcf_file

# pick some as "truth"
truth_svs = rd.sample(svs, int(len(svs)*.9))
writeSVs('truth.vcf', truth_svs, ref_seq)

# pick some as "calls"
calls_svs = rd.sample(svs, int(len(svs)*.9))
writeSVs('calls.vcf', calls_svs, ref_seq)

# pick regions around some variants as confident regions
reg_svs = rd.sample(svs, int(len(svs)*.7))
fl = 200
with open('conf.bed', 'wt') as outf:
    for reg in reg_svs:
        outf.write('S\t{}\t{}\n'.format(reg['pos']-fl,
                                        reg['pos']+reg['len']+fl))

# pick some around some variants as simple repeat
reg_svs = rd.sample(svs, int(len(svs)*.2))
fl = 20
with open('simprep.bed', 'wt') as outf:
    for reg in reg_svs:
        outf.write('S\t{}\t{}\n'.format(reg['pos']-fl,
                                        reg['pos']+reg['len']+fl))

from pysam import VariantFile
import sys
#annotate_vcf_repeatmasker.vcf <vcf_file> <new_vcf_file>
import subprocess

vfc_filename = sys.argv[1]
vfc_out_filename = sys.argv[2]

vcf_in_file = VariantFile(vfc_filename)
header = vcf_in_file.header
if 'REPEATTYPE' not in header.info:
    header.info.add('REPEATTYPE', number=1, type='String', description="Repeat classification of variant content")

temp_fasta = "temp.fasta"
with open(temp_fasta, 'w') as temp_fasta:
    for rec in vcf_in_file.fetch():
        temp_fasta.write(">" + rec.id + "\n")
        temp_fasta.write(rec.info['SEQ'] + "\n")
subprocess.call(["/programs/RepeatMasker_4-1-0/RepeatMasker", "temp.fasta"])

d = {}
with open("contigs.fasta.out", 'r') as repeat_master_out:
    repeat_master_out.readline()
    repeat_master_out.readline()
    repeat_master_out.readline()
    for l in repeat_master_out.readlines():
        split = l.split()
        id = split[4]
        type = split[10]
        if id not in d:
            d[id] = []
        d[id].append(type)
d_final = {}
for k,v in d.items():
    repeat_types = list(set(v))
    if len(repeat_types) > 1:
        d_final[k] = "Complex"
    else:
        d_final[k] = repeat_types[0]

vcf_out_file = VariantFile(vfc_out_filename, 'w', header=header)

for rec in vcf_in_file.fetch():
    if 'SVLEN' not in rec.info:
        if rec.id in d_final:
            rec.info['REPEATTYPE'] = d_final[rec.id]
        else:
            rec.info['REPEATTYPE'] = "NotMasked"

    vcf_out_file.write(rec)


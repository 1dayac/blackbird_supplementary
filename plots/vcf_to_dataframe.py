from pysam import VariantFile
import sys
#annotate_vcf_repeatmasker.vcf <vcf_file_tp> <vcf_file_fn> <Method Name> <Insetions/Deletions> <out dataframe>
import subprocess

vfc_filename_tp_base = sys.argv[1]
vfc_filename_fn = sys.argv[2]
method_name = sys.argv[3]
sv_type = sys.argv[4]
out_dataframe = sys.argv[5]

vcf_in_file1 = VariantFile(vfc_filename_tp_base)
correct_d = {}
for rec in vcf_in_file1.fetch():
    if rec.info['REPEATTYPE'] not in correct_d:
        correct_d[rec.info['REPEATTYPE']] = 0
    correct_d[rec.info['REPEATTYPE']] += 1

vcf_in_file1 = VariantFile(vfc_filename_fn)
wrong_d = {}

for rec in vcf_in_file1.fetch():
    if rec.info['REPEATTYPE'] not in wrong_d:
        wrong_d[rec.info['REPEATTYPE']] = 0
    wrong_d[rec.info['REPEATTYPE']] += 1

with open(out_dataframe, 'w') as out_handle:
    out_handle.write("Method,SV type,RepeatType,Called,NotCalled\n")
    for k,v in correct_d.items():
        correct = v
        if k in wrong_d:
            wrong = wrong_d[k]
        else:
            wrong = 0
        out_handle.write(",".join([method_name, sv_type, k, str(correct), str(wrong)]))
        out_handle.write('\n')

    for k,v in wrong.items():
        if k not in correct_d:
            out_handle.write(",".join([method_name, sv_type, k, str(0), str(v)]))
            out_handle.write('\n')

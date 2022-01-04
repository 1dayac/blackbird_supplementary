import sys
info_message = ""


with open(sys.argv[1], 'r') as input_vcf:
    with open(sys.argv[2], 'w') as output_vcf:
        for r in input_vcf.readlines():
            if not r.startswith("#"):
                r = r.strip()
                r += "\tGT\t./.\n"
                output_vcf.write(r)
                continue
            if r.startswith("##"):
                output_vcf.write(r)
                if r.startswith("##INFO=<ID=MATEID"):
                    output_vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"")
                continue
            if r.startswith("#"):
                r = r.strip()
                r += "\tFORMAT\tSAMPLE1\n"
                output_vcf.write(r)
                continue
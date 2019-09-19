import sys
from Bio import SeqIO, Seq

records_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))
chrom = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

sub_ref = records_dict[chrom][start:end]
records_dict[chrom + "_" + str(start) + "_" + str(end)] = sub_ref


with open(sys.argv[5], "w") as output_handle:
    for record in records_dict.values():
        if record.id == chrom + "_" + str(start) + "_" + str(end):
            SeqIO.write(record, output_handle, "fasta")

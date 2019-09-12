svs = []
region_defined = True
region_chr = "chr1"
region_start = 0
region_end = 300000000

class SV:
    def __init__(self, chrom, pos, length, svtype, repeat_type, seq = ""):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.checked = False
        self.svtype = svtype
        self.repeat_type = repeat_type
        self.seq = seq.upper()

def get_seq(line):
    start_pos = line.find("SEQ")
    ans = ""
    start_pos += 4
    while line[start_pos] != ';':
        ans += line[start_pos]
        start_pos += 1
    return ans


def get_len(line):
    start_pos = line.find("SVLEN")
    ans = ""
    start_pos += 6
    while line[start_pos].isdigit():
        ans += line[start_pos]
        start_pos += 1
    if ans == "":
        return 0
    return int(ans)

def get_svtype(line):
    start_pos = line.find("SVTYPE")
    ans = ""
    start_pos += 7
    while line[start_pos] != ';' and line[start_pos] != '\n':
        ans += line[start_pos]
        start_pos += 1
    if ans == "":
        return "UNKNOWN"
    return ans

def get_repeattype(line):
    start_pos = line.find("REPEAT_TYPE")
    ans = ""
    start_pos += 12
    while line[start_pos] != ';':
        ans += line[start_pos]
        start_pos += 1
    if ans == "":
        return "UNKNOWN"
    return ans


sv_dict = []

with open("results/chm1/pacbio.vcf", "r") as pacbio:
    for r in pacbio.readlines():
        if r.startswith("#"):
            continue
        splitted = r.split("\t")
        if region_defined:
            if region_chr != splitted[0] or (region_start >= int(splitted[1]) >= region_end) :
                continue
        sv = SV(splitted[0], int(splitted[1]), get_len(r), get_svtype(r), get_repeattype(r), get_seq(r))
        sv_dict.append(sv)





from Bio import SeqIO, Seq

records_dict = SeqIO.to_dict(SeqIO.parse("chr1.fa", "fasta"))

#with open("chr1.fa", "w") as output_handle:
#    for record in records_dict.values():
#        if record.id == "chr1":
#            SeqIO.write(record, output_handle, "fasta")

#records_dict = SeqIO.to_dict(SeqIO.parse("insertions_with_anchors.fasta", "fasta"))
#print(len(records_dict["chr1"]))
#records_dict["chr1"] = records_dict["chr1"][:100] + Seq.Seq("aaaaaaa") + records_dict["chr1"][100:]
#print(len(records_dict["chr1"]))
#exit()

shift = 0
for sv in sv_dict:
    print("Old " + str(len(records_dict["chr1"])))
    if sv.svtype == "INS":
        new_seq = records_dict["chr1"][: sv.pos + shift] + sv.seq + records_dict["chr1"].seq[sv.pos + shift:]
        shift += sv.length

    if sv.svtype == "DEL":
        new_seq = records_dict["chr1"][: sv.pos + shift] + records_dict["chr1"].seq[sv.pos + shift + sv.length:]
        shift -= sv.length

    records_dict["chr1"] = new_seq
    print("New " + str(len(records_dict["chr1"])))



with open("chr1.extended.fa", "w") as output_handle:
    for record in records_dict.values():
        if record.id == "chr1":
            SeqIO.write(record, output_handle, "fasta")



class SV:
    def __init__(self, chrom, pos, length, svtype, repeat_type, seq = ""):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.checked = False
        self.svtype = svtype
        self.repeat_type = repeat_type
        self.seq = seq.upper()


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



region_defined = True
region_chr = "chr1"
region_start = 0
region_end = 300000000

sv_dict = {}
repeat_type_dict = {}
repeat_type_found_dict = {}

total_length = 0
with open("results/simulated/pacbio.vcf", "r") as pacbio:
    for r in pacbio.readlines():
        if r.startswith("#"):
            continue
        splitted = r.split("\t")
        if region_defined:
            if region_chr != splitted[0] or (region_start >= int(splitted[1]) >= region_end) :
                continue
        sv = SV(splitted[0], int(splitted[1]), get_len(r), get_svtype(r), get_repeattype(r))
        total_length += sv.length
        if get_svtype(r) + " " + get_repeattype(r) not in repeat_type_dict:
            repeat_type_dict[get_svtype(r) + " " + get_repeattype(r)] = 0
        repeat_type_dict[get_svtype(r) + " " + get_repeattype(r)] += 1
        if splitted[0] not in sv_dict:
            sv_dict[splitted[0]] = []
        sv_dict[splitted[0]].append(sv)
print("Total PacBio length - " + str(total_length))
def Near(sv1, sv2):
    return abs(sv1.pos - sv2.pos) <= 100


for type in repeat_type_dict.keys():
    repeat_type_found_dict[type] = 0

manta = 0
with open("results/simulated/manta.vcf", "r") as manta_vcf:
    for r in manta_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, 0, "", "", "")
        found = False

        for sv in sv_dict[chrom]:
            if sv.checked:
                continue
            if Near(sv, new_sv):
                found = True
                sv.checked = True
                manta += 1
                break
print("Manta" + str(manta))

not_in_pacbio = 0
with open("results/simulated/blackbird.vcf", "r") as my_vcf:
    for r in my_vcf.readlines():
        if r.startswith("#"):
            continue
        chrom = r.split("\t")[0]
        pos = int(r.split("\t")[1])
        new_sv = SV(chrom, pos, get_len(r), get_svtype(r), "UNKNOWN")

        found = False
        if chrom not in sv_dict:
            sv_dict[chrom] = []

        for sv in sv_dict[chrom]:
            if sv.checked:
                continue
            if Near(sv, new_sv):
                found = True
                sv.checked = True
                repeat_type_found_dict[sv.svtype + " " + sv.repeat_type] += 1
                break
        if not found:
            not_in_pacbio += 1



print("Found only in Blackbird - " + str(not_in_pacbio))

total = 0
total_found = 0
for repeat in repeat_type_dict.keys():
    print(repeat + ": " + str(repeat_type_found_dict[repeat]) + "/" + str(repeat_type_dict[repeat]))
    total += repeat_type_dict[repeat]
    total_found += repeat_type_found_dict[repeat]


print("Total: " + str(total_found) + "/" + str(total))
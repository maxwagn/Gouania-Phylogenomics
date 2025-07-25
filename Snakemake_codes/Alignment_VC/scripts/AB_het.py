import sys
from Bio import bgzf
from pysam import VariantFile

bcf_input = sys.argv[1] #"20210222test.variants.raw.prefilter.chr1_test.bcf"
bcf_output = sys.argv[2]

bcf_in = VariantFile(bcf_input)  # auto-detect input format

bcf_in.header.add_meta('INFO', items=[('ID',"AB_het"), ('Number',1), ('Type','Float'), ('Description','Allelic Balance Heterozygotes')])

bcf_outfile = bgzf.BgzfWriter(bcf_output)
bcf_outfile.write(str(bcf_in.header).strip())

for line in bcf_in.fetch():
    line_list = str(line).rstrip()
    line_list = line_list.split('\t')
    DP_inx = line_list[8].split(":").index("DP")
    AD_inx = line_list[8].split(":").index("AD")
    out = str('\t'.join(line_list[:7])) + "\t"
    AB_Het = 0
    DP_Het = 0
    AD_Het = 0
    for indiv in line_list[9:]:
        if indiv.startswith("0/1") or indiv.startswith('1/0') \
        or indiv.startswith('1|0') or indiv.startswith('0|1'): 
            DP_Het = DP_Het + int(indiv.split(":")[DP_inx])
            AD_Het_list = indiv.split(":")[AD_inx].split(",")[1:] #all alt alleles
            AD_Het = AD_Het + sum([int(i) for i in AD_Het_list])
        else:
            pass
    if AD_Het > 0:
        AB_Het = AD_Het/DP_Het
        AB_str = ";AB_het=" + str(AB_Het)    
        info_line = line_list[7] + ";AB_het=" + str(AB_Het)
        out += info_line + "\t" + str('\t'.join(line_list[8:]))
    else:
        out += str(line_list[7]) + "\t" + str('\t'.join(line_list[8:]))

    bcf_outfile.write("\n" + out)

bcf_outfile.close()



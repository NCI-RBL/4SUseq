import subprocess
import json
import sys

# python {params.pyscript} {params.sample} {output.json} \
#     {input.rawlanestxt} \
#     {input.trimmedlanestxt} \
#     {input.fastuniqlanestxt} \
#     {input.postsecondarysupplementaryfilterbamflagstat} \
#     {input.postinsertionfilterbamflagstat} \
#     {input.postmapqwidowfilterbamflagstat} \
#     {input.pluspostmapqwidowfilterbamflagstat} \
#     {input.minuspostmapqwidowfilterbamflagstat} \
#     {input.plustoSNPcallingbamflagstat} \
#     {input.minustoSNPcallingbamflagstat} \
#     {input.mutatedbamflagstat} \
#     {input.unmutatedbamflagstat} \
#     {input.vcf} \
#     {input.plusvcf} \
#     {input.minusvcf}

def get_nreads_from_labels(filename,samplename):
	f=open(filename)
	for l in f.readlines():
		l=l.strip().split("\t")
		if l[0]=="#" and l[1]==samplename:
			return l[2]

def get_nfragments_from_flagstat(filename):
	nf=int(float(subprocess.check_output("grep properly "+filename+"|awk '{print $1/2}'",shell=True).strip()))
	return nf


def get_nmutations_from_vcf(filename):
	nm=int(float(subprocess.check_output("zcat "+filename+" | grep -v ^# |wc -l",shell=True).strip()))
	return nm



samplename=sys.argv[1]  #sample name
outfile=sys.argv[2]     #output json name.... sample.json
rawlanestxt=sys.argv[3] # raw lanes.txt
trimmedlanestxt=sys.argv[4] # trimmed lanes.txt
fastuniqlanestxt=sys.argv[5] # fastuniq lanes.txt
postsecondarysupplementaryfilterbamflagstat=sys.argv[6]
postinsertionfilterbamflagstat=sys.argv[7]
pluspostmapqwidowfilterbamflagstat=sys.argv[8]
minuspostmapqwidowfilterbamflagstat=sys.argv[9]
plustoSNPcallingbamflagstat=sys.argv[10]
minustoSNPcallingbamflagstat=sys.argv[11]
mutatedbamflagstat=sys.argv[12]
unmutatedbamflagstat=sys.argv[13]
vcf=sys.argv[14]
plusvcf=sys.argv[15]
minusvcf=sys.argv[16]

data=dict()
data['samplename']=samplename
data['nfragments']=dict()
data['nfragments']['raw'] = get_nreads_from_labels(rawlanestxt,samplename)
data['nfragments']['trim'] = get_nreads_from_labels(trimmedlanestxt,samplename)
data['nfragments']['fastuniq'] = get_nreads_from_labels(fastuniqlanestxt,samplename)
data['nfragments']['post_secondary_supplementary_filter'] = get_nfragments_from_flagstat(postsecondarysupplementaryfilterbamflagstat)
data['nfragments']['post_insertion_filter'] = get_nfragments_from_flagstat(postinsertionfilterbamflagstat)
data['nfragments']['post_mapq_widow_filter_plus_strand'] = get_nfragments_from_flagstat(pluspostmapqwidowfilterbamflagstat)
data['nfragments']['post_mapq_widow_filter_minus_strand'] = get_nfragments_from_flagstat(minuspostmapqwidowfilterbamflagstat)
data['nfragments']['to_SNP_calling_plus_strand'] = get_nfragments_from_flagstat(plustoSNPcallingbamflagstat)
data['nfragments']['to_SNP_calling_minus_strand'] = get_nfragments_from_flagstat(minustoSNPcallingbamflagstat)
data['nfragments']['mutated'] = get_nfragments_from_flagstat(mutatedbamflagstat)
data['nfragments']['unmutated'] = get_nfragments_from_flagstat(unmutatedbamflagstat)
data['nmutations'] = get_nmutations_from_vcf(vcf)
data['nmutations_plus'] = get_nmutations_from_vcf(plusvcf)
data['nmutations_minus'] = get_nmutations_from_vcf(minusvcf)
out_file = open(outfile, "w") 
json.dump(data, out_file, indent = 6) 
out_file.close()

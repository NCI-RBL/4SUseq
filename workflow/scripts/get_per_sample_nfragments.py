import subprocess
import json
import sys

def get_nreads_from_labels(filename,samplename):
	f=open(filename)
	for l in f.readlines():
		l=l.strip().split("\t")
		if l[0]=="#" and l[1]==samplename:
			return l[2]

samplename=sys.argv[1]  #sample name
outfile=sys.argv[2]     #output json name.... sample.json
data=dict()
data['samplename']=samplename
data['nfragments']=dict()
data['nfragments']['raw'] = get_nreads_from_labels("raw_fastq/lanes.txt",samplename)
data['nfragments']['trim'] = get_nreads_from_labels("trim/lanes.txt",samplename)
data['nfragments']['fastuniq'] = get_nreads_from_labels("fastuniq/lanes.txt",samplename)
data['nfragments']['post_secondary_supplementary_filter'] = int(float(subprocess.check_output("grep properly hisat2/"+samplename+".post_secondary_supplementary_filter.bam.flagstat|awk '{print $1/2}'",shell=True).strip()))
data['nfragments']['post_insertion_filter'] = int(float(subprocess.check_output("grep properly hisat2/"+samplename+".post_insertion_filter.bam.flagstat|awk '{print $1/2}'",shell=True).strip()))
data['nfragments']['post_mapq_widow_filter'] = int(float(subprocess.check_output("grep properly hisat2/"+samplename+".toSNPcalling.bam.flagstat|awk '{print $1/2}'",shell=True).strip()))
data['nmutations'] = int(float(subprocess.check_output("zcat vcf/"+samplename+".TtoC.vcf.gz | grep -v ^# |wc -l",shell=True).strip()))
out_file = open(outfile, "w") 
json.dump(data, out_file, indent = 6) 
out_file.close()

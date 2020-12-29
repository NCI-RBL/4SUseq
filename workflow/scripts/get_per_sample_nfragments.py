import subprocess
import json
import sys

samplename=sys.argv[1]  #sample name
outfile=sys.argv[2]     #output json name.... sample.json
data=dict()
data['samplename']=samplename
data['nfragments']=dict()
data['nfragments']['raw'] = int(subprocess.check_output("grep \b"+samplename+"\b raw_fastq/lanes.txt|grep ^#|awk '{print $3}'",shell=True).strip())
data['nfragments']['trim'] = int(subprocess.check_output("grep \b"+samplename+"\b trim/lanes.txt|grep ^#|awk '{print $3}'",shell=True).strip())
data['nfragments']['fastuniq'] = int(subprocess.check_output("grep \b"+samplename+"\b fastuniq/lanes.txt|grep ^#|awk '{print $3}'",shell=True).strip())
data['nfragments']['post_secondary_supplementary_filter'] = int(float(subprocess.check_output("grep properly hisat2/"+samplename+".post_secondary_supplementary_filter.bam.flagstat|awk '{print $1/2}'",shell=True).strip()))
data['nfragments']['post_insertion_filter'] = int(float(subprocess.check_output("grep properly hisat2/"+samplename+".post_insertion_filter.bam.flagstat|awk '{print $1/2}'",shell=True).strip()))
data['nfragments']['toSNPcalling'] = int(float(subprocess.check_output("grep properly hisat2/"+samplename+".toSNPcalling.bam.flagstat|awk '{print $1/2}'",shell=True).strip()))
data['nmutations'] = int(float(subprocess.check_output("grep -v ^# vcf/test.TtoC.vcf |wc -l",shell=True).strip()))
out_file = open(outfile, "w") 
json.dump(data, out_file, indent = 6) 
out_file.close()

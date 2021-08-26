import HTSeq
import sys
import argparse
import os
def get_sname(s):
	sname=s.name
	sname=sname.split()[0]
	return sname

parser = argparse.ArgumentParser(description='Filter FASTQ by readids')
parser.add_argument('--infq', dest='infq', type=str, required=True,
                    help='input FASTQ file')
parser.add_argument('--infq2', dest='infq2', type=str, required=True,
                    help='input FASTQ file')
parser.add_argument('--outfq', dest='outfq', type=str, required=True,
                    help='filtered output FASTQ file')
parser.add_argument('--outfq2', dest='outfq2', type=str, required=True,
                    help='filtered output FASTQ file')
parser.add_argument('--readids', dest='readids', type=str, required=True,
                    help='file with readids to keep (one readid per line)')
parser.add_argument('--complement', dest='complement', action='store_true', 
                    help='complement the readid list, ie., include readids NOT in the list')
args = parser.parse_args()
rids=set(map(lambda x:x.strip(),open(args.readids,'r').readlines()))
sequences = dict( (get_sname(s), s) for s in HTSeq.FastqReader(args.infq))
sequences2 = dict( (get_sname(s), s) for s in HTSeq.FastqReader(args.infq2))
if args.complement:
	rids=set(sequences.keys())-rids
outfqfilename=args.outfq
outfqfilename2=args.outfq2
dummy=outfqfilename.strip().split(".")
dummy2=outfqfilename2.strip().split(".")
if dummy[-1]=="gz":
	dummy.pop(-1)
	outfqfilename=".".join(dummy)	
if dummy2[-1]=="gz":
	dummy2.pop(-1)
	outfqfilename2=".".join(dummy2)	
outfqfile = open(outfqfilename,'w')
outfqfile2 = open(outfqfilename,'w')
for rid in rids:
	s=sequences[rid]
	s2=sequences2[rid]
	s.write_to_fastq_file(outfqfile)
	s2.write_to_fastq_file(outfqfile2)
outfqfile.close()
outfqfile2.close()
if dummy[-1]=="gz":	
	os.system("pigz -p4 -f "+outfqfilename)
if dummy2[-1]=="gz":	
	os.system("pigz -p4 -f "+outfqfilename2)

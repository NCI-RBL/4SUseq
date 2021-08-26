import HTSeq
import sys
import argparse
import os
parser = argparse.ArgumentParser(description='Filter FASTQ by readids')
parser.add_argument('--infq', dest='infq', type=str, required=True,
                    help='input FASTQ file')
parser.add_argument('--outfq', dest='outfq', type=str, required=True,
                    help='filtered output FASTQ file')
parser.add_argument('--readids', dest='readids', type=str, required=True,
                    help='file with readids to keep (one readid per line)')
parser.add_argument('--complement', dest='complement', action='store_true', 
                    help='complement the readid list, ie., include readids NOT in the list')
args = parser.parse_args()
rids=list(map(lambda x:x.strip(),open(args.readids,'r').readlines()))
outfqfile = open(args.outfq,'w')
for s in HTSeq.FastqReader(args.infq):
	sname=s.name
	sname=sname.split()[0]
	if args.complement:
		if not sname in rids:
			s.write_to_fastq_file(outfqfile)
	else:
		if sname in rids:
			s.write_to_fastq_file(outfqfile)
outfqfile.close()	

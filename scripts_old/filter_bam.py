import pysam,sys
import argparse

def get_ninsertions(ctuples):
	for t in ctuples:
		operation,length=t
		if operation==1: # operation == 1 is insertion (see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats)
			return length
	return 0

parser = argparse.ArgumentParser(description='filter PE bamfile by mapQ values')
parser.add_argument('-i',dest='inBam',required=True,help='Input Bam File')
parser.add_argument('-o',dest='outBam',required=True,help='Output Bam File')
parser.add_argument('-q',dest='mapQ',type=int,required=False,help='mapQ value ... default 6, maqQ of 5 or lower will be discarded',default=6)
parser.add_argument('-n',dest='ninsertion',type=int,required=False,help='number of insertions allowed ... default 6... reads with 5 or fewer insertions will be kept',default=6)
args = parser.parse_args()
samfile = pysam.AlignmentFile(args.inBam, "rb")
mapq=dict()
for read in samfile.fetch():
	if read.is_unmapped:
		continue
	if read.is_supplementary:
		continue
	if read.is_secondary:
		continue
	if read.is_duplicate:
		continue
	if read.is_proper_pair:
		if get_ninsertions(read.cigartuples) >= args.ninsertion:
			continue
		if read.mapping_quality < args.mapQ and read.query_name in mapq:
			del mapq[read.query_name]
		if read.mapping_quality >= args.mapQ:
			if not read.query_name in mapq:
				mapq[read.query_name]=0
			if read.is_read1:
				mapq[read.query_name]+=1
			if read.is_read2:
				mapq[read.query_name]+=1
samfile.close()
samfile = pysam.AlignmentFile(args.inBam, "rb")
pairedreads = pysam.AlignmentFile(args.outBam, "wb", template=samfile)
for read in samfile.fetch():
	if read.query_name in mapq:
		if mapq[read.query_name]!=2:
			continue
		pairedreads.write(read)
samfile.close()
pairedreads.close()

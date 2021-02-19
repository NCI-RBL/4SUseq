import pysam,sys
import argparse

parser = argparse.ArgumentParser(description='filter PE bamfile by readids')
parser.add_argument('-i',dest='inBam',required=True,help='Input Bam File')
parser.add_argument('-o',dest='outBam',required=True,help='Output Bam File')
parser.add_argument('-o2',dest='outBam2',required=False,help='Complement Output Bam File')
parser.add_argument('--readids', dest='readids', type=str, required=True,
                    help='file with readids to keep (one readid per line)')
args = parser.parse_args()

rids=set(map(lambda x:x.strip(),open(args.readids,'r').readlines()))
samfile = pysam.AlignmentFile(args.inBam, "rb")
alignments = dict()
for read in samfile.fetch():
	q=read.query_name
	q=q.strip()
	# print("##"+q+"##")
	try:
		alignments[q].append(read)
	except:
		alignments[q]=list()
		alignments[q].append(read)
# print(alignments)
count1=0
count2=0
outsamfile = pysam.AlignmentFile(args.outBam, "wb", template=samfile)
for rid in rids:
	# print(rid in alignments)
	try:
		for r in alignments[rid]:
			count1+=1
			outsamfile.write(r)
	except:
		continue
outsamfile.close()
if args.outBam2:
	outsamfile2 = pysam.AlignmentFile(args.outBam2, "wb", template=samfile)
	rids2=set(alignments.keys())
	complementary_rids=rids2-rids
	for rid in complementary_rids:
		for r in alignments[rid]:
			count2+=1
			outsamfile2.write(r)
	outsamfile2.close()
samfile.close()
# print(count1,count2)

#collapse nreads based on lanes and flowcells
#input --> lanes.txt
import fileinput
samples=dict()

for line in fileinput.input():
	l=line.strip().split()
	samplename=l[0]
	nreads=l[1]
	flowcell=l[2]
	lane=l[3]
	if not samplename in samples:
		samples[samplename]=dict()
		samples[samplename]["nreads"]=0
		samples[samplename]["flowcell_lane"]=list()
	samples[samplename]["nreads"]+=int(nreads)
	samples[samplename]["flowcell_lane"].append(flowcell+"_"+lane)

for k in sorted(samples):
	v=samples[k]
	print("%s\t%d\t%s"%( k,v["nreads"],",".join(v["flowcell_lane"])))

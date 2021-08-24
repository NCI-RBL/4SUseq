import json
import pandas
import glob
import sys

x=[]
for jsonfile in glob.glob("*.json"):
	with open(jsonfile) as f:
		x.append(json.load(f))

df=pandas.json_normalize(x,max_level=1)
df.to_csv(sys.argv[1],sep="\t",index=False)  #sys.argv[1] is the output filename

#!/usr/bin/env python

import time
import cPickle as pickle

start=time.clock()
data=dict()

FILE=open("lambda_virus.fa", 'r')
FILE.readline()
line=""
for subline in FILE.readlines():
	line+=subline.strip()
FILE.close()

for i in range(len(line)-7):
	word=line[i:i+8]
	if word in data:
		data[word]=data.get(word)+"\t"+str(i)
	else:
		data[word]=str(i)

OUT=open("dict.pickle", 'w')
pickle.dump(data, OUT)
OUT.close()

print "Time elapsed: %fs" % (time.clock()-start)

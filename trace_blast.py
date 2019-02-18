#!/usr/bin/env python

import time, string
import cPickle as pickle

def search(item):
	result1, result2=alignment(query_seq, genome_seq[item:item+len(query_seq)])
	result=""
	score = 0
	for k in range(0, len(result1)):
		if result1[k] == result2[k]:
			score += 1
			result+="|"
		else:
			result+=" "
	score = float(score)/len(result1)
	
	if score > 0.8:
		print "Query sequence found at "+str(item)+"-"+str(item+len(query_seq))+":\n"
#		print query_seq+"\n"+genome_seq[item:item+len(query_seq)]+"\n"
		print "Query:  "+result1+"\n        "+result+"\nGenome: "+result2+"\nIdentity score: "+str(score)+"\n"

def alignment(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    g = -3
    matrix = []
    trace = []
    for i in range(0, m+1):
        tmp = []
        for j in range(0, n+1):
            tmp.append(0)
        matrix.append(tmp)
        trace.append(tmp)
    for sii in range(0, m+1):
        matrix[sii][0] = sii*g
        trace[sii][0] = 0
    for sjj in range(0, n+1):
        matrix[0][sjj] = sjj*g
        trace[0][sjj] = 0
    for siii in range(1, m+1):
        for sjjj in range(1, n+1):
            if seq1[siii-1] == seq2[sjjj-1]:
                matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] +2, matrix[siii][sjjj-1] + g)
		if(matrix[siii][sjjj] == matrix[siii - 1][sjjj - 1] +2):
			trace[siii][sjjj]=3
		elif(matrix[siii][sjjj] == matrix[siii - 1][sjjj] +g):
			trace[siii][sjjj]=1
		else:
			trace[siii][sjjj]=2
            else:
                matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] -1, matrix[siii][sjjj-1] + g)
		if(matrix[siii][sjjj] == matrix[siii - 1][sjjj - 1] -1):
			trace[siii][sjjj]=3
		elif(matrix[siii][sjjj] == matrix[siii - 1][sjjj] +g):
			trace[siii][sjjj]=1
		else:
			trace[siii][sjjj]=2
    sequ1 = []
    sequ2 = []
    x=1
    y=1
    while x != m+1 and y != n+1:
        if trace[x][y]==3:
            sequ1.append(seq1[x-1])
            sequ2.append(seq2[y-1])
            x+=1
            y+=1
        elif trace[x][y]==1:
            sequ1.append('-')
            sequ2.append(seq2[y-1])
            y+=1
        else:
            sequ1.append(seq1[x-1])
            sequ2.append('-')
            x+=1
#    sequ1.reverse()
#    sequ2.reverse()
    result1 = string.join(sequ1, '')
    result2 = string.join(sequ2, '')

#    print str(len(result1))+"\t"+str(len(result2))+"\n"
#    if(len(result1)>len(result2)):
#	result2+='-'*(len(result1)-len(result2))
#    elif len(result1)<len(result2):
#        result1+='-'*(len(result2)-len(result1))

#    score = 0
#    for k in range(0, len(result1)):
#        if result1[k] == result2[k]:
#            score += 1
#    score = float(score)/len(result1)
    return result1, result2

if __name__ == "__main__":
	start=time.time()
	FILE=open("dict.pickle", 'r')
	data=pickle.load(FILE)
	FILE.close()
	global query_seq,genome_seq
	FILE=open("search1.txt", 'r')
	query_seq=FILE.read()
	query_seq=query_seq.strip()
	FILE.close()
	FILE=open("lambda_virus.fa", 'r')
	line=FILE.readline()
	genome_seq=""
	for line in FILE.readlines():
		genome_seq+=line.strip("\n")
	FILE.close()
	loci=set()
	for i in range(len(query_seq)-7):
		word=query_seq[i:i+8]
		if word in data:
			for item in data[word].split("\t"):
				loci.add(int(item)-i)
	#print "loci size: "+str(len(loci))+"\n"
	#print loci
	for item in loci:
		search(item)
	print "Time elapsed: %fs" % (time.time()-start)
#	score, result1, result2=alignment("CATGAGGTTGCCCCGTATTCAGTGTCGCTGATTTGTATTGTCTG", "CATGAGGTTGCCCCGTATTCAGTGTCGCTGATTTGTATTG")
#	print result1+"\n"+result2+"\n"+str(score)+"\n"

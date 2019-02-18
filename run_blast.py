#!/usr/bin/env python

import time, string
import cPickle as pickle

def search(item):
	result1, result2=alignment(query_seq, genome_seq[item:item+len(query_seq)])
	extend=0
	for i in range(len(result2)):
		if(result2[i] == '-'):
			extend+=1
		else:
			break
	start=item-extend
	result2=genome_seq[start:item]+result2[extend:]
	extend=0
	for i in range(len(result2)):
		if(result2[-i-1] == '-'):
			extend+=1
		else:
			break
	if extend != 0:
		result2=result2[:-extend]+genome_seq[item+len(query_seq):item+len(query_seq)+extend]
	result=""
	score = 0
	for k in range(len(result1)):
		if result1[k] == result2[k]:
			score += 1
			result+="|"
		else:
			result+=" "
	score = float(score)/len(result1)
	
	if score > 0.8:
		print "Query sequence found at "+str(start)+"-"+str(item+len(query_seq)+item+extend)+":\n"
		print "Query:  "+result1+"\n        "+result+"\nGenome: "+result2+"\nIdentity score: "+str(score)+"\n"

def alignment(seq1, seq2):
	m = len(seq1)
	n = len(seq2)
	g = -3
	matrix = []
	tmp = []
	for j in range(0, n+1):
		tmp.append(0)
	for i in range(0, m+1):
		matrix.append(tmp[:])
	for sii in range(0, m+1):
		matrix[sii][0] = sii*g
	for sjj in range(0, n+1):
		matrix[0][sjj] = sjj*g
	for siii in range(1, m+1):
		for sjjj in range(1, n+1):
			if seq1[siii-1] == seq2[sjjj-1]:
				matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] +2, matrix[siii][sjjj-1] + g)
			else:
				matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] -1, matrix[siii][sjjj-1] + g)
	sequ1 = []
	sequ2 = []
	while m > 0 and n > 0:
		if (matrix[m-1][n-1]>=matrix[m-1][n] and matrix[m-1][n-1]>=matrix[m][n-1]) or seq1[m-1] == seq2[n-1]:
			sequ1.append(seq1[m-1])
			sequ2.append(seq2[n-1])
			m -= 1
			n -= 1
		elif matrix[m][n-1]>=matrix[m-1][n-1] and matrix[m][n-1]>=matrix[m-1][n]:
			sequ1.append('-')
			sequ2.append(seq2[n-1])
			n -= 1
		else:
			sequ1.append(seq1[m-1])
			sequ2.append('-')
			m -= 1
	if m>n:
		for i in range(m-n):
			sequ1.append(seq1[m-n-1-i])
			sequ2.append('-')
	elif n>m:
		for i in range(n-m):
			sequ2.append(seq2[n-m-1-i])
			sequ1.append('-')

	sequ1.reverse()
	sequ2.reverse()
	result1 = string.join(sequ1, '')
	result2 = string.join(sequ2, '')
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
				ifPass=1
				for j in range(-2,3):
					if int(item)-i+j in loci:
						ifPass=0
				if ifPass:
					loci.add(int(item)-i)
	print "loci size: "+str(len(loci))+"\n"
	print loci
	for item in loci:
		search(item)
	print "Time elapsed: %fs" % (time.time()-start)

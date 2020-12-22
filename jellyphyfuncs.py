import os, sys
import time
import argparse
from collections import namedtuple
import numpy as np
import pandas as pd
from pathlib import Path
import glob
from scipy.spatial import distance

def build_nestdict(folder):
	#this takes the opened file and builds a dictionary then calls hash to hash it then adds it to nest then
	#exits function, back to main which opens the next file and repeats until no more files.
	nestdictstart=time.time()
	Node=namedtuple=('Node',['key','value'])
	nested_dict={}
	kmerlist=set()
	filecount=0
	table=[None] * 2**20
	files=folder.glob('./*dump.txt')
	for filepath in files:
		file=str(filepath)
		fastafile=file.split('/')[1]
		filename=fastafile.split('.')[0]
		print(filename)
		with open(file) as f:
			kcdict={}
			klist=f.read()
			lines=klist.splitlines()
			filecount+=1
			#dictstart=time.time()
			for line in lines:
				words = line.split()
				kmer = words[0] 
				kmerlist.add(kmer)
				count = int(words[1])##
				kcdict[kmer]=count
			#dictend=time.time()
			#print('Dictionary and set built/added to in: ',dictend-dictstart, "secs")
			#Hash each dictionary/genome, computes hash by dividing dict by size of the table
			index=hash(filename)% len(table)
			#stores hash key at index if empty, if not it chains the hash key
			if table[index] is None:
				table[index] = kcdict
			else:
				table[index]= Node(kcdict, table[index])
		nested_dict[filename]=kcdict
	nestdictend=time.time()
	print("Nested dictionary build and hash indexng done in: ",nestdictend-nestdictstart,"secs")
	return nested_dict, kmerlist,filecount

def set_and_add_zero(nested_dict,kmerlist):
	setstart=time.time()
	for genome,content in nested_dict.items():
		for x in kmerlist - content.keys():
			content[x] = 0
	setend=time.time()
	print("Set comparison to dictionary ran in: ",setend-setstart,"secs")
	totalkmers=len(kmerlist)
	print("Total number of k-mers: ",totalkmers)
	return nested_dict, totalkmers
 
def rel_freq(nested_dict,totalkmers):##int
#Relative frequency calculation and value replacement
	startrf=time.time()
	for genome,content in nested_dict.items():
		for key, value in content.items():
			content[key]= value/totalkmers
	endrf=time.time()
	print("Relative frequencies calculated and replaced in dictionaries in: ", endrf-startrf, "secs")
	return nested_dict 

def jensenshannon_dist(nested_dict):
##JensenShannonDivergence calculation pairwise and addition to nested list
	startcalc=time.time()#
	nesteddict2={}
	for genome, content in nested_dict.items():
		OGgenome=genome
		nesteddict2[OGgenome]={}
		p=np.array(list(content.values()))
		for genome, content in nested_dict.items():
			q=np.array(list(content.values()))
			jsdiv=distance.jensenshannon(p,q)**2#
			nesteddict2[OGgenome][genome]=jsdiv
	endcalc=time.time()
	print("Jensen-Shannon calculations and list builds finished in: ",endcalc-startcalc,"secs")
	return nesteddict2

def matrix_build(nesteddict2):
	startdf=time.time()
	df=pd.DataFrame(nesteddict2)
	enddf=time.time()
	print("Dataframe built in: ",enddf-startdf,"secs")
	return df

def main(inputDir, outputFile):
	start=time.time()
	folder=Path(inputDir)
	#Build dict/table and hash
	nested_dict,kmerlist,filecount=build_nestdict(folder)
	nested_dict,totalkmers=set_and_add_zero(nested_dict,kmerlist)
	#Call relative frequency function on nested dictionary
	nested_dict=rel_freq(nested_dict,totalkmers)
	#Calculate Jensen-Shannon divergence distance between all pairs of genomes and organise into a nested dictionary
	nesteddict2=jensenshannon_dist(nested_dict)
	#Call matrix function to convert nested dictionary of js values into a pandas dataframe
	df=matrix_build(nesteddict2)
	print("Total number of genomes: ",filecount)
	with open(outputFile, "w") as f:
		print(df, file=f)
	#outputFile.close()
	end=time.time()
	print("Program finished in: ",end-start,"secs")


if __name__=='__main__':

	parser = argparse.ArgumentParser(description='K-mer profile table generation')
	parser.add_argument('--inputDir', help='Input a directory containing jellyfish k-mer count files, with the -c column option and dump option, with filenames consisting of X characters and with fasta file endings e.g. NCYC100.fasta', required=True)
	parser.add_argument('--outputFile', help='Distance matrix output file')
	args = parser.parse_args()
	inputDir=args.inputDir
	outputFile=args.outputFile

	main(inputDir, outputFile)

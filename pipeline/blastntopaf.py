#!/usr/bin/env python3

import os
import argparse
import re
			
			
			
def blastntocigar(alignfile, query = ""):

	if len(query):
		with open( query, mode ='r') as f:
			queries = [x.splitlines() for x in f.read().split(">")[1:]]
			qnames = {x[0].split()[0]:"Query_"+str(i+1) for i,x in enumerate(queries)}
			queries = {"Query_"+str(i+1):[x[0].split()[0], len("".join(x[1:]))] for i,x in enumerate(queries)}	


	allaligns = []
	with open(alignfile, mode = 'r') as f:
				
		for line in f:
			
			if len(line) == 0 or line[0]=="@" or line.startswith("[ERROR]"):
				continue
			elements = line.strip().split()
			#name, qsize, qstart, qend, path, rsize, rstart, rend = elements[:8]
			rname, qname = elements[0], elements[2]
			if len(query):
				if qname in qnames:
					qname = qnames[qname]
				elif qname.startswith("BL_ORD_ID:"):
					qname = "Query_"+str(int(qname.split(":")[1])+1)		
				thequeryinfo = queries[qname]
				qname = thequeryinfo[0]
				qfullsize = thequeryinfo[1]
			else:
				qfullsize = 0

			identity =float( [x for i,x in enumerate(elements) if x[:5] == "PI:f:"][0][5:])
			match = int( [x for i,x in enumerate(elements) if x[:5] == "AS:i:"][0][5:] ) 
			
			strand, qstart, cigar = elements[1], int(elements[3])-1, elements[5]
			
			strand = '+' if strand == '0' else '-'
			rstart = [int(x[:-1]) for x in re.findall(r'^\d+H',cigar)]
			if len(rstart):
				rstart = rstart[0]
			else:
				rstart = 0
				
			Hsize =  sum([int(x[:-1]) for x in re.findall(r'\d+H',cigar)])
			rsize = sum([int(x[:-1]) for x in re.findall(r'\d+[MXI]',cigar)])
			rfullsize = rsize + Hsize
			qsize = sum([int(x[:-1]) for x in re.findall(r'\d+[SMXD]',cigar)])
			qend = qsize + qstart

			fullsize = sum([int(x[:-1]) for x in re.findall(r'\d+[A-Z=]',cigar)])
			
			
			rend = rstart + rsize
			if strand == '-':
				temp = rend
				rend = rfullsize  - rstart
				rstart = rfullsize - temp
				
			align = "\t".join(list(map(str,[qname , qfullsize,  qstart,qend,strand ,rname, rfullsize, rstart, rend, match, fullsize-Hsize, 60,identity, cigar])))
				
			print(align)
			
	return allaligns

	
def main(args):
	
	blastntocigar(args.input, args.query)
	
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	
	parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
	parser.add_argument("-q", "--query", help="path to query file", dest="query", type=str,default = "")
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()

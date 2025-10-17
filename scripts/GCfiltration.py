#!/usr/bin/env python3

import collections as cl
import math
import pandas as pd
import argparse
import re
import time
import gc
import sys
import os
def repekmers(kmer):
	
	tmers = cl.defaultdict(int)
	dmers = cl.defaultdict(int)
	lastlastx = ''
	lastx = ''
	for x in kmer:
		tmers[(lastlastx, lastx,x)]+=1
		dmers[(lastx,x)] += 1
		lastlastx = lastx	
		lastx = x
	if max(list(dmers.values())) > 9 or max(list(tmers.values())) > 6:
		return 1
	
	return 0


def maskkmers(infile, cutoff = 0.7):
	
	kmerlist = set()
	masklist = set()
	
	with open(infile, mode = 'r') as r:
		
		text = r.read().splitlines()
		
	for line in text:
		
		if len(line) == 0 or line[0] == '@' or line[0] == '>' :
			
			continue
		line = line.split()[0]	
		#if not re.findall(r'(.)\1{9,}',line):
		if not repekmers(line):
			gc = line.count("C") + line.count("G")
			if gc < cutoff*31 and gc > (1-cutoff)*31:
				
				kmerlist.add(kmerencode(line.strip()))
			else:
				masklist.add(kmerencode(line.strip()))
				
		
	return kmerlist, masklist



def filterkmers(infile, cutoff = 0.7):
	
	kmerlist = set()
	
	
	with open(infile, mode = 'r') as r:
		
		text = r.read().splitlines()
		
	for line in text:
		
		if len(line) == 0 or line[0] == '@' or line[0] == '>' :
			
			continue
		line = line.split()[0]	
		#if not re.findall(r'(.)\1{9,}',line):
		if not repekmers(line):
			gc = line.count("C") + line.count("G")
			if gc < cutoff*31 and gc > (1-cutoff)*31:
				
				kmerlist.add(kmerencode(line.strip()))
				
	return kmerlist

	
def makereverse(seq):
	
	tran=str.maketrans('ATCGatcg', 'TAGCtagc')
	
	return seq[::-1].translate(tran)

def kmerencode(kmer):
	
	kmerint = 0
	for base in kmer:
		
		kmerint *= 4
		if base=='a' or base =='A':
			
			kmerint += 0
			
		elif base=='t' or base =='T':
			
			kmerint += 3
			
		elif base=='c' or base =='C':
			
			kmerint += 1
		elif base=='g' or base =='G':
			
			kmerint += 2
		else:
			raise ValueError("non ATCG base")
			
	return kmerint

def kmerdecode(kmerint, size = 31):
	
	index = 0
	text = ['A' for a in range(size)]
	
	while kmerint:
		
		text[index] = "ACGT"[kmerint % 4]
		kmerint //= 4
		index += 1
		
	return "".join(text[::-1])


class GC_filter:
	
	def __init__(self, kmerlist  ,kmersize = 31, cutoff = 0.7, window = 100, ifmask = 0):
		
		self.cutoff = cutoff
		self.kmersize = kmersize
		self.window_size = window
		self.kmerlist = kmerlist
		self.maskedkmers = set()
		self.intoperator = 4**(kmersize-1)
		
		self.window = ""
		
		self.current_k = 0
		self.reverse_k = 0
		
		self.current_size = 0
		self.current_wsize = 0
		
		self.current_gc = 0
		
		self.translator = str.maketrans('ATCGatcg', 'TAGCtagc')
		
		
	def filterbyGC(self, ref):
		
		currentline = ""
		
		with open(ref,mode = 'r') as r:
			
			for line in r:
				
				
				line = line.strip()
				
				
				if len(line) ==0 or line[0]==">":
					
					self.current_k=0
					self.reverse_k=0
					self.current_size = 0
					self.current_wsize = 0
					self.current_gc = 0
					self.window = ""
					
					continue
				
				for char in line:
					
					char = char.upper()
					
					if char not in ['A','T','C','G']:
						
						self.current_k = 0
						self.reverse_k = 0
						self.current_size = 0
						self.current_wsize = 0
						self.current_gc = 0 
						self.window = "" 
						return
					
					
					if self.current_size >= self.kmersize:
						
						self.current_k %= self.intoperator						
						self.current_k <<= 2
						self.current_k += ['A','C','G','T'].index(char)
						
						self.reverse_k // 4
						self.reverse_k += (['A','C','G','T'].index(char)) * self.intoperator
						
					else:
						
						
						self.current_k <<= 2
						self.current_k += ['A','C','G','T'].index(char)
						self.reverse_k += (['A','C','G','T'].index(char)) * (4**self.current_size)
						
						self.current_size += 1
						
						
					if  self.current_wsize >= self.window_size:
						
						self.window = self.window[1:] + char
						
						if self.reverse_k%4 in [1,2]:
							
							self.current_gc-=1
							
					else:
						
						self.current_wsize += 1
						self.window += char
						
						
					if char in ['G','C']:
						
						self.current_gc+=1
						
						
					kmer = max(self.current_k, self.reverse_k)
					
					if self.current_size == self.kmersize and kmer in self.kmerlist:
						
						if self.current_gc < self.window_size*self.cutoff and self.current_gc > self.window_size*(1-self.cutoff):
							
							continue
						
						self.kmerlist.remove(kmer)
						
						
						
						
						
						
def main(args):
	
	
	infile = args.input
	
	outfile = args.output
	
	reffile = args.ref
	
	masked = set()
	if args.mask:
		kmerlist, masklist = maskkmers(infile,args.cutoff)
	else:
		kmerlist = filterkmers(infile,args.cutoff)
		
	gc.collect()
	time.sleep(1)
	
	data = GC_filter(kmerlist,cutoff = args.cutoff)
	if len(reffile):
		data.filterbyGC(reffile)
		if os.path.isfile(reffile+"._app_"):
			data.filterbyGC(reffile+"._app_")
			
	with open (outfile, mode = 'w') as w:
		
		for kmer in data.kmerlist:
			w.write(">\n"+kmerdecode(kmer)+"\n")
			
		if args.mask:
			masklist.update(kmerlist - data.kmerlist)
			for kmer in masklist:
				w.write(">\n" + kmerdecode(kmer).lower() + "\n")
				
	#print(",".join(list(map(lambda x: str(len(x)) if 1==1 else "0.0" if len(x)==0 else str(sum(x )/(len(x))),list(samples.values())))))
				
				
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program filter kmers by GC and repeats")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	parser.add_argument("-r", "--ref", help="path to reference file", dest="ref", type=str, default = "")
	parser.add_argument("-c", "--cutoff", help="GC filtration cutoff, bothside", dest="cutoff", type=float, default = 0.8)
	parser.add_argument("-m", "--mask", action="store_true",help="masking instead of filtering kmers")
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()

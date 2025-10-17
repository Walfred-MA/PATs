#!/usr/bin/env python3

import re
import pandas as pd
import os
import argparse
import sys
import collections as cl

class MaskChecker:
	
	def __init__(self, thefile):
	

		self.maskregion = cl.defaultdict(list)


		"""		
		with open(thefile, mode = 'r') as f:
			
			reads = [read.split() for read in f.read().splitlines()]
		

		for read in reads:
		
			self.maskregion[read[0]].append([int(read[-2]), int(read[-1])])
		"""

		with open(thefile, mode = 'r') as f:
			reads = f.read().split(">")[1:]
		self.names = [seq[0].split()[0] for seq in reads]
		reads = {seq[0].split()[0]:"".join(seq[1:]) for seq in reads}
		for name, seq in reads.items():
			
			self.maskregion[name] = [x.span() for x in re.finditer(r'[a-z]+',seq)]

		
	def MaskSize(self, name, region):
		
		maskedregion = [x for y in self.maskregion[name] for x in y] 
				
		addstart =0
		addend = 0
		
		overlap = []
		for index, coordi in enumerate(maskedregion):
			
			if coordi < region[0]:
				
				addend = 1
				
				continue
			
			
			if coordi >= region[1]:
				
				addstart = 1
				
				continue
			
			overlap.append((index,coordi))
									
		if len(overlap) == 0:
			
			if addstart and addend:
				
				overlap = [(0,region[0]),(0,region[1])]
			
			else:
		
				return 0 
			
		else:
			
			if overlap[0][0] % 2:
			
				overlap.insert(0, (0,region[0]))
		
			if overlap[-1][0] % 2 == 0:
			
				overlap.append((0,region[1]))
			
		
		overlapsize = -sum([overlap[index*2][1] - overlap[index*2+1][1]  for index in range(len(overlap)//2)])
				
		return overlapsize
			

def run(args):
	
	maskchecker = MaskChecker(args.ref)
	
	table = pd.read_csv(args.input, header = None, sep = "\t", usecols = range(14) )
	
	passfiler = []
	for index, row in table.iterrows():

		query = row.iat[0]

		if  query.count("_") > 1 and query.split("_")[-1].isdigit():
			start = int(query.split("_")[-2]	)
			table.at[index,0] = "_".join(query.split("_")[:-2])
			table.at[index,2] += start
			table.at[index,3] += start 

		name = row.iat[5]
		if name.startswith("BL_ORD_ID:"):
			name = maskchecker.names[int(name.split(":")[1])] 
		region = row[7:9].tolist()
		
		length = row.iat[6]

		qual = row.iat[11]
	
		simi = row.iat[12]	
		match = row.iat[9]
		alignsize = row.iat[10]
		mismatch = alignsize - match
		
		unmasksize = region[1] - region[0] - maskchecker.MaskSize(name, region)
	
		if ( unmasksize > min(length/3, 1000) or (region[1] - region[0]) > min(length/3, 3000) ) and simi >= 90 :
	
			passfiler.append(index)
	
	table_filter = table.iloc[passfiler]
	
	table_filter.to_csv(args.output,index=False, header=None, sep = "\t", mode='w')
	


def main():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",
						type=str)
	parser.add_argument("-r", "--ref", help="group infor", dest="ref",
										type=str)
	
	parser.set_defaults(func=run)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	
	main()

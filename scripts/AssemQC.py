#!/usr/bin/env python3


import os
import sys
import argparse
import collections as cl

mainchr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

class error_region:
	def __init__(self, loadfile, extend = 25000):
		
		self.extend = extend
		self.error_regions = cl.defaultdict(list)
		
		if len(loadfile) == 0:
			return
		
		with open(loadfile,mode = 'r') as f:
			
			for line in f:
				
				if len(line.strip())==0 or line[0] in ["@","#","=","\t"," "]:
					continue
				contig, start, end = line.split()[:3]
				self.error_regions[contig].append(sorted([int(start),int(end)]))
	
	def add(self, contig, region):
		
		self.error_regions[contig].append([region[0] - self.extend, region[1] + self.extend])
		
		
	def findoverlap(self, contig, region):
		
		error_oncontig = self.error_regions[contig]
		
		if len(error_oncontig) == 0:
			return False
		
		region_size = region[1]-region[0]
		for error_region in error_oncontig:
			
			allcordi = error_region + region
			
			overlap =   max(allcordi) - min(allcordi)  -   (error_region[1]-error_region[0]) - region_size
			if overlap <= 0:
				
				return True
			
		return False

def extendfilter(allcontigs,error_regions, ifextend):
	
	filtered = set()
	for contig, regions in allcontigs.items():
		
		keepfilter = 1
		while keepfilter:
			keepfilter = 0
			for region in regions:
				if error_regions.findoverlap(contig, region[:2]):
					if region[2] in filtered:
						continue
					keepfilter = 1
					filtered.add(region[2])
					error_regions.add(contig, region[:2])
					break
			if not ifextend:
				break
				
	return filtered

def getminedge(inputfile,scafsizes, whitelist):
	
	min_edge = 1000000000
	with open(inputfile, mode = 'r') as f:
		for line in f:
			if len(line) ==0 or line[0] != ">":
				continue
			header = line[1:]
			region = header.split()[1]
			contig = region.split(":")[0]
			haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1" if "#" not in contig else "_h".join(contig.split("#")[:2])
			if haplo in whitelist or contig in mainchr:
				coordi = region.split(":")[1]
				if coordi.endswith("-") or coordi.endswith("+"):
					coordi = coordi[:-1]
				coordi = coordi.split("-")
				coordi = ["-".join(coordi[:-1]),coordi[-1]]
				coordi = sorted([int(x) for x in coordi])
				
				scafsize = scafsizes.get(contig, 1000000000)
				
				edge = max(1000,min(coordi[0], scafsize - coordi[1]))
				min_edge = min(min_edge, edge)
	return min_edge


def runfilter(input, output, scafsizes, blacklist,error_regions, cutoff, blackcutoff, ifextend):
	
	
	Ngaps = cl.defaultdict(int)
	name = ""
	with open(input, mode = 'r') as f:
		for line in f:
			if line.startswith( ">"):
				name = line.split()[0][1:]
			Ngaps[name] += (line.count("N") + line.count("n"))
			
			
	passnum = 0
	totalnum = 0
	iffilter = 0
	filterlist = set()
	allcontigs = cl.defaultdict(list)
	with open(input, mode = 'r') as f:
		
		for line in f:
			if len(line.strip()) ==0:
				continue
			if line[0] == ">":
				
				iffilter = 1
				header = line[1:]
				name = header.split()[0]
				region = header.split()[1]
				contig = region.split(":")[0]
				
				haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1" if "#" not in contig else "_h".join(contig.split("#")[:2])
				
				coordi = region.split(":")[1]
				if coordi.endswith("-") or coordi.endswith("+"):
					coordi = coordi[:-1] 
				coordi = [coordi.split("-")[-1], "-".join(coordi.split("-")[:-1])]
				coordi = sorted([int(x) for x in coordi])
				
				scafsize = scafsizes.get(contig, 1000000000) 
				edge = min(coordi[0], scafsize - coordi[1]) 
				size = coordi[1]-coordi[0]
				totalnum += 1
				
				allcontigs[contig].append([coordi[0],coordi[1],name])
				
				if Ngaps[name] > 100 or (haplo in blacklist and edge < blackcutoff) or (haplo not in blacklist and edge < cutoff) or error_regions.findoverlap(contig,coordi) :
					iffilter = 1
					if ifextend:
						error_regions.add(contig,coordi)
					filterlist.add(name)
				else:
					iffilter = 0
					passnum += 1

	filtered = extendfilter(allcontigs,error_regions,ifextend)
	filterlist.update(filtered)
	with open(input,mode = 'r') as r, open(output,mode = 'w') as w, open(output+"_QCfiltered.fa",mode = 'w') as w2:
		
		for line in r:
			
			if line.startswith(">"):
				
				name= line.split()[0][1:]
				iffilter = 0
				if name in filterlist:
					iffilter = 1
					
			if not iffilter:
				w.write(line)
			else:
				w2.write(line)
	
	
	return passnum , totalnum 


def main(args):
	
	whitelist = ["CHM13_h1","CN00001_h1","CN00001_h2", "CN00002_h1", "CN00002_h2", "HG002_h1","HG002_h2"]
	blacklist = ['apr001_h1','apr001_h2','apr002_h1','apr002_h2','apr003_h1','apr003_h2','apr004_h1','apr004_h2','apr005_h1','apr005_h2','apr006_h1','apr006_h2','apr007_h1','apr007_h2','apr008_h1','apr008_h2','apr009_h1','apr009_h2','apr010_h1','apr010_h2','apr011_h1','apr011_h2','apr012_h1','apr012_h2','apr013_h1','apr013_h2','apr014_h1','apr014_h2','apr015_h1','apr015_h2','apr016_h1','apr016_h2','apr018_h1','apr018_h2','apr019_h1','apr019_h2','apr020_h1','apr020_h2','apr021_h1','apr021_h2','apr022_h1','apr022_h2','apr023_h1','apr023_h2','apr024_h1','apr024_h2','apr025_h1','apr025_h2','apr026_h1','apr026_h2','apr027_h1','apr027_h2','apr029_h1','apr029_h2','apr030_h1','apr030_h2','apr031_h1','apr031_h2','apr032_h1','apr032_h2','apr033_h1','apr033_h2','apr034_h1','apr034_h2','apr035_h1','apr035_h2','apr036_h1','apr036_h2','apr037_h1','apr037_h2','apr038_h1','apr038_h2','apr039_h1','apr039_h2','apr040_h1','apr040_h2','apr041_h1','apr041_h2','apr042_h1','apr042_h2','apr043_h1','apr043_h2','apr044_h1','apr044_h2','apr045_h1','apr045_h2','apr046_h1','apr046_h2','apr047_h1','apr047_h2','apr048_h1','apr048_h2','apr049_h1','apr049_h2','apr050_h1','apr050_h2','apr051_h1','apr051_h2','apr52_h1','apr52_h2','aprf_h1','aprf_h2','aprm_h1','aprm_h2','aprs_h1','aprs_h2','GW00001_h1','GW00001_h2','GW00002_h1','GW00002_h2','GW00003_h1','GW00003_h2','GW00004_h1','GW00004_h2','GW00005_h1','GW00005_h2','GW00006_h1','GW00006_h2','GW00007_h1','GW00007_h2','GW00008_h1','GW00008_h2','GW00009_h1','GW00009_h2','GW00010_h1','GW00010_h2','GW00011_h1','GW00011_h2','GW00012_h1','GW00012_h2','GW00013_h1','GW00013_h2','GW00014_h1','GW00014_h2','GW00015_h1','GW00015_h2','GW00016_h1','GW00016_h2','GW00017_h1','GW00017_h2','GW00018_h1','GW00018_h2','GW00019_h1','GW00019_h2','GW00020_h1','GW00020_h2','GW00021_h1','GW00021_h2','GW00022_h1','GW00022_h2','GW00023_h1','GW00023_h2','GW00024_h1','GW00024_h2','GW00025_h1','GW00025_h2','GW00026_h1','GW00026_h2','GW00027_h1','GW00027_h2','GW00028_h1','GW00028_h2','GW00029_h1','GW00029_h2','GW00030_h1','GW00030_h2','GW00031_h1','GW00031_h2','GW00032_h1','GW00032_h2','GW00033_h1','GW00033_h2','GW00034_h1','GW00034_h2','GW00035_h1','GW00035_h2','GW00036_h1','GW00036_h2','GW00037_h1','GW00037_h2','GW00038_h1','GW00038_h2','GW00039_h1','GW00039_h2','GW00040_h1','GW00040_h2','GW00041_h1','GW00041_h2','GW00042_h1','GW00042_h2','GW00043_h1','GW00043_h2','GW00044_h1','GW00044_h2','GW00045_h1','GW00045_h2','GW00046_h1','GW00046_h2','GW00047_h1','GW00047_h2','GW00048_h1','GW00048_h2','GW00049_h1','GW00049_h2','GW00050_h1','GW00050_h2','GW00051_h1','GW00051_h2','GW00052_h1','GW00052_h2','GW00053_h1','GW00053_h2','GW00054_h1','GW00054_h2','GW00055_h1','GW00055_h2','GW00056_h1','GW00056_h2','GW00057_h1','GW00057_h2','HG03125_h1','HG03125_h2','NA12878_h1','NA12878_h2','NA18990_h1','NA18990_h2','NA19007_h1','NA19007_h2','NA19066_h1','NA19066_h2','NA19079_h1','NA19079_h2','NA19088_h1','NA19088_h2','ksa001_h1','ksa001_h2','ksa002_h1','ksa002_h2','ksa003_h1','ksa003_h2','ksa004_h1','ksa004_h2','ksa005_h1','ksa005_h2','ksa006_h1','ksa006_h2','ksa007_h1','ksa007_h2','ksa008_h1','ksa008_h2','ksa009_h1','ksa009_h2','NA12889_h1','NA12889_h2','NA12890_h1','NA12890_h2','NA12891_h1','NA12891_h2','NA12892_h1','NA12892_h2']
	
	scafsizes = {}
	if len(args.scaf):
		with open(args.scaf, mode = 'r') as f: 
				for info in f:
					if len(info):
						info = info.strip().split()
						scafsizes[info[0]] = int(info[1])
						
	if len(args.ref):
		refs = args.ref.split(",")
		for ref in refs:
			with open(ref+".fai", mode = 'r') as f:
				for info in f:
					if len(info):
						info = info.strip().split()
						scafsizes[info[0]] = int(info[1])
						
	if len(args.query):
		with open(args.query, mode = 'r') as f:
			for line in f:
				line = line.split()
				
				with open(line[1]+".fai", mode = 'r') as f:
					for info in f:
						if len(info):
							info = info.strip().split()
							scafsizes[info[0]] = int(info[1])
							
	error_regions = error_region(args.error)
	
	minedge = getminedge(args.input, scafsizes, whitelist)
	runfilter(args.input, args.output, scafsizes, blacklist, error_regions, min(25000, minedge//2), 25000, args.ifextend)
	
	
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to input data file", dest="output", type=str, required=True)
	parser.add_argument("-e", "--error", help="path to input data file",dest="error", type=str,default = "")
	parser.add_argument("-x", "--ifextend", help="path to input data file",dest="ifextend", type=int,default = 1)
	parser.add_argument("-q", "--querypath", help="path to input data file",dest="query", type=str, default = "")
	parser.add_argument("-s", "--scaf", help="path to input data file",dest="scaf", type=str, default = "")
	parser.add_argument("-c", "--cut", help="path to scaffold file",dest="cut", type=int, default = 25000)
	parser.add_argument("-r", "--ref", help="path to ref file",dest="ref", type=str, default = "")
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
	

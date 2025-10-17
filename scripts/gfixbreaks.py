#!/usr/bin/env python3

import os
import threading
import argparse
import concurrent.futures
import datetime
import collections as cl
import re
import glob
import json
import hashlib 
from statistics import median

script_folder = os.path.dirname(os.path.abspath(__file__))
cutoffdistance = 30000
cache = ""

def process_locus2(region, haplopath, folder, outputfile):
	
	newname,  contig, start, end, ifexon, mapinfo = region
	
	if "#" not in contig: 
		haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
	else:
		haplo = "_h".join(contig.split("#")[:2])
		
	path = haplopath[haplo]
	
	outputfile0 = os.path.join(folder, f"{newname}.fa")
	start = max(0,int(start))
	header = ">{}\t{}:{}-{}\t{}\t{}".format(newname, contig, start, end, ifexon, mapinfo)
	
	cmd1 = f'echo "{header}" >> {outputfile} && samtools faidx {path} {contig}:{start+1}-{end} | tail -n +2 >> {outputfile}'
	os.system(cmd1)
	
	return outputfile0


def makenewfasta2(regions, queryfile, outputfile, folder, threads):
	
	haplopath = dict()
	with open(queryfile, mode = 'r') as f:
		for line in f:
			line = line.strip().split()
			haplopath[line[0]] = line[1]
			
	os.system("rm {} || true ".format(outputfile))
	
	alloutputs = []
	
		#with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
			#futures = []
	for region in regions:
		
			# Submit jobs to executor
		process_locus2( region, haplopath, folder, outputfile)
			#futures.append(executor.submit(process_locus2, region, haplopath, folder, outputfile))
		#alloutputs.append(os.path.join(folder, f"{region[0]}.fa"))
		# Wait for all threads to complete
		#concurrent.futures.wait(futures)
		
def makegraph(fastafile, folder, threads):
	
	#os.system(f"{script_folder}/kmernorm -i {fastafile} -o {fastafile}_norm.gz -w 0 -m 1")
	#os.system(f"rm {fastafile}_part*.fa || true");
	#os.system(f"python {script_folder}/normpartition.py -i {fastafile} -n {fastafile}_norm.gz -o {fastafile}_part")
	
	#files = glob.glob(f"{fastafile}_part*.fa")
	
	files = []
	if not len(files):
		files = [fastafile]
		
	for file in files:
		pass
		os.system("python {}/graphmake.py -i {}  -d {} -t {} -l 0 ".format(script_folder, file, folder, threads))
		
	return files

def process_locus(index, line, haplopath, folder, outputfile):
	line = line.strip().split()
	contig = line[0]
	
	if "#" not in contig: 
		haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
	else:
		haplo = "_h".join(contig.split("#")[:2])
		
	path = haplopath[haplo]
	strd = "-i" if line[1] == "-" else ""
	
	outputfile0 = os.path.join(folder, f"loci_{index}.fa")
	
	line[2] = max(0,int(line[2]))
	header = ">loci_{}\t{}:{}-{}{}\t{}".format(index, contig, line[2], line[3], line[1], line[4])
	
	if os.path.isfile(path):
		cmd1 = f'echo "{header}" > {outputfile0} && samtools faidx {path} {contig}:{line[2]+1}-{line[3]} {strd} | tail -n +2 >> {outputfile0}'
		os.system(cmd1)
		
	return outputfile0


def makenewfasta(locifile, queryfile, outputfile, folder, threads):
	
	haplopath = dict()
	with open(queryfile, mode = 'r') as f:
		for line in f:
			line = line.strip().split()
			haplopath[line[0]] = line[1]
			
	os.system("echo > "+outputfile)
	
	alloutputs = []
	with open(locifile, mode='r') as f:
		with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
			futures = []
			for index, line in enumerate(f, start=1):
				# Submit jobs to executor
				#process_locus(index, line, haplopath, folder, outputfile)
				futures.append(executor.submit(process_locus, index, line, haplopath, folder, outputfile))
				alloutputs.append(os.path.join(folder, f"loci_{index}.fa"))
			# Wait for all threads to complete
			concurrent.futures.wait(futures)
			
	for result in alloutputs:
		outputfile0 = result
		cmd2 = f"cat {outputfile0} >> {outputfile}"
		os.system(cmd2)
		
		
def findbreaks(inputfile, assemfile, outputfile):
	
	assemsizes = cl.defaultdict(lambda: 10000000000)
	if len(assemfile):
		with open(assemfile, mode = 'r') as f:
			for line in f:
				path = line.split()[-1]
				with open(path+".fai", mode = 'r') as f:
					for line in f:
						line = line.split()
						assemsizes[line[0]] = int(line[1])
				
	ifexons = set()
	names = []
	contigs = cl.defaultdict(list)
	name_tocontig = cl.defaultdict(list)
	index = -1
	with open(inputfile, mode = 'r') as f:
		
		for line in f:
			
			if line.startswith(">"):
				
				index += 1
				name, locus = line[1:].split()[:2]
				names.append(name)
				exon = line.split()[-2]
				if exon == "Exon":
					ifexons.add(name)
					
				strd = "+"
				contig, region = locus.split(":")
				if region[-1] in ["+","-"]:
					region = region[:-1]
					strd = locus[-1]
					
				region = region.split("-")
				region = [int(region[0]), int(region[1]), name] if strd == "+" else [int(region[1]), int(region[0]), name]
				
				contigs[contig].append(region)
				
	allsize =0 
	new_contigs = cl.defaultdict(list)
	for contig, regions in contigs.items():
		
		new_regions = []
		
		strd = sorted(regions, key = lambda x:  -abs(x[0]-x[1]))[0]
		
		strd = '+' if strd[1] > strd[0] else '-'
		
		regions_sort = sorted(regions) if contig[-1] == '+' else sorted(regions, reverse = 1)
		
		lastend = regions_sort[0]
		lastregion = regions_sort[0]
		new_regions = [[regions_sort[0]]]
		
		allsize += abs(regions_sort[0][1]-regions_sort[0][1])
		for region in regions_sort[1:]:
			if max(region[0],region[1], lastend[0],lastend[1]) - min(region[0],region[1], lastend[0],lastend[1]) - abs(region[1]-region[0] ) - abs(lastend[1] - lastend[0]) < cutoffdistance:
				
				new_regions[-1].append(region)
			else:
				new_regions.append([region])
				
			allsize += abs(region[0]-region[1])
			lastend = region
			
		new_contigs.update({(contig+strd,i):regions for i,regions in enumerate(new_regions)})
		
		
	blacklist = set()
	with open(outputfile, mode = 'w') as w:
		index = 1
		for (scarf, index), regions in new_contigs.items():
			locus = "{}".format(scarf[:-1])
			start = min(regions[0][0],regions[0][1],regions[-1][0],regions[-1][1])
			end = max(regions[0][0],regions[0][1],regions[-1][0],regions[-1][1])
			names = ";".join([region[-1] for region in regions])
			
			start = max(0,start-cutoffdistance)
			end = min(end+cutoffdistance, assemsizes[scarf[-1]])
			if start == 0 or end == assemsizes[scarf[:-1]]:
				blacklist.add(f"loci_{index}")
				
			w.write("{}\t{}\t{}\t{}\t{}\n".format(locus,scarf[-1],start,end,names))
			
			index += 1
			
	return blacklist

def polishranges(ranges):
	
	if len(ranges) <= 1:
		return [ranges]
	
	size = abs(ranges[0][2] - ranges[0][1])
	lastend = ranges[0][2]
	breaks = set()
	for i,(index,start,end) in enumerate(ranges[1:]):
		gap = abs(start-lastend)
		if gap > size:
			breaks.add(i+1)
			size = end - start
		else:
			size += end - start - gap
		lastend = end
		
	size = abs(ranges[-1][2] - ranges[-1][1])
	lastend = ranges[-1][1]
	
	for i,(index, start,end) in enumerate(ranges[::-1][1:]):
		gap = abs(end-lastend)
		if gap > size:
			breaks.add(len(ranges) - 2 - i)
			size = end - start
		else:
			size += end - start - gap
		lastend = start
		
	eachbreak = []
	allbreaks = []
	for i,x in enumerate(ranges):
		if i in breaks:
			allbreaks.append(eachbreak)
			eachbreak = []
		eachbreak.append(x)
	allbreaks.append(eachbreak)
	
	return allbreaks

def Overlap(allaligns, allow_gap = 300):
	
	all_coordinates = [x for y in allaligns for x in y[:2]]
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	allgroups = []
	current_group = set([])
	current_group_size = 0 
	current_group_start = 0
	
	number_curr_seq = 0
	coordinate = 0
	current_group = [0,0]
	last_coordinate = -allow_gap-1
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		
		if index %2 == 0 :
			
			number_curr_seq += 1
			
			if number_curr_seq == 1:
				
				current_group_start = coordinate
				
				if coordinate - last_coordinate> allow_gap :
					
					allgroups.append(current_group)
					
					current_group = [coordinate,coordinate]					
		else:
			number_curr_seq -= 1
			current_group[1] = coordinate
			
		last_coordinate = coordinate
		
	allgroups.append(current_group)
	
	allgroups = [x for x in allgroups if x[1] - x[0] > 100]
	
	return allgroups

def qposi_togposi(cigar_str, find_qposes):
	
	l = len(find_qposes)
	if l == 0:
		
		return []
	
	allcigars = re.findall(r'[><][^><]+', cigar_str)
	
	curr_rposi = 0
	curr_qposi = 0
	new_rposi = 0
	new_qposi = 0
	
	curr_strd = 1
	curr_strd_str = ">"
	curr_size = 0
	
	find_qposes = find_qposes + [0]
	find_qpos_index = 0
	find_qpos = find_qposes[0]
	find_rposes = []
	
	for cigars in allcigars:
		
		path, cigar  = cigars.split(":")
		
		if len(find_rposes):
			
			find_rposes.append((pathname,curr_rposi,curr_qposi-1,curr_strd_str ))
			
		pathname = "_".join(path.split("_")[:-2])[1:]
		
		curr_strd_str = cigars [0]
		curr_strd = 1 if curr_strd_str  =='>' else -1
		
		curr_rposi,curr_qposi = path.split("_")[-2:]
		
		curr_rposi,curr_qposi = int(curr_rposi), int(curr_qposi)
		
		if len(find_rposes):
			find_rposes.append((pathname,curr_rposi,curr_qposi,curr_strd_str ))
			
		cigars = re.findall(r'\d+[M=XIHD]', cigars)
		
		if cigars[0][-1] == 'H':
			cigars = cigars[1:]
		if cigars[-1][-1] == 'H':
			cigars = cigars[:-1]
			
		new_rposi, new_qposi = curr_rposi,curr_qposi
		
		for cigar in cigars:
			
			curr_size = int(cigar[:-1])
			char = cigar[-1]
			
			if char == '=' or char ==  'M':
				
				new_rposi = curr_rposi + curr_strd * curr_size
				new_qposi = curr_qposi + curr_size
				
			elif char == 'I':
				
				new_qposi = curr_qposi + curr_size 
				
			elif char == 'D' or char == 'H':
				
				new_rposi = curr_rposi + curr_strd * curr_size
				
			elif char == 'X':
				
				new_rposi = curr_rposi + curr_strd * curr_size
				new_qposi = curr_qposi + curr_size
				
			elif char == 'S':
				
				if curr_qposi == 0:
					curr_qposi = curr_size
					
				curr_size = 0
				continue
			
			while find_qpos_index<l and find_qpos <= new_qposi:
				
				if char == '=' or char ==  'M' or char =='X' :
					
					find_rposes.append((pathname, curr_rposi + curr_strd * (find_qpos-curr_qposi), find_qpos, curr_strd_str ))
				else:
					
					find_rposes.append((pathname,curr_rposi, find_qpos, curr_strd_str ))
					
				find_qpos_index += 1
				find_qpos = find_qposes[find_qpos_index]
				
			if find_qpos_index == l:
				
				return find_rposes
			
			curr_size = 0
			curr_rposi = new_rposi
			curr_qposi = new_qposi
			
			
	while find_qpos_index<l:
		
		find_qpos_index += 1
		find_rposes.append((pathname, curr_rposi, find_qpos, curr_strd_str ))
		
	return find_rposes


def gposi_toqposi(cigar_str, find_rposes):
	
	l = len(find_rposes)
	if l == 0:
		
		return []
	
	find_rposes_dict = cl.defaultdict(list)
	
	cigars = re.findall(r'\d+[M=XHDI]', cigar_str)
	
	curr_qposi = 0
	new_qposi = 0
	curr_rposi = 0
	new_rposi = 0
	
	find_rposes = find_rposes+[-1]
	
	find_rpos_index = 0
	find_rpos = find_rposes[find_rpos_index] 
	curr_strd = 1
	find_qposes = []
	for cigar in cigars:
		
		curr_size = int(cigar[:-1])
		char = cigar[-1]
		
		if char == '=' or char ==  'M':
			
			new_rposi = curr_rposi + curr_size 
			new_qposi = curr_qposi + curr_size * curr_strd
		elif char == 'I':
			
			new_qposi = curr_qposi + curr_size * curr_strd
			
		elif char == 'D' or char == 'H':
			
			new_rposi = curr_rposi + curr_size 
			
		elif char == 'X':
			
			new_rposi = curr_rposi + curr_size 
			new_qposi = curr_qposi + curr_size * curr_strd
			
			
		elif char == 'S':
			
			if curr_qposi == 0:
				curr_qposi = curr_size 
				
			curr_size = 0
			continue
		
		while find_rpos_index<l and find_rpos <= new_rposi + (1-curr_strd)-1//2 :
			
			
			if char == '=' or char ==  'M' or char =='X' :
				
				find_qposes.append(curr_qposi + curr_strd * (find_rpos-curr_rposi))
			else:
				
				find_qposes.append(curr_qposi)
				
			find_rpos_index += 1
			find_rpos = find_rposes[find_rpos_index]
			
		if find_rpos_index == l:
			return find_qposes
		
		curr_size = 0
		curr_rposi = new_rposi
		curr_qposi = new_qposi
		
		
	while find_rpos_index<l:
		
		find_rpos_index += 1
		find_qposes.append(curr_qposi)
		
		
	return find_qposes


def getspan_ongraph(name, path, cigar, rposis, qposis):
	
	cigar = cigar.translate(str.maketrans('', '', 'ATCGatcgNn-'))
	
	cigars = re.findall(r'[><][^><]+', cigar)
	pathes = re.findall(r'[><][^<>]+',path)
	
	
	rposis = [ [int(x.split("_")[0]),int(x.split("_")[1])] for x in rposis.split(";")]
	
	qposis = [ [int(x.split("_")[0]),int(x.split("_")[1])] for x in qposis.split(";")]
	
	newcigars= "".join(["{}_{}_{}:{}".format(thename,rposi[0],qposi[0],cigar.split(":")[1]) if thename[0] == '>' else "{}_{}_{}:{}".format(thename,rposi[1],qposi[0],cigar.split(":")[1]) for thename,cigar,rposi,qposi in zip(pathes,cigars,rposis, qposis)])
	
	span = list(zip(pathes, [tuple([int(a) for a in x]) for x in rposis], [tuple([int(a) for a in x]) for x in qposis]))
	
	span = [x if x[0][0] =='>' else tuple([x[0],x[1],sorted(x[2],reverse = 1)])  for x in span]
	
	
	return newcigars, span


def readgraph(graphfile):
	
	contigspan = dict()
	pathsize = dict()
	graphinfo = dict()
	with open(graphfile, mode = 'r') as f:
		
		for line in f:
			
			line = line.strip().split()
			
			if len(line) < 5:
				continue
			
			name, path, cigar, rposis, qposis = line 
			
			graphinfo[name],contigspan[name]  = getspan_ongraph(name, path, cigar, rposis, qposis)
			
			
	return graphinfo, contigspan



def readcontigsinfo(inputfile, breakfile, refprefix = "NC_0609"):
	
	genetoscaff = dict()
	with open(inputfile, mode = 'r') as f:
		for line in f:
			if line.startswith(">"):
				line = line.strip().split() + ["",""]
				name = line[0][1:]
				locus = line[1]
				strd = locus[-1]
				contig, coordi = locus.split(":")
				if coordi[-1] in ['-','+']:
					coordi = coordi[:-1]
				start, end = map(int, coordi.split("-"))
				start, end = min(start,end), max(start,end)
				start, end = max(0, start),  end
				genetoscaff[name] = [contig, strd, start, end, line[2], line[3]]
				
				
	locuslocations = dict()
	genetolocus = cl.defaultdict(list)
	with open(breakfile, mode = 'r') as f:
		
		scarf_index = 0
		lastend = -30000
		laststart = -30000
		lastcontig = ""
		lastgenes = set()
		index = 0
		for index, line in enumerate(f):
			
			index += 1
			
			line = line.strip().split()
			
			locusname = f"loci_{index}" 
			contig, strd, start, end = line[0], line[1], int(line[2]), int(line[3])
			
			if start > lastend or contig != lastcontig:
				scarf_index += 1
				
			lastcontig = contig
			lastend = end
			
			names = line[-1].split(";")
			
			lastgenes = set(names)
			
			if len(names) == 1:
				pass 
				
			if refprefix in contig:
				locusname = "Ref_"+locusname
				
			locuslocations[locusname ] = [contig+"_"+str(scarf_index), strd, start, end]
			
			for name in names:
				if strd == '+':
					qstart = genetoscaff[name][2] - start
					qend = genetoscaff[name][3] - start
				else:
					qstart = end - genetoscaff[name][3]
					qend = end  - genetoscaff[name][2]
					qstart,qend = min(qstart,qend), max(qstart,qend)
					
				genetolocus[name].append([locusname, qstart, qend] )
				
				genetoscaff[name][0] = contig+"_"+str(scarf_index)
				
			genetoscaff[locusname] = [contig+"_"+str(scarf_index), strd, start, end]
			
			
	return genetoscaff, genetolocus, locuslocations

def combinebreaks_eachpath(names, breaks):
	
	all_coordinates = breaks
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	break_groups = []
	break_group = []
	
	break_group_index = set([])
	break_group_indexs = []
	
	last_coordinate = 0
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		
		if coordinate - last_coordinate < 100:
			
			break_group.append(coordinate)
			break_group_index.add(names[index//2])
			
		else:
			
			break_groups.append(break_group)
			break_group_indexs.append(break_group_index)
			
			break_group = [coordinate]
			break_group_index = set([names[index//2]])
			
		last_coordinate = coordinate
		
	break_groups.append(break_group)
	break_group_indexs.append(break_group_index)
	
	extendgroup = cl.defaultdict(set)
	
	break_count = cl.defaultdict(int)
	break_uniformed = dict()
	break_tocontignames = cl.defaultdict(set)
	for group,names in zip(break_groups,break_group_indexs):
		
		if len(group) == 0:
			continue
		themedian = int(median(group))
		break_tocontignames[themedian] = names
		
		for coordi in group:
			break_uniformed[coordi] = themedian
			
	return break_uniformed, break_tocontignames


def combinebreaks(gbreaks):
	
	coordi_onpath = cl.defaultdict(list)
	name_onpath = cl.defaultdict(list)
	for name, breaks in gbreaks.items():
		
		breaks = [ x for y in breaks for x in y]
		
		for (path, gcoordi, pcoordi, strd,  genename) in breaks:
			
			coordi_onpath[path].append(gcoordi)
			name_onpath[path].append(name)
			
	break_uniformed_new = dict()
	for path, breaks in coordi_onpath.items():
		
		break_uniformed, break_tocontignames = combinebreaks_eachpath(name_onpath[path], breaks)
		
		for posi0, posi1 in break_uniformed.items():
			
			break_uniformed_new[(path,posi0)] = posi1
			
	return break_uniformed_new



def findvalidbreaks(break_contignames, contigspans):
	
	allcontigs_spans = [z  for y in contigspans.values() for x in y for z in [x[0]-10,x[1]+10]]
	
	allnames =  [z  for name, y in contigspans.items() for x in y for z in [name,name]] 
	
	l1 =  len(allcontigs_spans)
	
	allcoordis = allcontigs_spans+ list(break_contignames.keys())
	
	allcoordis_sort = sorted(list(range(len(allcoordis))), key = lambda x: allcoordis[x])
	
	locus_span = 0
	segment_span = 0
	
	locus_lastposi = cl.defaultdict(int)
	
	allbreaks_found = []
	current_names = set()
	lastbreak = 0
	
	breaks_allspancontigs = cl.defaultdict(list)
	for sortindex, index in enumerate(allcoordis_sort):
		
		if index < l1:
			if index % 2 == 0:
				segment_span += 1
				current_names.add(allnames[index])
			else:
				segment_span -= 1
				try:
					current_names.remove(allnames[index])
				except:
					pass
					
		else:
			coordi = allcoordis[index]
			common = break_contignames[coordi] & set(current_names)
			uncommon = set(current_names) ^ common
			breaks_allspancontigs[coordi].extend(current_names)
			
			totalbreak = len(common)
			totalcontigs = len(current_names)
			
			ifref1,ifref2 = any(x for x in common if x.startswith('Ref_')), any(x for x in uncommon if x.startswith('Ref_'))
			ifref = 10
			if ifref1 and not ifref2 :
				ifref = 1
				
			elif ifref2 and not ifref1:
				ifref = -1
			else:
				ifref = 0
				
				
			allbreaks_found.append([coordi,totalbreak, totalcontigs, ifref, uncommon , common]) 
			
			
			
	locus_addbreaks = cl.defaultdict(set)
	locus_removebreaks = cl.defaultdict(set)
	
	refbreaks = dict()
	novelbreaks = set()
	
	for abreak in allbreaks_found:
		
		coordi, totalbreak, totalcontigs, ifref, uncommon, common = abreak
		
		#mindis = min([abs(x-coordi) for x in existingbreak]+[10000000])
		ifuse = 0 
		
		if ifref == 1:
			refbreaks[coordi] = 1
			
		elif ifref == -1:
			refbreaks[coordi] = -1
			
		elif ifref == 0:
			refbreaks[coordi] = 0
			
			
		elif totalbreak > 0.5*totalcontigs:
			novelbreaks.add((totalbreak, coordi))
			
	return refbreaks, novelbreaks, breaks_allspancontigs


class graphDB:
	
	def __init__(self):
		self.breaks = cl.defaultdict(list)
		self.blocks = cl.defaultdict(list)
		self.blocksize = cl.defaultdict(list)
		self.blockset = cl.defaultdict(set)
		self.allgenes = cl.defaultdict(list)
		self.break_uniformed = cl.defaultdict(list)
		self.chunkindex_togene = cl.defaultdict(set)
		self.genechunks_index = cl.defaultdict(list)  
		self.validregions = cl.defaultdict(list)      
		self.largechunks = set()
		self.smallchunks = set()
		
	def save_json(self, path):
		data = {
		"breaks": self.breaks,
		"blockset": {
			f"{k[0]}|{k[1]}|{k[2]}": v for k, v in self.blockset.items()
		},
		"blocksize": self.blocksize,
		"chunkindex_togene": {
			str(k): [[list(t[0]), t[1]] for t in v]
			for k, v in self.chunkindex_togene.items()
		},
		"genechunks_index": {
			",".join(map(str, k)): v for k, v in self.genechunks_index.items()
		},
		"validregions": self.validregions,
		}
		with open(path, "w") as f:
			json.dump(data, f)
			
			
	@classmethod
	def load_json(cls, path):
		with open(path) as f:
			data = json.load(f)
			
		obj = cls()
		obj.breaks = cl.defaultdict(list, data["breaks"])
		
		obj.blockset = cl.defaultdict(set, {
		(k.split("|")[0], int(k.split("|")[1]), int(k.split("|")[2])): v
		for k, v in data["blockset"].items()
		})
		
		obj.blocksize = cl.defaultdict(list, {
		int(k): v for k, v in data["blocksize"].items()
		})
		obj.chunkindex_togene = cl.defaultdict(set, {
		int(k): {(tuple(x[0]), x[1]) for x in v}
		for k, v in data["chunkindex_togene"].items()
		})
		obj.genechunks_index = cl.defaultdict(list, {
		tuple(map(int, k.split(","))): v for k, v in data["genechunks_index"].items()
		})
		obj.validregions = cl.defaultdict(list, data.get("validregions", {}))
		
		return obj
	
	def load_breaks(self, break_uniformed):
		
		self.break_uniformed = break_uniformed
		
	def load_refbreaks(self, refbreaks):
		
		refbreaks = [(x[0],self.break_uniformed.get( (x[0],x[1]), x[1] ) ) for x in refbreaks]
		
		for abreak in refbreaks:
			
			self.breaks[abreak[0]].append(abreak[1])
			
		for path, breaks in self.breaks.items():
			
			self.breaks[path] = sorted(list(set(breaks)))
			
	def load_refchunks(self):
		
		index_offsite = 0
		for path, breaks in self.breaks.items():
			
			size = int(path.split('_')[-1]) - int(path.split('_')[-2])
			
			breaks = [0]+breaks+[size+1]
			blocks = [(i+index_offsite, path, breaks[i], breaks[i+1]) for i in range(len(breaks)-1) if breaks[i+1] - breaks[i] > 100]
			
			self.blocks[path].extend(blocks)
			
			for i, block in enumerate(blocks):
				self.blockset[tuple([block[1],block[2],block[3]])] = i+index_offsite
				self.blocksize[i+index_offsite] = block[3] -  block[2]
			index_offsite += len(blocks)
			
		return self
	
	def overlap_chunks(self, chunks):
		
		chunks = sorted(chunks, key = lambda x: x[3] + x[4])
		chunks_new = []
		for chunk in chunks:
			
			
			block_index = self.blockset.get(tuple(chunk[:3]),-1)
			
			if block_index == -1:
				
				overlaps = []
				for block,index in self.blockset.items():
					
					if block[0] != chunk[0]:
						continue
					
					overlap = abs(block[2] - block[1]) + abs(chunk[2]-chunk[1]) - max(block[2], block[1], chunk[2], chunk[1]) + min(block[2], block[1], chunk[2], chunk[1])
					
					if overlap > 0.5 * max(abs(block[2] - block[1]),abs(chunk[2]-chunk[1])):
						
						overlaps.append((overlap, index))
						
				overlaps = sorted(overlaps, reverse = 1)
				
				if len(overlaps) and overlaps[0][0] > 100:
					block_index = overlaps[0][1]
					
			chunks_new.append([block_index] + chunk)
			
		return chunks_new
	
	def load_refgenes(self, loci):
		
		def cutspan_bybreaks(spans_flat, breaks):
			
			path_spans = cl.defaultdict(list)
			for (path, gposi, qposi, strd) in spans_flat:
				
				path_spans[path].append((self.break_uniformed.get((path,gposi), gposi),qposi,strd))
				
				
			allchunks = []
			for path, spans_onpath in path_spans.items():
				
				breaks_onpath = breaks[path]
				
				spans_onpath = [sorted([spans_onpath[2*i],spans_onpath[2*i+1]], key = lambda x: x[0]) for i in range(len(spans_onpath)//2)]
				spans_onpath = [(x[0],x[1],y[0][2]) for y in spans_onpath for x in y]
				
				allgcoordis = [x[0] for x in spans_onpath]
				allqcoordis = [(x[1],1 if x[2] == '>' else -1) for x in spans_onpath]
				
				l1 = len(allgcoordis)
				
				allgcoordis =  allgcoordis  + breaks_onpath
				allgcoordis_sortindex = sorted(list(range(len(allgcoordis))), key = lambda x: allgcoordis[x])
				
				current_segs = []
				
				lastqposi = cl.defaultdict(list)
				lastbreak = cl.defaultdict(int)
				for x in allgcoordis_sortindex:
					
					coordi = allgcoordis[x]
					if x < l1:
						if x %2 == 0:
							current_segs.append(x//2)
							lastbreak[x//2] = coordi
							lastqposi[x//2] = allqcoordis[x]
						else:
							allchunks.append([path, lastbreak[x//2], coordi, lastqposi[x//2][0], allqcoordis[x][0]])
							
							current_segs.remove(x//2)
					else:
						for breakindex in current_segs:
							qposi = lastqposi[breakindex][0] + lastqposi[breakindex][1]  * (coordi - lastbreak[breakindex] )
							
							allchunks.append([path, lastbreak[breakindex], coordi, lastqposi[breakindex][0], qposi])
							lastbreak[breakindex] = coordi
							lastqposi[breakindex] = [qposi, lastqposi[breakindex][1]]
							
							
			chunks = sorted([x for x in allchunks if x[2] - x[1] > 100], key = lambda x: x[3])
			
			return chunks
		
		self.genechunks = cl.defaultdict(list)
		self.genechunks_index = cl.defaultdict(list)
		self.chunkindex_togene = cl.defaultdict(set)
		for genes, spans, cigarstr in loci:
			
			for gene in genes:
				
				genename = gene[0][-1]
				
				gene = [tuple([x[0],x[1],x[2],x[3]]) for x in gene]
				
				chunks = cutspan_bybreaks(gene, self.breaks)
				
				chunks = [x for x in chunks if abs(x[2] - x[1]) > 100]
				
				chunks_new = self.overlap_chunks( chunks)
				
				self.genechunks[genename] = chunks_new
				
				self.genechunks_index[tuple([x[0] for x in chunks_new if x[0] >=0 ])].append(genename)
				self.genechunks_index[tuple([x[0] for x in chunks_new[::-1] if x[0] >=0 ])].append(genename)
				
				for i,x in enumerate(chunks_new):
					self.chunkindex_togene[x[0]].add((tuple([x[0] for x in chunks_new if x[0] >=0 ]), i))
					self.chunkindex_togene[x[0]].add((tuple([x[0] for x in chunks_new[::-1] if x[0] >=0 ]), len(chunks_new)-1-i))
					
		self.genes = list(self.genechunks_index.keys())
		self.genenum = len(self.genechunks_index)
		
		return self
	
	def getchunks_fromspan(self, spans, cigarstr):
		
		
		allchunks = []
		
		cigars = re.findall(r'[><][^><]+', cigarstr)
		
		
		path_spans = []
		
		for (path, grange, qrange), cigar in zip(spans,cigars):
			
			path_spans.append([grange, qrange, path[1:], path[0], cigar])
			
		pathsizes = dict()
		for grange, qrange, path, strd, cigar in path_spans:
			
			qstart = min(qrange)
			
			breaks_onpath = self.breaks[path]
			
			pathsize = int(path.split('_')[3]) - int(path.split('_')[2])
			
			#allgcoordis = [self.break_uniformed.get((path, x),x) for y in spans_onpath for x in y[0]]
			allbreaks = [x for x in breaks_onpath if x >= grange[0] and x < grange[1]]
			
			allgcoordis = [grange[0]] + allbreaks + [grange[1]]
			
			if strd == "<":
				breaks_qposi = gposi_toqposi(cigar, [pathsize - x for x in allgcoordis][::-1])
				breaks_qposi = breaks_qposi[::-1]
			else:
				breaks_qposi = gposi_toqposi(cigar, allgcoordis)
				
			breaks_qposi = [x + qstart for x in breaks_qposi]
			
			
			lastgposi, lastqposi = allgcoordis[0], breaks_qposi[0]
			for gposi, qposi in zip(allgcoordis[1:], breaks_qposi[1:]):
				
				if abs(qposi - lastqposi)>1000:
					allchunks.append([path,lastgposi, gposi, lastqposi, qposi])
				lastgposi, lastqposi =  gposi,qposi
				
				
		allchunks = sorted([[x[0]]+sorted(x[1:3])+sorted(x[3:]) for x in allchunks if abs(x[2] - x[1]) > 100], key = lambda x: x[3]+x[4])
		
		allchunks_index = self.overlap_chunks(allchunks)
		
		return allchunks_index
	
	def update_database(self, allchunks, usedchunks, results):
		
		index_tochunkindex = dict()
		
		newchunks = []
		for i,chunk in enumerate(allchunks):
			if chunk[0] == -1:
				
				findoverlap = 0
				for pastnewchunk in newchunks:
					overlap =  abs(chunk[3]-chunk[2]) +  abs(pastnewchunk[3]-pastnewchunk[2])  - max(pastnewchunk[3], pastnewchunk[2], chunk[3] , chunk[2]) + min(pastnewchunk[3], pastnewchunk[2], chunk[3] , chunk[2])
					if overlap > 0.5 * max(abs(chunk[3]-chunk[2]),abs(pastnewchunk[3]-pastnewchunk[2])):
						newindex = chunk[0]
						findoverlap = 1
						break
					
				if findoverlap == 0:
					newindex = len(self.blockset)
					self.blocks[chunk[1]].append([chunk[2],chunk[3]])
					self.blockset[(chunk[1],chunk[2],chunk[3])] = newindex
					self.blocksize[newindex] = abs(chunk[3] - chunk[2])
					index_tochunkindex[i] = newindex
					chunk[0] = newindex
					newchunks.append(chunk)
				else:
					index_tochunkindex[i] = newindex
					chunk[0] = newindex
					
		novels = []
		novel = []
		lastqend = -1000000
		for i,x in enumerate(allchunks):
			if i == -1:
				qstart,qend = x[4],x[5]
				if abs(qstart - lastqend) < 30000:
					novel.append([i,qstart,qend])
				else:
					if len(novel):
						novels.append(novel)
					novel = [[i,qstart,qend]]
				lastqend = qend
				
		if len(novel):
			novels.append(novel)
			
		novels = [x for y in novels for x in polishranges(y) if len(x)]
		
		for result in results:
			chunks_new = tuple([index_tochunkindex.get(x,y) for x,y in zip(result[0],result[1])])
			if chunks_new in self.genechunks_index or len(chunks_new) <= 1:
				continue
			
			genename = "New_"+str(len(self.genechunks))
			
			self.genechunks[genename] = chunks_new
			self.genechunks_index[chunks_new].append(genename)
			self.genechunks_index[chunks_new[::-1]].append(genename)
			
			for i,x in enumerate(chunks_new):
				self.chunkindex_togene[x].add((chunks_new, i))
				self.chunkindex_togene[x].add((chunks_new[::-1], len(chunks_new)-1-i))
				
		for chunks_new in novels:
			
			chunks_new = [x[0] for x in chunks_new]
			chunks_new = tuple(range(min(chunks_new),max(chunks_new)+1))
			
			if chunks_new in self.genechunks_index or len(chunks_new) <= 1:
				continue
			
			genename = "Novel_"+str(len(self.genechunks))
			
			self.genechunks[genename] = chunks_new
			self.genechunks_index[chunks_new].append(genename)
			self.genechunks_index[chunks_new[::-1]].append(genename)
			
			for i,x in enumerate(chunks_new):
				self.chunkindex_togene[x].add((chunks_new, i))
				self.chunkindex_togene[x].add((chunks_new[::-1], len(chunks_new)-1-i))
				
		self.genes = set(list(self.genechunks_index.keys()))
		self.genenum = len(self.genechunks_index)
		
		
	def overlap_genes(self, spans, cigarstr, ifupdate = 0):

		def determinegenes(allchunks,passed_genes):
			
			passed_genes = sorted([x for x in passed_genes if len(x[2])], key = lambda x: max( [  y[2] for y in x[2] ] ), reverse =1 )
			
			
			selected_genes = []
			exclude = set()
			for gene in passed_genes:
				
				chunk_usedindex = list(range(min(gene[0]),max(gene[0])+1))
				
				
				if exclude & set(chunk_usedindex):
					continue
				
				selected_genes.append(gene)
				
				if len(gene[0]):
					exclude.update(chunk_usedindex)
					
			for ichunk, chunk in enumerate(allchunks):
				
				chunkindex = chunk[0]
				qstart,qend = chunk[4],chunk[5]
				
				if chunkindex  >= 0 and ichunk not in exclude:
					
					gene = [[ichunk], [chunkindex], [], qstart, qend , [ [qstart, qend]]]
					selected_genes.append(gene)
					
			return selected_genes,exclude
		
		def extendgenes(current_genes, passed_genes, ichunk, chunkindex, qstart, qend, ifupdate = 1):
			
			current_genes = sorted(current_genes, key = lambda x: len(x[1]), reverse =1 )
			
			continued_genes = []
			
			mingap = 1000000000
			ifmatch = 0
			for index,gene in enumerate(current_genes):
				
				insersize = qstart - gene[4]
				
				if insersize > 30000:
					passed_genes.append(gene)
					continue
				
				involves = gene[2]
				
				involves_new = []
				for (involvegene, aligns, score) in involves:
					
					alignindex = aligns[-1]
					
					if chunkindex in involvegene[(alignindex+1):]:
						
						newindex = alignindex+1 +involvegene[(alignindex+1):].index(chunkindex)
						
						delsize = sum([self.blocksize[involvegene[x]] for x in range(alignindex+1, newindex)])
						
						newscore = qend - qstart - 2*insersize - 2*delsize
						
						score += newscore
						
						mingap = min(mingap, insersize+delsize)
						
						if delsize < 30000 and newscore >= 0:
							involves_new.append((involvegene, aligns+[newindex], score))
							
				if len(involves_new):
					
					if not ifupdate:
						gene_old = [x if type(x) != type([]) else [y for y in x] for x in gene]
						continued_genes.append(gene_old)
						
					gene[0].append(ichunk)
					gene[1].append(chunkindex)
					gene[2] = involves_new
					gene[4] = qend
					gene[5].append([qstart, qend])
					continued_genes.append(gene)
					
				elif len(gene[0]) >0:
					
					if ifupdate:
						continued_genes.append(gene)
					else:
						passed_genes.append(gene)
						
						
			return continued_genes, mingap
		
		current_genes = cl.defaultdict(list)
		passed_genes = []
		allchunks = self.getchunks_fromspan(spans, cigarstr)
		allchunks = sorted(allchunks, key = lambda x: x[4] + x[5])
		
		lastchunk = -100000
		for ichunk, chunk in enumerate(allchunks):
			
			chunkindex = chunk[0]
			qstart,qend = chunk[4],chunk[5]
			
			allinvolve_genes = self.chunkindex_togene.get(chunkindex,[])
			
			if len(allinvolve_genes) == 0 and len(current_genes) == 0:
				continue
			
			continued_genes, mingap = extendgenes(current_genes, passed_genes, ichunk, chunkindex, qstart, qend, ifupdate)
			
			current_genes = continued_genes
			
			new_gene = [[ichunk], [chunkindex], [], qstart, qend , [ [qstart, qend]]]
			
			allinvolve_genes = self.chunkindex_togene.get(chunkindex,[])
			for (gene,alignindex) in allinvolve_genes:
				
				new_gene[2].append((gene,[alignindex],qend - qstart))
				
			current_genes.append(new_gene)
			
		passed_genes.extend(current_genes)
		
		results, used_chunks = determinegenes(allchunks,passed_genes)
		
		for result in results:
			result[2] = [(x[2],self.genechunks_index.get(x[0], [""])[0]) for x in result[2]]
			if len(result[2]):
				result[2] = result[2][0]
				
		results = sorted(results, key = lambda x: tuple(x[0]))
		
		return results 
		
		
	
def getblocks(gbreaks, break_uniformed, contigspans, genetoscaff, graphinfo, blacklists):
	
	saveorload = 0
	
	allloci = [x for x in gbreaks.keys() if x.startswith("Ref_")] + [x for x in gbreaks.keys() if not x.startswith("Ref_")]
	if saveorload == 0:
	
		thegraph = graphDB()
		thegraph.load_breaks(break_uniformed)
	
		refgenes = []
		for name, gbreak in gbreaks.items():
			
			if name.startswith("Ref_"):
					
				thegraph.load_refbreaks([(x[0],x[1]) for y in gbreak for x in y])
					
				refgenes.append([gbreaks[name], contigspans[name], graphinfo[name]])
					
		thegraph.load_refchunks().load_refgenes(refgenes)
	
		results = dict()
		for loci in allloci:
			x = 1 if loci not in blacklists else 0
			results[loci] = thegraph.overlap_genes(contigspans[loci], graphinfo[loci], ifupdate=x)
			
			
		for loci in allloci:
			results[loci] = thegraph.overlap_genes(contigspans[loci], graphinfo[loci])
	
	"""
	if saveorload == 1:
		thegraph = graphDB()
		thegraph = thegraph.load_json(cache)
	
		results = dict()
		for loci in allloci:
			results[loci] = thegraph.overlap_genes(contigspans[loci], graphinfo[loci])
	"""
	
	return results,thegraph


def breaks_ongraph(genetoscaff, genetolocus, graphinfo, contigspan, blacklists):
	
	gbreaks = cl.defaultdict(list)

	segmentinfo =  cl.defaultdict(list)
	for name, infos in genetolocus.items():
		
		for info in infos:
			
			locus,start,end = info
		
			if locus.startswith("Ref_"):
				graphinfo[locus] = graphinfo[locus[4:]]
			
				if locus[4:] in contigspan:
					contigspan[locus] = contigspan[locus[4:]]
					del contigspan[locus[4:]]
						
			if locus not in graphinfo:
				continue
		
			coordinates = qposi_togposi(graphinfo[locus], [start,end-1])
		
			gbreaks[locus].append([list(x)+[name] for x in coordinates])
		
			segmentinfo[locus].append([coordinates[0],coordinates[1], name])
				
	gbreaks = {locus:sorted(regions, key = lambda x: x[0][2]) for locus, regions in gbreaks.items() }

	breaks_andspans = cl.defaultdict(list)
	for name, data in gbreaks.items():
		breaks_andspans[name] = [x for x in data]
	for name, data in contigspan.items():
		breaks_andspans[name].append([[x[0][1:],x[1][y],x[2][y],x[0][0],""] for x in data for y in [0, 1]])
		
	break_uniformed = combinebreaks(breaks_andspans)

	results,thegraph = getblocks(gbreaks, break_uniformed, contigspan, genetoscaff, graphinfo, blacklists)

	return results,thegraph

def assemblysmall(allresults, cutoff = 20000):
	
	def ifalign(thelist,x):
		
		l = len(x)
		for i in range(len(thelist)):
			if x == thelist[i:(i+l)]:
				return 1
		return 0
	
	smallchunks = set()
	bindrecords = cl.defaultdict(int)
	iffix = set()
	for loci,results in allresults.items():
		
		if len(results) == 1:
			continue
		
		for i,result in enumerate(results):
			
			size = result[4] - result[3]
			
			if size < cutoff:
				iffix.add(loci)
				smallchunks.add(tuple(sorted(result[1])))
				
				if i == 0:
					
					bindindex2 = tuple(sorted([result[1][-1], results[i+1][1][0]]))
					bindrecords[bindindex2] += 1
					continue
				elif i == len(results)-1:
					bindindex1 = tuple(sorted([result[1][0], results[i-1][1][-1]]))
					
					bindrecords[bindindex1] += 1
					continue
				else:
					bindindex1 = tuple(sorted([result[1][0], results[i-1][1][-1]]))
					bindindex2 = tuple(sorted([result[1][-1], results[i+1][1][0]]))
					if abs(results[i+1][4] - results[i+1][3]) >= abs(results[i-1][4] - results[i-1][3]):
						bindrecords[bindindex2] += 1
					else:
						bindrecords[bindindex1] += 1
						
	for loci,results in allresults.items() :
		if loci not in iffix:
			continue
		connects = [0 for x in results]
		
		for i,result in enumerate(results):
			
			if tuple(sorted(result[1])) in smallchunks:
				
				if i == 0:
					connects[i] = 1
					continue
				elif i == len(results)-1:
					connects[i-1] = 1
					continue
				else:
					bindindex1 = tuple(sorted([result[1][0], results[i-1][1][-1]]))
					bindindex2 = tuple(sorted([result[1][-1], results[i+1][1][0]]))
					if bindrecords[bindindex2] > bindrecords[bindindex1]:
						connects[i] = 1
					else:
						connects[i-1] = 1
						
		lastindex = 0
		newresults = [results[0]]
		for i,result in enumerate(results[1:]):
			
			if connects[i] == 1:
				
				newresults[-1][0].extend(result[0])
				newresults[-1][1].extend(result[1])
				newresults[-1][3] = min(newresults[-1][3], result[3])
				newresults[-1][4] = max(newresults[-1][4], result[4])
				newresults[-1][5].extend(result[5])
				
			else:
				newresults.append(result)
				
		results = newresults
		
		allresults[loci] = results
		
	for loci,results in allresults.items() :
		
		laststart, lastend  = results[0][3], results[0][4]
		for i,result in enumerate(results[1:]):
			
			start, end = result[3], result[4]
			
			distance = max(laststart, lastend,start, end ) - min(laststart, lastend,start, end ) - abs(laststart-  lastend) - abs(start - end)
			if distance > 1 and distance < 1000:
				if laststart < start:
					allresults[loci][i][4] = start
				else:
					allresults[loci][i-1][4] = laststart
					
			laststart, lastend = start, end
			
	return allresults, iffix

def clean_nonoverlap(breaksoncontigs, locuslocations, geneoncontigs):
	
	genelabs = cl.defaultdict(lambda: [0,[]])
	for loci, breaks in breaksoncontigs.items():
		
		scaflocus = locuslocations[loci]
		
		lastend = 0
		for segments in breaks:
			
			for i,segment in enumerate(segments[5]):
				
				start,end = min(segment[0],segment[1]),max(segment[0],segment[1])
				
				if segment[0] < max(start,lastend):
					segment[0] = max(start,lastend)
				if segment[1] < max(end,lastend):
					segment[0] = max(end,lastend)
					
				lastend = max(end,lastend)
				
		oldsegments = geneoncontigs[scaflocus[0]]
		
		for segment in breaks:
			
			subsegments_indexs, subsegments_ranges = segment[1],  segment[5]
			
			for segmentindex, subsegments_range in zip(subsegments_indexs, subsegments_ranges):
				
				start, end = subsegments_range[0], subsegments_range[1]
				
				if scaflocus[1] == '+':
					scafstart = scaflocus[2] + start
					scafend = scaflocus[2] + end
				else:
					scafstart = scaflocus[3] - end
					scafend = scaflocus[3] - start
					
				overlap = sorted([( (x[1] - x[0]) + (scafend - scafstart) - max(x[1],scafend,x[0],scafstart) + min(x[1],scafend,x[0],scafstart), x) for i,x in enumerate(oldsegments)], reverse=1)[0]   
				
				if overlap[0] > genelabs[segmentindex][0]:
					genelabs[segmentindex] = overlap
					
					
	regions = []
	for loci, breaks in breaksoncontigs.items():
		
		scaflocus = locuslocations[loci]
		
		newsegments = []
		for i,segment in enumerate(breaks):
			
			subsegments_indexs, subsegments_ranges = segment[1],  segment[5]
			
			start = 100000000000000
			end = -1
			subsegments_indexs_new = []
			subsegments_ranges_new = []
			for segmentindex, subsegments_range in zip(subsegments_indexs, subsegments_ranges):
				
				size = subsegments_range[1] - subsegments_range[0]
				overlap = genelabs[segmentindex] 
				
				if len(overlap) and overlap[0] > min(1000, 0.5*size):   
					start = min(start, subsegments_range[0])
					end = max(end, subsegments_range[1])
					
					subsegments_indexs_new.append(segmentindex)
					subsegments_ranges_new.append(subsegments_range)
					
			if end == -1:
				continue
			
			newsegment = segment
			
			newsegment[1] =  subsegments_indexs_new
			newsegment[5] =  subsegments_ranges_new
			
			newsegment[3] = start
			newsegment[4] = end
			
			newsegments.append(newsegment)
			
		breaksoncontigs[loci] = newsegments
		
	return breaksoncontigs, oldsegments, genelabs


def annotate_regions(breaksoncontigs, locuslocations, genetoscaff):
	
	geneoncontigs = cl.defaultdict(list)
	for genename, info in genetoscaff.items():
		
		if len(info) > 4:
			geneoncontigs[info[0]].append(info[2:]+[genename])
			
			groupprefix = genename.split('_')[0]
			
	novel_overlap = [x[2][1]  for loci, breaks in breaksoncontigs.items() for x in breaks if len(x[2]) ]
	novel_overlap = list(set(novel_overlap))
	
	breaksoncontigs, oldsegments, genelabs = clean_nonoverlap(breaksoncontigs, locuslocations, geneoncontigs)
	
	breaksoncontigs = {k:v for k,v in breaksoncontigs.items() if len(v)}
	
	
	iffix = 1
	while iffix:
		breaksoncontigs, iffix = assemblysmall(breaksoncontigs)
		
		
	regions = []
	for loci, breaks in breaksoncontigs.items():
		
		scaflocus = locuslocations[loci]
		
		for segment in breaks:
			
			start, end = segment[3], segment[4]
			
			if scaflocus[1] == '+':
				scafstart = scaflocus[2] + start
				scafend = scaflocus[2] + end
			else:
				scafstart = scaflocus[3] - end
				scafend = scaflocus[3] - start
				
			overlap = sorted([( (x[1] - x[0]) + (scafend - scafstart) - max(x[1],scafend,x[0],scafstart) + min(x[1],scafend,x[0],scafstart), x) for i,x in enumerate(oldsegments)], reverse=1)[0]   
			
			regions.append([scaflocus[0], int(scafstart), int(scafend)] + overlap[1][2:]+[loci,start,end])
			
	new_regions = cl.defaultdict(list)
	haplo_counter = cl.defaultdict(int)
	orginal_loci = dict()
	for region in regions:
		
		contig, start, end, ifexon, mapinfo,oldname,locus,localstart,localend = region
		
		haplo = "_h".join(contig.split("#")[:2]) if "#" in contig else "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
		
		haplo_counter[haplo] += 1
		
		newname = groupprefix+"_"+ haplo +"_"+ str(haplo_counter[haplo])
		
		#header = "{}\t{}:{}-{}\t{}\t{}".format(newname, contig, start, end, ifexon, mapinfo  )
		
		new_regions[contig].append([newname,  contig, start, end, ifexon, mapinfo])
		
		orginal_loci[newname] = (locus,localstart,localend)
		
	new_regions = [x[:1]+["_".join(x[1].split("_")[:-1])]+x[2:] for y in new_regions.values() for x in y]
	
	return new_regions, orginal_loci


def coordinate_uniform(cbreaks_sort):
	
	uniform_coordis = []

	currbreaks =  []
	lastbreak = -100
	for coordi in cbreaks_sort:
		if (coordi - lastbreak) < 100:
			currbreaks.append(coordi)
		else:
			if len(currbreaks):
				uniform = cl.Counter(currbreaks).most_common(1)[0][0]
				uniform_coordis.extend([uniform] * 1)
					
			currbreaks = [coordi]
			
		lastbreak = coordi
		
	if len(currbreaks):
		uniform = cl.Counter(currbreaks).most_common(1)[0][0]  
		uniform_coordis.extend([uniform] * 1)
		
		
	return uniform_coordis



def locatebreaksoncontigs(breaksonscaf, locuslocations, genetoscaff):
	
	oldbreaksoncontig = cl.defaultdict(list)
	for gene,obreaks in genetoscaff.items():
		
		scaf, strd, start, end = obreaks[:4]
	
		oldbreaksoncontig[scaf].extend([start, end])
		
	results = cl.defaultdict(list)
	for scafford, cbreaks in breaksonscaf.items():
		
		if len(cbreaks) == 0:
			continue
	
		cbreaks_sort_index = sorted(range(len(cbreaks)), key = lambda x: cbreaks[x])
	
		cbreaks_sort = sorted([cbreaks[x] for x in cbreaks_sort_index])
	
		cbreaks_sort_uniform = coordinate_uniform(cbreaks_sort)
	
		oldmin = min(oldbreaksoncontig[scafford]+[-100]) 
		oldmax = max(oldbreaksoncontig[scafford]+[100000000000])
	
		cbreaks_sort_uniform_rangeindex = [i for i,x in enumerate(cbreaks_sort_uniform) if x > oldmin + 100 and x< oldmax-100 ]
	
		if len(cbreaks_sort_uniform_rangeindex) == 0:
			
			if min(cbreaks_sort_uniform_rangeindex) > 0:
				cbreaks_sort_uniform_rangeindex = [min(cbreaks_sort_uniform_rangeindex)-1] + cbreaks_sort_uniform_rangeindex
			if max(cbreaks_sort_uniform_rangeindex) < len(cbreaks_sort_uniform_rangeindex) - 1:
				cbreaks_sort_uniform_rangeindex =  cbreaks_sort_uniform_rangeindex + [max(cbreaks_sort_uniform_rangeindex)+1]
					
		cbreaks_sort_uniform = [cbreaks_sort_uniform[x] for x in cbreaks_sort_uniform_rangeindex]
	
		if len(cbreaks_sort_uniform ) == 0:
			cbreaks_sort_uniform = oldbreaksoncontig[scafford]
		elif  len(cbreaks_sort_uniform ) <= 1 or oldmin < min(cbreaks_sort_uniform) - 1000 or oldmax > max(cbreaks_sort_uniform) + 1000:
			cbreaks_sort_uniform = [x for x in oldbreaksoncontig[scafford] if x < min(cbreaks_sort_uniform) -1000 ] + cbreaks_sort_uniform + [x for x in oldbreaksoncontig[scafford] if x > max(cbreaks_sort_uniform) + 1000 ]
			
			#cbreaks_sort_uniform = sorted(oldbreaksoncontig[scafford])
			
		results[scafford] = cbreaks_sort_uniform 
			
			
	return results

def MergeChunks(regions, mergedistance):
	
	regions_sort = sorted([x for x in regions], key = lambda x: (x[1], x[2]))

	lastcontig = ""
	lastend = -mergedistance - 1

	merged_regions = []
	current_region = []
	for region in regions_sort:
		
		newname,  contig, start, end, ifexon, mapinfo = region
		if contig != lastcontig:
			lastend = -mergedistance - 1
			
		if start - lastend < mergedistance:
			current_region[3] = end
			current_region[4] = "Exon" if ifexon == "Exon" else current_region[4]
			
		else:
			merged_regions.append(current_region)
			current_region = [x for x in region]
			
		lastend = end
		lastcontig = contig
		
		
	merged_regions = merged_regions[1:]
	merged_regions.append(current_region)


	if len([x for x in merged_regions if len(x)]) == 0:
		return 0, []

	maxsize = max([x[3]-x[2] for x in merged_regions])
	if maxsize < mergedistance:
		
		regions = merged_regions
		
	return maxsize, regions

def uniformbreaks(inputfile, haplomergefile, graphfile, refcontig, blacklists):
	
	genetoscaff, genetolocus, locuslocations = readcontigsinfo(inputfile, haplomergefile, refcontig)

	graphinfo, contigspan = readgraph(graphfile)

	results, thegraph = breaks_ongraph(genetoscaff, genetolocus, graphinfo, contigspan, blacklists)

	regions, orginal_loci = annotate_regions(results,  locuslocations, genetoscaff)
	
	#normalized = json.dumps(sorted(regions, key=lambda x: tuple(map(str, x))), sort_keys=True)
	#hashvalue = hashlib.sha256(normalized.encode()).hexdigest()
	
	if cache:
		allgraphspans = cl.defaultdict(list)
		for region in regions:

			name,  contig, start, end = region[:4]
			
			locus,localstart,localend = orginal_loci[name]
			
			graphcoordi = qposi_togposi(graphinfo[locus], [localstart,localend-1])
			
			for index in range(len(graphcoordi)//2):
				
				path, start = graphcoordi[2*index][:2]
				path, end = graphcoordi[2*index+1][:2]
				
				allgraphspans[path].append(sorted([start, end]))
				
		merged_span = dict()
		for chr, spans in allgraphspans.items():
			
			merged_span[chr] = [x for y in Overlap(spans) for x in y[:2]]
			
		thegraph.validregions = merged_span
		
		thegraph.save_json(cache)
		

	return regions

def main(args):
	
	if len(args.folder)==0: 
		folder = args.output + "_temp"+ datetime.datetime.now().strftime("%y%m%d_%H%M%S/").replace("%","_")
		os.system("rm -rf {} || true".format(folder))
		os.makedirs(folder,exist_ok=True)
		
	else:
		folder = args.folder
		os.system("mkdir {} || true".format(folder))
		
	blacklists = findbreaks(args.input, args.query, args.output+"_loci.txt")
	haplomergefile = args.output+"_loci.txt.fasta"
	makenewfasta(args.output+"_loci.txt", args.query, haplomergefile, folder ,args.threads)
	os.system(f"{script_folder}/KmerStrd -i {haplomergefile} -o {haplomergefile}_ && mv {haplomergefile}_ {haplomergefile}")
	makegraph(haplomergefile, folder, args.threads)

	if args.cache:
		global cache
		cache = haplomergefile + "_graphcache.json"
		
	regions = uniformbreaks(args.input, args.output+"_loci.txt", haplomergefile+"_allgraphalign.out", args.ref, blacklists)

	maxsize, regions = MergeChunks(regions, args.mergeall)

	if len(regions) == 0:
		print(f"gfixbreak.py Error, Missing regions: {args.input}")
		
		
	makenewfasta2(regions, args.query, args.output, folder, args.threads)
	if len(args.folder)==0:
		os.system("rm  -rf {} || true".format(folder))
			
	#os.system(f" cat {haplomergefile} | grep \"^>\" > {haplomergefile}_   &&  mv {haplomergefile}_ {haplomergefile} ")
	os.system(f"rm {haplomergefile}_norm.gz || true ")
	return

def run():
	"""
			Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-n", "--norm", help="path to input data file",dest="norm", type=str, required="")
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	parser.add_argument("-t", "--threads", help="path to output file", dest="threads",type=int, default = 1)
	parser.add_argument("-q", "--query", help="path to output file", dest="query",type=str, default = True)
	parser.add_argument("-d", "--folder", help="path to output file", dest="folder",type=str, default = "")   
	parser.add_argument("-m", "--mergeall", help="path to output file", dest="mergeall",type=int, default = 80000)
	parser.add_argument("-r", "--ref", help="prefix of ref contig", dest="ref",type=str, default = "NC_0609")
	parser.add_argument("--cache", help="save caches", action="store_true")

	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()

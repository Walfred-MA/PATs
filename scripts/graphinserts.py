#!/usr/bin/env python3

import os
import argparse
import re
import collections as cl

pathminsize = 300

class kmeraligns:
	
	def __init__(self, cigars, ksize = 30, cutoff = 50):
		
		self.cigars = re.findall(r'\d+[a-zA-Z=]', cigars)
		self.ksize = ksize
		self.cutoff = cutoff
		self.min = -cutoff
		self.max = cutoff*3
		self.consensus_f = []
		self.consensus_r = []
	
	def getalignsegs(self):
		
		lastzero = 0
		self.consensus_f = [[0,1]]
		rposi = 0
		qposi = 0 
		lastvar = 0
		lastscore = 0
		score = 0 
		currseg = 1
		for index, cigar in enumerate(self.cigars):
			
			lastscore = score
			thesize = int(cigar[:-1])
			thetype = cigar[-1]
			if thetype == "=" or thetype == "M":
				
				qposi += thesize
				rposi += thesize
					
				score += thesize - self.ksize
				score = min(self.max,max(self.min, score))
				
			elif thetype == "X":
				
				qposi += thesize
				rposi += thesize
					
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			elif thetype == "H":
				
				rposi += thesize
				continue
			
			elif thetype == "I" :
				
				rposi += thesize
				
				score -= ( min(self.ksize, qposi - lastvar)  )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			elif thetype == "D":
				
				qposi += thesize
					
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			if lastscore * score < 0 or score == 0:
				
				lastzero = index
				
			if currseg <= 0 and score >= self.cutoff:
				
				self.consensus_f.append([lastzero,lastzero])
				currseg = 1
				
			if currseg > 0 and score <= 0:
				
				
				self.consensus_f[-1][1] = index 
				currseg = -1
				
		if currseg > 0:
			if score >= lastscore:
				self.consensus_f[-1][1] = index + 1
			else:
				self.consensus_f[-1][1] = index 
		
		self.consensus_f = [x for x in self.consensus_f if x[1] > x[0]]
		
		lastzero = len(self.cigars)
		self.consensus_r = [[lastzero-1,lastzero]]
		rposi = 0
		qposi = 0 
		lastvar = 0
		lastscore = 0
		score = 0 
		currseg = 1
		for index, cigar in enumerate(self.cigars[::-1]):
			
			index = len(self.cigars) - index
			lastscore = score
			thesize = int(cigar[:-1])
			thetype = cigar[-1]
			
			if thetype == "=" or thetype == "M":
				
				qposi += thesize
				rposi += thesize
					
				score += thesize - self.ksize
				score = min(self.max,max(self.min, score))
				
			elif thetype == "X":
				
				qposi += thesize
				rposi += thesize
				
				if score > 0 and thesize > self.cutoff:
					score = 0
					
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			elif thetype == "H":
				
				rposi += thesize
				continue
			
			elif thetype == "I" :
				
				rposi += thesize
				
				score -= ( min(self.ksize, qposi - lastvar)  )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			elif thetype == "D":
				
				qposi += thesize
				
				if score > 0 and thesize > self.cutoff:
					score = 0
					
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			if lastscore * score < 0 or score == 0:
				
				lastzero = index
				
			if currseg <= 0 and score >= self.cutoff:
				
				self.consensus_r.append([lastzero,lastzero])
				currseg = 1
				
			if currseg > 0 and score <= 0:
				
				self.consensus_r[-1][0] = index 
				currseg = -1
		
		if currseg > 0:
			if score >= lastscore:
				self.consensus_r[-1][0] = 0
			else:
				self.consensus_r[-1][0] = 1
		
		self.consensus_r = [x for x in self.consensus_r if x[1] > x[0]]
		
		
		return self
		
	def combineseg(self):
		
		consensus = self.consensus_f + self.consensus_r
		consensus_pos = [a for b in consensus for a in b]
		
		sortindexes = sorted(range(len(consensus_pos)), key = lambda x: consensus_pos[x])
		
		num_curr = 0
		self.consensus_comb = []
		for index, sortindex in enumerate(sortindexes):
			
			if sortindex % 2 == 0:
				num_curr += 1
				
				if num_curr == 1:
					
					gapsize = 10000000
					if len(self.consensus_comb):
						lastend = self.consensus_comb[-1][-1]
						gapsize = sum([int(self.cigars[x][:-1]) for x in range(lastend,sortindex) if self.cigars[x][-1]])
					
					if gapsize > self.cutoff:
						self.consensus_comb.append([consensus_pos[sortindex], consensus_pos[sortindex]])
					
			else:
				num_curr -= 1
				
				if num_curr == 0:
					
					self.consensus_comb[-1][-1] = consensus_pos[sortindex]
		
		return self
	
	def trim(self):
		
		self.consensus_trimed = []
		for segment in self.consensus_comb:
			
			lastvar = 0
			score = -1
			lastzero = -1
			rposi = 0
			qposi = 0 
			lastvar = 0
			lindex = segment[1]-segment[0]
			rindex = segment[1]-segment[0]
			lpass = 0
			for index in range(segment[0],segment[1]):
				
				lastscore = score
				cigar = self.cigars[index] 
				thesize = int(cigar[:-1])
				thetype = cigar[-1]
				
				if thetype == "=" or thetype == "M":
					
					if score < 0 and thesize > self.cutoff:
						score = 0
						
					score += thesize - self.ksize
					score = min(self.max,max(self.min, score))

					if lastscore < 0 and score >= 0:
						
						lastzero = index					

					if score > self.ksize:
						lindex = lastzero if lastscore >= 0 else index
						lpass = 1
						break
					
				elif thetype == "X":
					
					if score > 0 and thesize > self.cutoff:
						score = 0
						
					score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
					
				elif thetype == "H":
					
					rposi += thesize
					
				
				elif thetype == "I" :
					
					rposi += thesize
					
					score -= ( min(self.ksize, qposi - lastvar)  )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
					
				elif thetype == "D":
					
					qposi += thesize
						
					score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
			
			score = -1
			lastzero = -1
			rposi = 0
			qposi = 0 
			lastvar = 0
			rpass = 0
			for index in list(range(segment[0],segment[1]))[::-1]:
				
				lastscore = score
				cigar = self.cigars[index] 
				thesize = int(cigar[:-1])
				thetype = cigar[-1]
				
				if thetype == "=" or thetype == "M":
					
					if score < 0 and thesize > self.cutoff:
						score = 0
						
					score += thesize - self.ksize
					score = min(self.max,max(self.min, score))
					
					if lastscore < 0 and score >= 0:
						
						lastzero = index

					if score >= self.ksize:
						rindex = lastzero if lastscore >= 0 else index
						rpass = 1
						break
					
				elif thetype == "X":
					
					score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
					
				elif thetype == "H":
					
					rposi += thesize
					
				
				elif thetype == "I" :
					
					rposi += thesize
					
					score -= ( min(self.ksize, qposi - lastvar)  )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
					
				elif thetype == "D":
					
					qposi += thesize
						
					score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
			
			if rindex >= lindex and lpass and rpass:
				self.consensus_trimed.append([lindex, rindex+1])
		
		return self
	def coordinate(self):
		
		self.segment_coordinates = [[] for x in self.consensus_comb]
		
		allstartends = {}
		for index,segment in enumerate(self.consensus_trimed):
			
			allstartends[2*segment[0]] = index
			allstartends[2*segment[1]-1] = index
			
		peakscore = 0
		score = 0
		rposi = 0
		qposi = 0 
		lastvar = 0
		for index, cigar in enumerate(self.cigars):
			
			if 2*index in allstartends:
				self.segment_coordinates[allstartends[2*index]].append(qposi)
				score = 0
				peakscore = 0
			
			thesize = int(cigar[:-1])
			thetype = cigar[-1]
			
			if thetype == "=" or thetype == "M":
				
				qposi += thesize
				rposi += thesize
				score += thesize - self.ksize
				score = max(self.min, score)
				peakscore = max(peakscore, score)
				
			elif thetype == "X":
				
				qposi += thesize
				rposi += thesize
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = max(self.min, score)
				
				lastvar = qposi
				
			elif thetype == "H":
				
				rposi += thesize
			
			elif thetype == "I" :
				
				rposi += thesize
				
				score -= ( min(self.ksize, qposi - lastvar)  )
				score = max(self.min, score)
				
				lastvar = qposi
				
			elif thetype == "D":
				
				qposi += thesize
				
				score -= ( min(self.ksize, qposi - lastvar)  )
				score = max(self.min, score)
				
				lastvar = qposi
				
			if 2*index + 1 in allstartends:
				
				
				segindex = allstartends[2*index + 1]
				self.segment_coordinates[segindex].append(qposi)
				self.segment_coordinates[segindex].append(peakscore)
				cigar_segs = "".join([self.cigars[x] for x in range(*self.consensus_trimed[segindex])])
				
				self.segment_coordinates[segindex].append(cigar_segs)
				
		
		return self
	

def findkmeraligns(cigars, ksize = 30, cutoff = 100):

	consensus = kmeraligns(cigars, ksize, cutoff).getalignsegs().combineseg().trim().coordinate()

	output = [x for x in consensus.segment_coordinates if len(x) and x[1] - x[0] > cutoff and x[2] > cutoff]
	
	return  output

def makereverse(seq):
	
	tran=str.maketrans('ATCGatcg', 'TAGCtagc')
	
	return seq[::-1].translate(tran)

def pathtoseq(path, pathstrs):
	
	pathes = re.findall(r'[><]s\d+', path)
	
	fullseq = ""
	for eachpath in pathes:
		
		strd = eachpath[0]
		
		seq = pathstrs[eachpath[1:]]
		
		if strd == "<":
			
			seq = makereverse(seq)
			
		fullseq += seq
		
	return fullseq

def unalignregion(allaligns, allow_gap = pathminsize):
	
	all_coordinates = [x for y in allaligns for x in y[:2]]
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	allgroups = []
	current_group_coordi = []
	
	number_curr_seq = 0
	last_coordinate = 0
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		
		if index %2 == 0 :
			
			if number_curr_seq == 0 and coordinate - last_coordinate>allow_gap :
				
				allgroups.append([last_coordinate, coordinate])
				
			number_curr_seq += 1
			
		else:
			
			number_curr_seq -= 1
			
		last_coordinate = coordinate
		
		
	return allgroups


def findaligns(cigars, posi):
	
	cigars = re.findall(r'\d+[a-zA-Z=]', cigars)

	newcigars = []
	
	refposi = 0 
	aligns = []
	newref = ""
	for cigar in cigars:
		
		thesize = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == "=" or thetype == "M":
			
			if thesize > pathminsize:
				aligns.append([posi, posi+thesize])
				
			posi += thesize
			refposi += thesize
			
		elif thetype == "X" or thetype == "S":
			
			posi += thesize
			refposi += thesize
			
		elif thetype == "D" or thetype == "H":
			
			refposi += thesize
			
		elif thetype == "I":
			
			posi += thesize
			
			
	return aligns


def insertdistract(alignfile, queryseq, output, singlefile = 1, length = pathminsize):
	
	with open(queryseq, mode = 'r') as f:
		reads = [read.splitlines() for read in f.read().split(">")[1:]]
		reads = {read[0].split()[0]:"".join(read[1:]) for read in reads}
		
		
	insertfiles = []
	foundseqs = []
	insert_count = 0
	
	
	with open(alignfile, mode = 'r') as f:
		
		for line in f:
			
			if len(line) == 0 or line[0]=="@":
				continue
			
			elements = line.strip().split()
			#name, qsize, qstart, qend, path, rsize, rstart, rend = elements[:8]

			cigar = elements[5]

			identity =float( [x for i,x in enumerate(elements) if x[:5] == "PI:f:"][0][5:])
			match = int( [x for i,x in enumerate(elements) if x[:5] == "AS:i:"][0][5:] ) 
			
			if identity < 90 or match < pathminsize or elements[4] != '255':
				continue
			
			#foundseqs[name].append([int(qstart), int(qend)])
			alignsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XD]', cigar)]+[0])
			#index, cigar = [(i,x) for i,x in enumerate(elements) if x[:5] == "cg:Z:"][0]
			
			#aligns = findkmeraligns(cigar)
			aligns = [[0, alignsize]]
			foundseqs.extend([[x[0]+int(elements[3])-1, x[1]+int(elements[3])-1] for x in aligns])  
			
	for seqname, seq in reads.items():
		
		unaligned = unalignregion([[0,0]]+foundseqs+[[len(seq)-1, len(seq)-1]])
	
		for gap in unaligned:
			
			insertname = output+"_i{}.fasta".format(insert_count)
			insert_count += 1
			insertfiles.append(insertname) 
			
			if singlefile:
				outpath = output
				flag = 'a'
			else:
				outpath = insertname
				flag = 'w'
			
			start = max(0,gap[0]-60)
			end = min( len(seq), gap[1]+60)
			size = end - start
            
			outseq = seq[start:end]
			unmasked =   sum(1 for c in outseq if c.isupper())
			masked = len(outseq) - unmasked
			if len(re.findall(r'[^ATCGatcg]', outseq)) or unmasked + masked < pathminsize:
				continue
	
			with open(outpath, mode = flag) as w:
				
				start = max(0,gap[0]-60)
				end = min( len(seq), gap[1]+60)

				size = end - start
				if size < length:

					if start == 0:

						end += length - size
				
					else:
						
						start -= min(start, length - size)
				outname = ">{}_{}_{}".format(seqname,start, end)
					
				w.write(">{}_{}_{}\n{}".format(seqname,start, end, seq[start:end])+"\n")
				
	return insertfiles


def main(args):
	
	
	insertdistact(args.input, args.query , args.output, args.singlefile, args.length)
	
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	
	parser.add_argument("-q", "--query", help="path to output file", dest="query", type=str,required=True)  
	parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
	parser.add_argument("-s", "--singlefile", help="path to output file", dest="singlefile", type=int, default = 1)
	parser.add_argument("-l", "--length", help="path to output file", dest="length", type=int, default = 300)
	
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()

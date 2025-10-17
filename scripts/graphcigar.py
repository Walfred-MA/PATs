#!/usr/bin/env python3

import os
import argparse
import re
import bisect

pathminsize = 30

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
						lindex = lastzero if lastscore >=0 else index
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
						rindex = lastzero if lastscore >=0 else  index
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
	
	
def findkmeraligns(cigars, ksize = 30, cutoff = 300):
	
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

def graphicalign(allaligns, allow_gap = pathminsize):
	
	all_coordinates = [x for y in allaligns for x in y[:2]]
	scores = [y[2] for y in allaligns]
	lineindex = [y[3] for y in allaligns]
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	current_alignscores = []
	allgroups = [[0,0,0,0]]
	
	ended_indexes = []
	
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		score = scores[index//2]
		
		if index %2 == 0 :
			
			
			if len(current_alignscores) ==0 or score > current_alignscores[-1][0]:
				
				if len(current_alignscores):
					allgroups[-1][-1] = coordinate
				allgroups.append([score, index//2, coordinate, coordinate])
				
			bisect.insort(current_alignscores, [score, -(index//2)])
			
		else:
			
			allgroups[-1][-1] = coordinate
			
			highest_index = -current_alignscores[-1][1]
			
			current_alignscores.remove([scores[index//2], -(index//2)])
			
			if  index//2 == highest_index  and len(current_alignscores):
				
				allgroups.append([current_alignscores[-1][0], -current_alignscores[-1][1]]+[coordinate, coordinate])
				
				
	allgroups = [[lineindex[x[1]]]+x[2:] for x in allgroups[1:]]
	
	last_coordinate = - allow_gap - 1
	last_index = -1
	
	
	new_allgroups = []
	for segment in allgroups:
		
		index,start, end = segment[0], segment[1], segment[2]
		
		if index == last_index  :
			new_allgroups[-1][-1] = end
		else:
			new_allgroups.append(segment)
			
		last_index = index 
		last_end = end
		
	allgroups = [x for x in new_allgroups if x[2] - x[1] > allow_gap]
	
	
	
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
			
			if thesize > 30:
				aligns.append([posi, posi+thesize])
				
			posi += thesize
			refposi += thesize
			
		elif thetype == "X" or thetype == "S":
			
			posi += thesize
			refposi += thesize
			
		elif thetype == "I" or thetype == "H":
			
			refposi += thesize
			
		elif thetype == "D":
			
			posi += thesize
			
	return aligns

def cigar_findqrange(cigars, qstart , qend):
	
	cigars = re.findall(r'\d+[a-zA-Z=]', cigars)
	
	
	refposi = 0
	posi = 0
	newcigars = []
	refstart = -1
	refend = -1
	newref = ""
	for cigar in cigars:
		
		
		lastposi = posi
		
		thesize = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == "=" or thetype == "M":
			
			posi += thesize
			refposi += thesize
			
			
		elif thetype == "X" or thetype == "S":
			
			posi += thesize
			refposi += thesize
			
		elif thetype == "I" or thetype == "H":
			
			refposi += thesize
			
		elif thetype == "D":
			
			posi += thesize
			#newref += "-"*thesize
			
		if posi <= qstart:
			
			continue
		
		if lastposi <= qstart:
			
			if thetype in ['M','=', 'X','S', 'I','H']:
				
				refstart = refposi - thesize + qstart - lastposi
			else:
				refstart = refposi
				
			corrsize = min(qend,posi) - qstart
			newcigars.append("{}{}".format(corrsize,thetype))
			
			
		if posi >= qend:
			
			if thetype in ['M','=', 'X','S', 'I','H']:
				
				refend = refposi - posi + qend
			else:
				refend = refposi
				
			if lastposi > qstart:
				corrsize = qend - max(qstart, lastposi)
				
				newcigars.append("{}{}".format(corrsize,thetype))
				
			break
		
		if lastposi > qstart:
			
			newcigars.append(cigar)
			
			
	return  [refstart, refend, "".join(newcigars)]


def cigarextend(cigars, refseq, queryseq):
	
	cigars = re.findall(r'\d+[a-zA-Z=]', cigars)
	
	refposi = 0
	posi = 0
	newcigars = []
	
	newref = ""
	for cigar in cigars:
		
		
		thesize = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == "=" or thetype == "M":
			
			
			matchormismatch = 0
			last_i = 0
			newcigar = []
			for i, (r,q) in enumerate(zip(refseq[refposi:(refposi+thesize)], queryseq[posi:(posi+thesize)])):
				
				if matchormismatch ==0 and r.upper() != q.upper():
					
					newcigar.append("{}=".format(i-last_i))
					last_i = i
					matchormismatch = 1
					check = 1
					
				elif matchormismatch ==1 and r.upper() == q.upper():
					
					newcigar.append("{}X{}".format(i-last_i, queryseq[(posi+last_i):(posi+i)]))
					last_i = i
					matchormismatch = 0
					
					
			i = thesize
			if matchormismatch ==0:
				newcigar.append("{}=".format(i-last_i))
			else:
				newcigar.append("{}X{}".format(i-last_i,queryseq[(posi+last_i):(posi+i)]))
				
			cigar = "".join(newcigar)
			
			
			newcigars.append(cigar)
			
			posi += thesize
			
			#newref += refseq[refposi:(refposi+thesize)]
			refposi += thesize
			
			
		elif thetype == "X" or thetype == "S":
			
			newcigars.append(cigar+queryseq[posi:(posi+thesize)])
			posi += thesize
			
			#newref += "X"*thesize
			refposi += thesize
			
		elif thetype == "I" or thetype == "H":
			
			newcigars.append(cigar)
			
			refposi += thesize
			
		elif thetype == "D":
			
			newcigars.append(cigar+queryseq[posi:(posi+thesize)])
			
			posi += thesize
			#newref += "-"*thesize
			
			
	return "".join(newcigars)


def cigarcheck(cigars, refseq, queryseq):
	
	cigars = re.findall(r'\d+[a-zA-Z=]', cigars)
	
	refposi = 0
	posi = 0
	newcigars = []
	
	newref = ""
	for cigar in cigars:
		
		
		thesize = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == "=" or thetype == "M":
			
			
			matchormismatch = 0
			last_i = 0
			newcigar = []
			for i, (r,q) in enumerate(zip(refseq[refposi:(refposi+thesize)], queryseq[posi:(posi+thesize)])):
				
				if r.upper() != q.upper():
					
					print("error",r,q, refposi, posi, thesize)
					print(refseq[refposi:(refposi+thesize)])
					print(queryseq[posi:(posi+thesize)])
					
				#print((r,q) )
				if matchormismatch ==0 and r.upper() != q.upper():
					
					newcigar.append("{}=".format(i-last_i))
					last_i = i
					matchormismatch = 1
					check = 1
					
				elif matchormismatch ==1 and r.upper() == q.upper():
					
					newcigar.append("{}X{}".format(i-last_i, queryseq[(posi+last_i):(posi+i)]))
					last_i = i
					matchormismatch = 0
					
					
			i = thesize
			if matchormismatch ==0:
				newcigar.append("{}=".format(i-last_i))
			else:
				newcigar.append("{}X{}".format(i-last_i,queryseq[(posi+last_i):(posi+i)]))
				
			cigar = "".join(newcigar)
			
			
			newcigars.append(cigar)
			
			newref += refseq[refposi:(refposi+thesize)]
			
			posi += thesize
			
			refposi += thesize
			
			
		elif thetype == "X" or thetype == "S":
			
			newcigars.append(cigar+queryseq[posi:(posi+thesize)])
			posi += thesize
			
			newref += "X"*thesize
			refposi += thesize
			
		elif thetype == "I" or thetype == "H":
			
			
			newcigars.append(cigar)
			
			refposi += thesize
			
		elif thetype == "D":
			
			newcigars.append(cigar+queryseq[posi:(posi+thesize)])
			posi += thesize
			
			newref += "-"*thesize
			
			
	return "".join(newcigars)

def overlap(x,y):
	
	return max(x[1],x[0],y[1],y[0]) - min(x[1],x[0],y[1],y[0]) - abs(x[1] - x[0]) - abs(y[1]-y[0])

def cigarleftcut(cigars, cutsize, qstart = 0):
	
	curr_strd = 1
	curr_rposi = 0
	curr_qposi = 0
	new_rposi = 0
	new_qposi = 0   
	index = 0
	cigarfix = []
	allpaths = []
	for (cigar,seq) in cigars:
		
		curr_size = int(cigar[:-1])
		char = cigar[-1]
		
		if char == '=' or char ==  'M':
			new_rposi = curr_rposi + curr_size
			new_qposi = curr_qposi + curr_size
		elif char == 'D':
			new_rposi = curr_rposi + curr_size
		elif char == 'X':
			new_rposi = curr_rposi + curr_size
			new_qposi = curr_qposi + curr_size
		elif char == 'I':
			new_qposi = curr_qposi + curr_size
			
		if new_rposi > cutsize:
			cigarfix = ["{}{}".format(new_rposi - cutsize,char)]
			if new_qposi > curr_qposi:
				curr_qposi += cutsize - curr_rposi
			index += 1
			break
		
		curr_qposi = new_qposi
		curr_rposi = new_rposi
		index += 1
	if curr_qposi > 0:
		cigarfix = [("{}I".format(curr_qposi),"{}_{}".format(qstart, qstart+curr_qposi))] + cigarfix
	newcigars = cigarfix+cigars[index:]
	
	return newcigars

def connectbreaks(cigars,pathes,qposis,rposis):
	
	lastcigar = ""
	lastrposi = (0,0)
	lastqposi = (0,0)
	lastpath = ""
	
	newcigars = []
	newpathes = []
	newqposis = []
	newrposis = []
	for cigar,path,qposi, rposi in zip(cigars,pathes,qposis,rposis):
		
		if path[0] == ">":
			rgap = rposi[0] - lastrposi[1]
		else:
			rgap = lastrposi[0] - rposi[1]
			
		if path == lastpath and qposi[1] > lastqposi[1]+300 and abs( rgap ) <300 and abs(qposi[0] - lastqposi[1]) < 300:
			
			if rgap > 0:
				cigarfix = "{}D".format(abs(rgap))
				cigar = [cigarfix] + cigar
			elif rgap < 0:
				cigar = cigarleftcut(cigar, abs(rgap), qposi[0])
			newcigars[-1].extend(cigar)
			
			newqposis[-1] = (newqposis[-1][0],qposi[1])
			if path[0] == ">":
				newrposis[-1] = (newrposis[-1][0],rposi[1])
			else:
				newrposis[-1] = (rposi[0], newrposis[-1][1])
		else:
			newcigars.append(cigar)
			newpathes.append(path)
			newqposis.append(qposi)
			newrposis.append(rposi)
			
		lastpath = path
		lastqposi = qposi
		lastrposi = rposi
		
	for i,(cigar,path,qposi, rposi) in enumerate(zip(newcigars,newpathes,newqposis,newrposis)):
		
		leftH,rightH = [],[]
		
		length = int(path.split("_")[-1]) - int(path.split("_")[-2])
		
		if path[0] == '>':
			if rposi[0] > 0:
				leftH = [["{}H".format(rposi[0]), ""]]
			if length -rposi[1] > 0:
				rightH = [["{}H".format(length - rposi[1]), ""]]
		else:
			if rposi[0] > 0:
				rightH = [["{}H".format(rposi[0]), ""]]
			if length - rposi[1] > 0:
				leftH = [["{}H".format(length - rposi[1]), ""]]
				
		newcigars[i] = leftH + cigar + rightH
		
	return newcigars, newpathes, newqposis, newrposis



def getspan_ongraph(pathes, cigars):
	
	qposi = 0 
	qposis = []
	rposis = []
	newcigars = []
	for pathname, cigar in zip(pathes,cigars):
		
		cigar0 = re.findall(r'(\d+[M=XHDI])([A-Za-z]*)', cigar)
		
		gstart,gend = 0,0
		if cigar0[0][0][-1] == 'H':
			gstart = int(cigar0[0][0][:-1])
			cigar0 = cigar0[1:]
		if cigar0[-1][0][-1] == 'H':
			gend = int(cigar0[-1][0][:-1])
			cigar0 = cigar0[:-1]
			
		if cigar[0]=='>':
			gend = gstart + sum([int(x[0][:-1]) for x in cigar0 if x[0][-1] in "M=XD"]+[0])
		else:
			gstart = gend
			gend = gstart + sum([int(x[0][:-1]) for x in cigar0 if x[0][-1] in "M=XD"]+[0])
			
		newcigars.append(cigar0)
		rposis.append((gstart,gend))
		
		qposi_new = qposi + sum([int(x[0][:-1]) for x in cigar0 if x[0][-1] in "M=XI"]+[0])
		qposis.append((qposi, qposi_new))
		
		qposi = qposi_new
		
		
	newcigars,pathes,qposis,rposis = connectbreaks(newcigars,pathes,qposis,rposis)
	
	
	return newcigars, pathes,rposis, qposis

def alignment_polish(name, path, cigar, gregion, qregion, qseq, rseq):
	if cigar == "":
		return name, path, cigar, gregion, qregion
	
	
	cigars = re.findall(r'[><][^><]+', cigar)
	pathes = re.findall(r'[><][^<>]+',path)
	
	cigars, pathes, gregions, qregions  = getspan_ongraph(pathes, cigars)
	for cigar,qregion in zip(cigars,qregions):
		
		for i, cigar0 in enumerate(cigar):
			if type(cigar0) == type(""):
				cigar0 = (cigar0, "")
            
			if cigar0[0].endswith('I') and "_" in cigar0[1] :
				start = qregion[0]
				end = int(cigar0[1].split('_')[1]) + start
				cigar[i] = cigar0[0] + qseq[start:end]
			else:
				cigar[i] = cigar0[0] + cigar0[1]
				
	newcigars= "".join(["{}:{}".format(thename,"".join(cigarstr)) if thename[0] == '>' else "{}:{}".format(thename,"".join(cigarstr)) for thename,cigarstr,rposi,qposi in zip(pathes,cigars,gregions, qregions)])
	
	gregions = ";".join(["{}_{}".format(x[0],x[1]) for x in gregions])
	qregions = ";".join(["{}_{}".format(x[0],x[1]) for x in qregions])
	return name, "".join(pathes), newcigars, gregions, qregions 





def graphcigar(graphfile, queryfile, alignfile):
	
	with open(queryfile, mode = 'r') as f:
		reads = [read.splitlines() for read in f.read().split(">")[1:]]
		queryname = reads[0][0].split()[0]
		queryseq= "".join(reads[0][1:]) 
		
	with open(graphfile, mode = 'r') as f:
		refs = [read.splitlines() for read in f.read().split(">")[1:]]
		refs = {read[0].split()[0]:"".join(read[1:]) for read in refs}
		
		
	allaligns = []
	with open(alignfile, mode = 'r') as f:
		
		lines = f.read().splitlines()
		
		
	for index, line in enumerate(lines):
		
		if len(line) == 0 or line[0]=="@":
			continue
		
		elements = line.strip().split()
		#name, qsize, qstart, qend, path, rsize, rstart, rend = elements[:8]
		
		identity =float( [x for i,x in enumerate(elements) if x[:5] == "PI:f:"][0][5:])
		match = int( [x for i,x in enumerate(elements) if x[:5] == "AS:i:"][0][5:] ) 
		if identity < 90 or match < 100 or elements[4] != '255':
			continue
		
		strand, qstart, cigar = elements[1], int(elements[3])-1, elements[5]
		
		#aligns = findkmeraligns(cigar)
		alignsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XD]', cigar)]+[0])
		aligns = [[0, alignsize, match]]
		allaligns.extend(  [[x[0]+ qstart, x[1] + qstart]+[ x[2],index] for x in aligns])
		
	if len(allaligns) == 0:
		return queryname, "", "", "", ""
	
	
	alignregions = [x for x in graphicalign(allaligns) if x[2]- x[1]]
	
	pathes = []
	for region in alignregions:
		
		line = lines [region[0]] 
		
		elements = line.strip().split()
		rname, strand, qstart, cigar = elements[0], elements[1], int(elements[3]) -1, elements[5]
		strand = 1 if strand == '0' else -1
		
		
		newcigar = cigar_findqrange(cigar, region[1] - qstart, region[2] - qstart) 
		
		refseq = refs[rname]
		qseq = queryseq
		if strand == -1:
			refseq = makereverse(refseq)
			
		rseq = refseq[newcigar[0]:]
		qseq = qseq[ region[1]:]
		
		newcigar[-1] = cigarextend(newcigar[-1], rseq, qseq )
		
		
		pathes.append(region + [rname, strand] + newcigar )
		
		
	fullpath = ""
	fullcigar = []
	reftemplate = []
	last_rname = ""
	last_qend = 0
	last_rend = 0
	last_strand = ""
	last_rclip = ""
	
	qstarts = []
	pathrange = []
	for path in pathes:
		
		index, qstart, qend, rname, strand, rstart, rend, cigar = path
		
		if rend == -1:
			rend = len(refs[rname])
			
		qstarts.append("{}_{}".format(qstart,qend))
		cigar = cigar.replace("H","I")
		
		if last_rname == rname and last_strand == strand and rstart >= last_rend:
			
			
			fullcigar[-1] = ""
			if rstart - last_rend:
				fullcigar.append(str(rstart - last_rend) + "I")
				
				
			if qstart - last_qend:
				fullcigar.append(str(qstart - last_qend) + "D" + queryseq[last_qend:qstart])
				
			fullcigar.extend(re.findall(r'\d+[a-zA-Z=][actgACTG]*', cigar))
			
			pathrange[-1][1] = rend
			qstarts = qstarts[:-1]
			if len(qstarts):
				qstarts[-1] = "{}_{}".format(qstarts[-1].split("_")[0],qend) 
				
		
			"""
			if strand == 1:
				reftemplate.append(refs[rname][last_rend:rend])
			else:
			reftemplate.append(makereverse(refs[rname])[last_rend:rend])
			"""
		
		else:
			if strand == 1: 
				fullcigar.append(">")
			else:
				fullcigar.append("<")
				
			if rstart:
				fullcigar.append(str(rstart) + "H")
				pass
				
			if qstart - last_qend:
				fullcigar.append(str(qstart - last_qend) + "D" + queryseq[last_qend:qstart])
				
		
			fullcigar.extend(re.findall(r'\d+[a-zA-Z=][actgACTG]*', cigar))
		
			if strand == 1:
				fullpath += ">{}".format(rname)
			else:
				fullpath += "<{}".format(rname)
				
				
			pathrange.append([rstart, rend])
		
			"""
			if strand == 1:
				reftemplate.append(refs[rname])
			else:
				reftemplate.append(makereverse(refs[rname]))
			"""
		
		if len(refs[rname]) - rend:
			fullcigar.append(str(len(refs[rname]) - rend) + "H")
			
		last_qend = qend
		last_rend = rend
		last_rname = rname
		last_strand = strand
	
	qstarts.append("{}_{}".format(last_qend,len(queryseq)))
	fullcigar = "".join(fullcigar)
	#reftemplate = "".join(reftemplate)
	#cigarcheck(fullcigar, reftemplate , queryseq)
	fullcigar = fullcigar.replace('D', ';').replace('I', 'D').replace(';', 'I')
	#print(queryname, fullpath, fullcigar, ";".join([str(x[0])+"_"+str(x[1]) for x in pathrange]), ";".join(list(map(str,qstarts))))
	return alignment_polish( queryname, fullpath, fullcigar, ";".join([str(x[0])+"_"+str(x[1]) for x in pathrange]), ";".join(list(map(str,qstarts))), queryseq, refseq )



def main(args):
	
	queryname, fullpath, fullcigar, pathranges, qranges= graphcigar(args.ref, args.query, args.input)
	with open(args.output, mode = 'w') as w:
		w.write("{}\t{}\t{}\t{}\t{}\n".format(queryname, fullpath, fullcigar, pathranges, qranges))
		
		
		
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	
	parser.add_argument("-r", "--ref", help="path to output file", dest="ref", type=str,required=True) 
	parser.add_argument("-q", "--query", help="path to output file", dest="query", type=str,required=True)  
	parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()

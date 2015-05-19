#!/usr/bin/env python

import sys, math, getopt, random, string, os, re
import networkx as nx

class Tree:
	def __init__(self, tree_file, genomeToLocusFile):
		self.tree_file = tree_file
		self.genomeToLocusFile = genomeToLocusFile
		self.genomeToLocus = {}
		self.locusToGenome = {}
		self.tree_string = ""
		self.tree = nx.Graph()
		self.ancestral = []
		self.rooted_tree = nx.DiGraph()
		
	def readTree(self):
		tree = open(self.tree_file,'r').read()
		#~ print tree
		trees = re.split("\)[\d\.]*\d+:",tree)
		jstring = "):"
		#~ print len(trees),trees
		tree = jstring.join(trees)
		#~ print tree
		#~ sys.exit()
		self.tree_string = tree
		
	def readGenomeToLocusFile(self):
		tags = open(self.genomeToLocusFile,'r').readlines()
		for t in tags:
			t = t.rstrip()
			l = t.split()
			if not len(l) > 1:
				continue
			self.genomeToLocus[l[0]] = l[1]
			self.locusToGenome[l[1]] = l[0]
		
	def codeGenomeID(self, genome):
		tag = ''
		if genome in self.genomeToLocus:
			tag = self.genomeToLocus[genome]
		else:
			tag = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
			while tag in self.locusToGenome:
				tag = ''.join(random.choice(string.ascii_uppercase) for x in range(3))
		#~ print tag, genome
		self.genomeToLocus[genome] = tag
		self.locusToGenome[tag] = genome
		return tag
		
	def reregisterGenomeID(self, id, newChildren):
		oldGenome = self.locusToGenome[id]
		newGenome = ";".join(newChildren)
		del self.genomeToLocus[oldGenome]
		#~ print id, newGenome
		self.genomeToLocus[newGenome] = id
		self.locusToGenome[id] = newGenome
		
		
	def parseTree(self):
		myTreeString = self.tree_string[0:-1]
		myExcisableRegions = self.indexParens(myTreeString)
		max_weight = 0.0
		while len(myExcisableRegions) > 0:
			myExcisableRegions = sorted(myExcisableRegions, key=lambda reg: reg[3], reverse=True)
			newSpeciesLocus = ""
			#~ print "excisable", len(myExcisableRegions)
			for mer in myExcisableRegions:
				region = mer[2][1:-1]
				#~ print len(myTreeString), len(mer[2])

				genomes = region.split(",")
				child_nodes = []
				#~ print "genomes", len(genomes)
				#~ print genomes
				for g in genomes:
					#~ print g
					nome = g.split(":")[0]
					dist = float(g.split(":")[1])
					locus = ""
					if nome in self.locusToGenome:
						locus = nome
					else:
						locus = self.codeGenomeID(nome)

					child_nodes.append((locus,dist))
				same_length = 0
				if len(child_nodes) ==1:
					print "child nodes ==1"
					sys.exit()
				if len(mer[2]) == len(myTreeString) and len(genomes)==2:
					same_length =1
					tnodes = []
					for cn in child_nodes:
						self.tree.add_node(cn[0])
						tnodes.append(cn[0])
					self.tree.add_edge(tnodes[0],tnodes[1],weight=child_nodes[0][1])
					myTreeString = region
					#~ print "BREAKING"
					break
					
				has_dist = 0
				while len(child_nodes) > 1 and same_length == 0:
					#~ print "WHILE", len(child_nodes)
					if len(child_nodes)>2 and child_nodes[0][1]>0.0:
						has_dist = 1
						#~ print "has_dist", has_dist
						break
					tkids = []
					#~ print len(child_nodes)
					tkids.append(child_nodes.pop(0))
					tkids.append(child_nodes.pop(0))
					#~ print len(child_nodes)
					tnodes = []
					for tk in tkids:
						self.tree.add_node(tk[0])
						tnodes.append(tk[0])
					tnodes.sort()
					newSpeciesLocus = self.codeGenomeID(";".join(tnodes))
					#~ print tkids,";".join(tnodes), newSpeciesLocus
					self.tree.add_node(newSpeciesLocus)
					weight = 0.0
					for tk in tkids:
						weight = tk[1]
						self.tree.add_edge(newSpeciesLocus,tk[0],weight=tk[1])
						#~ print tk[0],self.tree.edge[tk[0]]
					#~ print newSpeciesLocus, self.tree.edge[newSpeciesLocus]
					child_nodes.append((newSpeciesLocus,weight))
				if has_dist ==1:
					#~ print myTreeString
					min_pair = (child_nodes[0], child_nodes[1])
					min_dist = child_nodes[0][1]+child_nodes[1][1]
					#~ print min_pair, min_dist
					subPaths = {}
					for c in child_nodes:
						#~ print c
						if self.tree_string.find(self.locusToGenome[c[0]])>-1:
							#~ self.genomes.append(self.locusToGenome[c[0]])
							self.tree.add_node(c[0])
						paths = nx.shortest_path_length(self.tree,c[0],None,'weight')
						longest_path = 0.0
						#~ print paths
						for p in paths:
							if p == c[0]:
								continue
							if paths[p] > longest_path:
								longest_path = paths[p]
								#~ print longest_path
						longest_path += c[1]
						subPaths[c[0]] = longest_path
					for c in child_nodes:
						for n in child_nodes:
							if c==n:
								continue
							myDist = subPaths[c[0]]+subPaths[n[0]]
							#~ print c[0], subPaths[c[0]], n[0], subPaths[n[0]], myDist
							if myDist < min_dist:
								min_dist = myDist
								min_pair = (c, n)
								#~ print min_pair, min_dist
					tkids = []
					#~ print len(child_nodes), child_nodes
					tkids.append(child_nodes.pop(child_nodes.index(min_pair[0])))
					tkids.append(child_nodes.pop(child_nodes.index(min_pair[1])))
					#~ print len(child_nodes), child_nodes
					tnodes = []
					for tk in tkids:
						tnodes.append(tk[0])
					tnodes.sort()
					newSpeciesLocus = self.codeGenomeID(";".join(tnodes))
					#~ print tkids,";".join(tnodes), newSpeciesLocus
					self.tree.add_node(newSpeciesLocus)
					weight = child_nodes[0][1]
					for tk in tkids:
						self.tree.add_edge(newSpeciesLocus,tk[0],weight=tk[1])
					self.tree.add_edge(newSpeciesLocus,child_nodes[0][0],weight=child_nodes[0][1])
					
					newSpeciesLocus = child_nodes[0][0]
					#~ newSpeciesLocus = "("+",".join(forReplace)+")"
				#~ print newSpeciesLocus
				myTreeString = myTreeString.replace(mer[2],newSpeciesLocus)
				#~ print myTreeString
				child_nodes = []
			myExcisableRegions = self.indexParens(myTreeString)
			
	def trimTaxa(self, taxaToKeep):
		ancestral = set([])
		for n in self.tree.nodes():
			n_nome = self.locusToGenome[n]
			if n_nome.find(";")>-1:
				ancestral.add(n)
				continue
			if not n_nome in taxaToKeep:
				#~ print "removed", n_nome
				self.tree.remove_node(n)
		for n in self.tree.nodes():
			n_nome = self.locusToGenome[n]
			if n_nome.find(";") > -1 and len(self.tree[n])<2:
				self.tree.remove_node(n)
		badAncestral = 1
		while (badAncestral):
			badAncestral = 0
			for n in self.tree.nodes():
				if n in ancestral and len(self.tree[n])==2:
					badAncestral+=1
					myedges = []
					myweights = []
					for e in self.tree[n]:
						myedges.append(e)
						myweights.append(self.tree[n][e]['weight'])
					newweight = myweights[0]+myweights[1]
					#~ print myedges,myedges[0],myedges[1], newweight
					self.tree.add_edge(myedges[0],myedges[1],weight=newweight)
					self.tree.remove_node(n)
		
	def rootByMidpoint(self):
		paths = nx.shortest_path_length(self.tree,None,None,'weight')
		path_count = nx.shortest_path_length(self.tree,None,None,None)
		#~ print len(paths)
		longest_path = 0.0
		longest_path_count = 0
		longest_pair = None
		path_dist = {}
		for p in paths:
			path_dist[p] = {}
			for r in paths[p]:
				if r==p:
					path_dist[p][p] = 0.0
					continue
				if r in path_dist[p]:
					if paths[p][r] < path_dist[p][r]:
							path_dist[p][r] = paths[p][r]
				else:
					path_dist[p][r] = paths[p][r]
				if paths[p][r] > longest_path:
					longest_pair = (p, r)
					longest_path = paths[p][r]
					longest_path_count = path_count[p][r]
				elif paths[p][r] == longest_path:
					if path_count[p][r] > longest_path_count:
						longest_path_count = path_count[p][r]
						longest_pair = (p,r)					
				elif longest_pair == None:
					longest_pair = (p,r)
					longest_path_count = path_count[p][r]
		path = nx.shortest_path(self.tree, longest_pair[0],longest_pair[1])
		mid = longest_path/2.0
		cur_node = path.pop()
		cur_length = 0.0
		root_edge = None
		while path:
			nextNode = path.pop()
			cur_length+= self.tree[cur_node][nextNode]['weight']
			if cur_length > mid:
				root_edge = (cur_node,nextNode)
				break
			cur_node = nextNode
		if root_edge == None:
			#this basically means all of the sequences are identical
			root_edge = longest_pair
		print root_edge, mid*2.0
		re_weight = (self.tree.edge[root_edge[0]][root_edge[1]]['weight'])/2.0
		print "root edge", root_edge
		#~ print root_edge[0],self.tree.edge[root_edge[0]]
		#~ print root_edge[1],self.tree.edge[root_edge[1]]
		self.tree.remove_edge(root_edge[0], root_edge[1])
		self.tree.add_node("root")

		mod_root_edge = []
		for re in root_edge:
			#~ print re, len(self.tree.edges(re)), self.tree.edges(re)
			if not len(self.tree.edges(re))==1:
				#~ print "edge from root to", re
				self.tree.add_edge("root",re,weight=re_weight)
				#~ print re, self.tree.edge[re]
				mod_root_edge.append(re)
			else:
				#~ print "single", re, self.tree.edge[re]
				new_rooter = ""
				for e in self.tree.edge[re]:
					if not e == re:
						new_rooter = e
				myWeight = re_weight+self.tree.edge[new_rooter][re]['weight']
				self.tree.remove_edge(new_rooter, re)
				self.tree.remove_node(re)
				self.tree.add_edge("root",new_rooter,weight=myWeight)
				#~ print "remove", re
				#~ print "new rooter", new_rooter, self.tree.edge[new_rooter]
				mod_root_edge.append(new_rooter)
		print "mod_edge",mod_root_edge
		self.rootTree(mod_root_edge)

	def rootTree(self, root_edge):
		self.rooted_tree = self.tree.to_directed()
		paths = nx.shortest_path(self.rooted_tree,source="root")
		for e in self.tree.edges():
			#~ print e
			e0 = len(paths[e[0]])
			e1 = len(paths[e[1]])
			if e0<e1:
				self.rooted_tree.remove_edge(e[1],e[0])
			else:
				self.rooted_tree.remove_edge(e[0],e[1])
		
		#~ print self.rooted_tree.edge['UIZ']
		for n in self.rooted_tree.nodes():
			if len(self.rooted_tree.out_edges(n)) > 0:
				self.ancestral.append(n)
		print "ancestral", len(self.ancestral)
		print "all nodes", len(self.tree.nodes())
		
	def indexParens(self, tree_string):
		left_parens = []
		right_parens = []
		parens = []
		left = tree_string.find("(",0)

		while left < len(tree_string)-1 and left > -1:
			left_parens.append(left)
			parens.append((left,"("))
			left = tree_string.find("(",left+1)

		right = tree_string.find(")",0)
		while right < len(tree_string) and right > -1:
			right_parens.append(right)
			parens.append((right,")"))
			right = tree_string.find(")",right+1)
		parens = sorted(parens, key=lambda tup: tup[0])
		myRegions = []
		i = 1
		while i<len(parens):
			p = parens[i]
			if p[1] == ")" and parens[i-1][1] == "(":
				l = parens[i-1][0]
				r = p[0]+1
				length = r-l
				myRegions.append((l,p[0],tree_string[l:r], length))
			i+=1
		return myRegions
		
	def writeLocusTagFile(self):
		tag_out = open(self.genomeToLocusFile, 'w')
		for g in self.genomeToLocus:
			line = "\t".join([g, self.genomeToLocus[g]])+"\n"
			tag_out.write(line)
		tag_out.close()
		print "Wrote locus tags to locus_tag_file.txt"
		
		
	def toNewickLabels(self):
		graph = self.rooted_tree.copy()
		leaf = set([])
		unproc = set([])
		text = {}
		subs = nx.connected_component_subgraphs(self.tree)
		#~ print "subs",len(subs)
		#~ sys.exit()
		len_succ = set([])
		for n in graph.nodes():
			mySucc = graph.successors(n)
			if len(mySucc) == 0:
				leaf.add(n)
				text[n] = n
			elif len(mySucc) == 1:
				#~ print n, mySucc
				#~ print graph.edge[n]
				#~ for m in mySucc:
					#~ print graph.edge[m]
				sys.exit()
			else:
				len_succ.add(len(graph.successors(n)))
				unproc.add(n)
		last_string = ""
		while len(unproc) > 0:
			#~ print len(unproc), unproc
			for_del = set([])
			for u in unproc:
				u_succ = set(graph.successors(u))
				if len(u_succ.intersection(leaf))>1:
					my_text = []
					for i in u_succ.intersection(leaf):
						my_text.append(text[i]+":"+str(graph.edge[u][i]['weight']))
						leaf.remove(i)
					text[u] = "("+",".join(my_text)+")"+u
					last_string = text[u]
					for_del.add(u)
					leaf.add(u)
			unproc -= for_del
			#~ print len(leaf)
		#~ return "("+last_string+");"
		return last_string+";"

	def toNewickNoLabels(self):
		graph = self.rooted_tree.copy()
		leaf = set([])
		unproc = set([])
		text = {}
		subs = nx.connected_component_subgraphs(self.tree)
		#~ print "subs",len(subs)
		#~ sys.exit()
		len_succ = set([])
		for n in graph.nodes():
			mySucc = graph.successors(n)
			if len(mySucc) == 0:
				leaf.add(n)
				text[n] = n
			elif len(mySucc) == 1:
				#~ print n, mySucc
				#~ print graph.edge[n]
				#~ for m in mySucc:
					#~ print graph.edge[m]
				sys.exit()
			else:
				len_succ.add(len(graph.successors(n)))
				unproc.add(n)
		last_string = ""
		while len(unproc) > 0:
			#~ print len(unproc), unproc
			for_del = set([])
			for u in unproc:
				u_succ = set(graph.successors(u))
				if len(u_succ.intersection(leaf))>1:
					my_text = []
					for i in u_succ.intersection(leaf):
						my_text.append(text[i]+":"+str(graph.edge[u][i]['weight']))
						leaf.remove(i)
					text[u] = "("+",".join(my_text)+")"
					last_string = text[u]
					for_del.add(u)
					leaf.add(u)
			unproc -= for_del
			#~ print len(leaf)
		return last_string+";"


	def calcMostEdgesToLeaves(self,unprocN,leaf,TG):
		mostLeaves = 0
		retNode = None
		for n in unprocN:
			succ = set(TG.successors(n))
			e_inter = succ.intersection(leaf)
			e_count = len(e_inter)
			#~ print n,e_count, e_inter
			e_count = 0
			for e in TG[n]:
				for l in leaf:
					if e == l:
						e_count += 1
			if e_count > mostLeaves:
				#~ print "edge to leaf",n, e_count, mostLeaves
				mostLeaves = e_count
				retNode = n
		return (retNode,mostLeaves)
		
	def makeClassSplits(self):
		ancestor_to_children = {}
		for a in self.ancestral:
			ancestor_to_children[a] = set([])
			leaves = []
			kids = self.getLeaves(a, leaves)
			for k in kids:
				ancestor_to_children[a].add(self.locusToGenome[k])
		return ancestor_to_children
		
	def getLeaves(self, ancestral_node, leaves):
		for c in self.rooted_tree.successors(ancestral_node):
			if len(self.rooted_tree.out_edges(c)) == 0:
				leaves.append(c)
			else:
				self.getLeaves(c, leaves)
		return leaves
			
	def generateAncestryMatrix(self, out_file):
		ancToLeaf = {}
		myLeaves = set([])
		for a in self.ancestral:
			leaves = []
			leaves = self.getLeaves(a, leaves)
			ancToLeaf[a] = leaves
			#~ print a, len(leaves), leaves
			for l in leaves:
				myLeaves.add(l)
		myLeaves = sorted(myLeaves)
		out  = open(out_file,'w')
		header = ""
		for m in myLeaves:
			header = header + "\t"+self.locusToGenome[m]
		out.write(header+"\n")
		for a in self.ancestral:
			line = "node_"+a
			for m in myLeaves:
				if m in ancToLeaf[a]:
					line = line+"\t1"
				else:
					line = line+"\t0"
			out.write(line+"\n")
		out.close()
			
		
			
		
def usage():
	print """
	python tree_poke.py [options]
	-t, --tree [file]
		where [file] is a species tree in newick format without bootstrap values
	-n, --node
		if flagged, ancestral nodes will be labelled.  Mappings will be stored in the locus file
	-l, --locus [file]
		where [file] is the location of the locus mappings, stored in locus_tag_file.txt by default
	-r, --root
		if flagged, tree will be rooted by the midpoint based on branch length distances
		CURRENTLY ALL TREES WILL BE ROOTED. OOPS?
	-k, --keep [file]
		where [file] is a one taxa per line file of taxa to retain in the tree.  The topology of the existing tree will be retained.
		No labels will be retained or reported.
	-m, --matrix [file]
		where [file] will be the output location of a binary matrix that contains ancestral to leaf node mappings.
	-o, --outTree [file]
		where [file] will be the output location of your modified tree. Tree will still be printed to standard out.
	-h, --help
		prints this and exits
	"""
def main(argv):
	tree_file = ""
	locus_tag = "locus_tag_file.txt"
	root = 1 #always on... need to develop
	labels = 0
	taxa_to_keep = ""
	matrix = ""
	outTree = ""
	
	try:
		opts, args = getopt.getopt(argv, "t:l:k:m:o:rnh",["tree=","locus=","keep=","matrix=","outTree=""root","node","help"])

	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-t","--tree"):
			tree_file = arg
		elif opt in ("-l","--locus"):
			locus_tag = arg
		elif opt in ("-r","--root"):
			root = 1
		elif opt in ("-n","--node"):
			labels = 1
		elif opt in ("-k","--keep"):
			taxa_to_keep = arg
		elif opt in ("-m","--matrix"):
			matrix = arg
		elif opt in ("-o","--outTree"):
			outTree = arg
		elif opt in ("-h","--help"):
			sys.exit(usage())
	
	myTree = Tree(tree_file, locus_tag)
	if locus_tag in os.listdir("."):
		myTree.readGenomeToLocusFile()
	myTree.readTree()
	myTree.parseTree()
	if len(taxa_to_keep)>0:
		keepers =  open(taxa_to_keep,'r').readlines()
		keep = []
		for k in keepers:
			k = k.rstrip()
			keep.append(k)
		myTree.trimTaxa(keep)
		labels = 0
	if root>0:
		myTree.rootByMidpoint()
		
	if len(matrix) > 0:
		myTree.generateAncestryMatrix(matrix)
	myNewTree = ""
	if labels > 0:
		myNewTree = myTree.toNewickLabels()
	else:
		myNewTree = myTree.toNewickNoLabels()
	#~ print myNewTree

	for l in myTree.locusToGenome:
		l_nome = myTree.locusToGenome[l]
		if l_nome.find(";")>-1:
			continue
		myNewTree = myNewTree.replace(l, l_nome)
	if len(outTree)>0:
		mout = open(outTree,'w')
		mout.write(myNewTree)
		mout.close()
	else:
		print myNewTree
	myTree.writeLocusTagFile()
	
if __name__ == "__main__":
	if len(sys.argv) == 1:
		sys.exit(usage())
	else:
		main(sys.argv[1:])
#!/usr/local/bin/python3
"""
Input:  *BLAST results run through blast_formatter with 
			-outfmt '6 qacc sacc gi staxid pident length mismtach gapopen qstart qend sstart send evalue bitscore'
			NOTE: Order matters in this output

		*names.dmp and nodes.dmp from new_taxdump.tar.gz

		*OTU table from DAMe's tabulateSumaclust.py 

Output: *A combined table with the top BLAST hit for each OTU
			- Top hit(s) are selected from hits with the highest bitscore AND percent identity
			- Consensus taxonomy is created by looking for lowest taxonomy rank where all top hits agree on assignment
				and subsequent ranks are changed to "NA"
			- PID cutoffs used are species: 98.0%, family: 96.5%, order: 95.0%
				(these values can be edited by changing the variables in the consensusHit() function)

"""
import sys
import os
import argparse
import copy

#Argument parsing
parser = argparse.ArgumentParser("BLAST OTU merging tool")
#NCBI Files
parser.add_argument("-n", "--names", dest="names", type=str, help="NCBI names.dmp file.", required=True)
parser.add_argument("-s", "--nodes", dest="nodes", type=str, help="NCBI nodes.dmp file.", required=True)
parser.add_argument("-m", "--merged", dest="merged", type=str, help="NCBI merged.dmp file.", required=True)
#User generated files
parser.add_argument("-b", "--fblast", dest="fblast", type=str, help="Formatted blast table", required=True)
parser.add_argument("-t", "--otutab", dest="otutab", type=str, help="OTU table from tabulateSumaclust.py", required=True)
#Optional
parser.add_argument("-a", "--alltax", dest="alltax", action="store_true", help="Outputs separate file with taxonomy of all hits")
parser.add_argument("-x", "--taxout", dest="taxout", type=str, help="The all taxonomy output filename.", required=False, default="out_alltax.tsv")
parser.add_argument("-o", "--output", dest="output", type=str, help="The output filename.", required=False, default="out.tsv")
args = parser.parse_args()
if os.path.isfile(args.output):
	print("Warning: Output file {} already exists. Overwriting file".format(args.output), file=sys.stderr)
if args.alltax and os.path.isfile(args.taxout):
	print("Warning: All taxonomy output file {} already exists. Overwriting file".format(args.taxout), file=sys.stderr)

#Traces the path of a given taxid to root node, compiling taxonomy names along the way and returning a dict
def parseToRoot(ti):
	taxonomy = {
		"phylum" : "NA", 
		"class" : "NA",
		"order" : "NA",
		"family" : "NA",
		"genus" : "NA",
		"species" : "NA"
	}
	#While not at the root node
	while (ti != "1"):
		#Adds names to taxonomy dictionary while tracing rootward
		try:
			level_name = name_dict[ti]
			(ti,level) = parent_dict[ti]
			if level in taxonomy:
				taxonomy[level] = level_name
		#Catches if a taxid isn't present in names.dmp or nodes.dmp 
		except KeyError as e:
			print("Taxid: {} not found in either names.dmp or nodes.dmp, looking in merged.dmp".format(ti), file=sys.stderr)
			#Looks through merged.dmp for 
			for tokens in lineSplitter(args.merged):
				found = False
				#If the taxid is found
				if tokens[0] == ti:
					ti = tokens[2]
					found = True
					print("Taxid: {} match found in merged.dmp".format(ti), file=sys.stderr)
					break
			if not found:
				print("Taxid: {} not found in merged.dmp either! Try downloading an updated taxdump".format(ti), file=sys.stderr)
				break
			
	return taxonomy

#Looks for agreement in taxonomy of multiple hits to return the correct taxonomy
def consensusHits(hl):
	comp = hl[0][3]
	conl = ""
	#Starting with species, checkes taxonomy level of each maxhit see if they all agree on assignment (and aren't NA)
	for l in tax_order:
		#This one-liner is a bit much. Maybe break into a couple steps for legibility?
		#Fancy map(lambda) function just returns element 3 (taxid) of each list in a list of lists
		if fulltax_dict[comp][l] != "NA" and all(fulltax_dict[t][l] == fulltax_dict[comp][l] for t in map(lambda v: v[3], hl)):
			conl = l
			break

	#Changes taxonomy dictionary to NAs up until the consensus level
	returnd = copy.copy(fulltax_dict[comp])
	for l in tax_order:
		if conl == l:
			break
		returnd[l] = "NA"

	return returnd

#Checks for PID cutoffs
def pidCutoff(pid, at):

	if pid < order_cutoff:
		for l in tax_order[:4]:
			at[l] = "NA"
	elif pid < family_cutoff:
		for l in tax_order[:2]:
			at[l] = "NA"
	elif pid < species_cutoff:
		at["species"] = "NA"

	return at

#Fancy generator function for tokenizing file lines
def lineSplitter(f):
	with open(f) as inputf:
		for line in inputf:
			yield line.rstrip().split("\t")

#Due to dictionary's unordered nature, this is necessary a couple places
tax_order = ["species", "genus", "family", "order", "class", "phylum"]

#If PID is less than these the values, the associated rank and lower will be changed to "NA"
species_cutoff = float(98.0)
family_cutoff = float(96.5)
order_cutoff = float(95.0)


#Makes a dict of name and taxid {taxid : name}
name_dict = {}
for tokens in lineSplitter(args.names):
	#Ignores synonyms, blast names, etc.
	if tokens[6] == "scientific name":
		name_dict[tokens[0]] = tokens[2]
print("Done loading taxid to name data - {}".format(len(name_dict)), file=sys.stderr)

#Makes the parent node dictionary {taxid : (parent_taxid, taxonomic_level)}
parent_dict = {}
for tokens in lineSplitter(args.nodes):
	parent_dict[tokens[0]] = (tokens[2], tokens[4])
print("Done loading parent data - {}".format(len(parent_dict)), file=sys.stderr)

'''
NOTE: The next couple steps could be combined to improve runtime by only accessing fblast once,
but for legibility's sake I've split them into two reading steps (three dependent on the "-a" optional flag). 
This part of the pipeline definitely isn't the bottleneck, speed gains were relatively negligible from reading once.
'''

#Makes a dictionary of max hit scores for each qacc {qacc : (maxbitscore, maxPID)}
maxhitscore_dict = {}
for tokens in lineSplitter(args.fblast):

	#If this qacc isn't in the dictionary, add it
	if tokens[0] not in maxhitscore_dict:
		maxhitscore_dict[tokens[0]] = (float(tokens[-1]), float(tokens[4]))
	#Aboslute monster of a conditional statement. In plain english, it's checking the dictionary to determine:
	#If this qacc has a larger bit score OR if the bit score is the same and the PID is larger, then update dictionary
	elif float(tokens[-1]) > maxhitscore_dict[tokens[0]][0] or (float(tokens[-1]) == maxhitscore_dict[tokens[0]][0] and float(tokens[4]) > maxhitscore_dict[tokens[0]][1]):
		maxhitscore_dict[tokens[0]] = (float(tokens[-1]), float(tokens[4]))

#Makes YET ANOTHER dictionary of BLAST hits that are equal to maxhitscore {qacc : list(blast_result)}
maxhits_dict = {}
for tokens in lineSplitter(args.fblast):

	#If it is one of the maxhitscore results
	if float(tokens[-1]) == maxhitscore_dict[tokens[0]][0] and float(tokens[4]) == maxhitscore_dict[tokens[0]][1]:
		#Adds full line to a list in the dictionary
		if tokens[0] not in maxhits_dict:
			maxhits_dict[tokens[0]] = [tokens]
		else:
			maxhits_dict[tokens[0]].append(tokens)

#Makes a dictionary of dictionaries (I <3 dictionaries) of taxids and their full taxonomy {taxid : {full taxonomy}}
fulltax_dict = {}
#If the optional flag is set for the taxonomy of each hit to be returned
if args.alltax:
	allhits = []
	taxidset = set()
	for tokens in lineSplitter(args.fblast):
		allhits.append(tokens)
		taxidset.add(tokens[3])
	#Adds full taxonomy dictionary for each taxid to fulltax_dict
	for t in taxidset:
		fulltax_dict[t] = parseToRoot(t)
	#Writes out the file with all taxonomy information
	with open(args.taxout, "w") as outputf:
		#Header
		outputf.write("qacc\tsacc\tgi\tstaxid\tpiden\tlengt\tmismtach\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore")
		for l in tax_order[::-1]:
			outputf.write("\t{}".format(l))
		outputf.write("\n")
		for hit in allhits:
			#BLAST hit info
			outputf.write(hit[0])
			for elem in hit[1:]:
				outputf.write("\t{}".format(elem))
			#Taxonomy
			for l in tax_order[::-1]:
				outputf.write("\t{}".format(fulltax_dict[hit[3]][l]))
			outputf.write("\n")

#If all taxonomy file isn't needed, only gets fulltax info of maxhits and doesn't bother accessing fblast again
else:
	taxidset = set()
	for q in maxhits_dict.values():
		for h in q:
			taxidset.add(h[3])
	#Adds full taxonomy dictionary for each taxid to fulltax_dict
	for t in taxidset:
		fulltax_dict[t] = parseToRoot(t)

#Makes a dictionary (last one I promise) to associate the formatted BLAST hit with OTU table
final_dict = {}
for q in maxhits_dict:

	#If there's only one hit with maxhitscore, assigns it's taxonomy
	if len(maxhits_dict[q]) == 1:
		append_tax = fulltax_dict[maxhits_dict[q][0][3]]
	#Else it gets a consensus taxonomy
	else:
		append_tax = consensusHits(maxhits_dict[q])

	#Checks for pid cutoffs
	append_tax = pidCutoff(float(maxhits_dict[q][0][4]), append_tax)

	#Adds all BLAST info from first maxhit (somewhat arbitrary choice) to dictionary
	final_dict[q] = maxhits_dict[q][0]
	#Removes OTU name (since it will be added from OTU table)
	final_dict[q].pop(0)
	#Adds the number of max hits as a column right before taxonomy info
	final_dict[q].append(len(maxhits_dict[q]))
	#Adds taxonomy info
	for l in tax_order[::-1]:
		final_dict[q].append(append_tax[l])

#Final appending to OTU table and writing
with open(args.otutab) as inputf, open(args.output, "w") as outputf:
	
	#Header
	outputf.write(inputf.readline().rstrip())
	outputf.write("\tsacc\tgi\tstaxid\tpiden\tlengt\tmismtach\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tnummaxhits")
	for l in tax_order[::-1]:
		outputf.write("\t{}".format(l))
	outputf.write("\n")

	#Subsequent rows
	for line in inputf:
		#Gets OTU name from OTU table
		otu_name = line.split("\t")[0]
		#Adds OTU table info
		outputf.write(line.rstrip())
		
		#Try block necessary because some OTUs won't return blast hits
		try:
			#Looks up OTU name in final_dict to write
			for item in final_dict[otu_name]:
				outputf.write("\t{}".format(item))
			outputf.write("\n")
		except KeyError as e:
			print("{} returned no BLAST hits".format(otu_name), file=sys.stderr)
			outputf.write("\n")

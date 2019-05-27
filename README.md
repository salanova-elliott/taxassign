# taxassign
Assigns taxonomy from BLAST searches and 

## Input:  
* BLAST results run through `blast_formatter` with 
			- outfmt '6 qacc sacc gi staxid pident length mismtach gapopen qstart qend sstart send evalue bitscore'
			> NOTE: Order matters in this output
* names.dmp and nodes.dmp from new_taxdump.tar.gz
* OTU table from DAMe's tabulateSumaclust.py 

## Output: 

* A combined table with the top BLAST hit for each OTU
  * Top hit(s) are selected from hits with the highest bitscore AND percent identity
  * Consensus taxonomy is created by looking for lowest taxonomy rank where all top hits agree on assignment
				and subsequent ranks are changed to "NA"
  * PID cutoffs used are species: 98.0%, family: 96.5%, order: 95.0%
				(these values can be edited by changing the variables in the consensusHit() function)

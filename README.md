# Genome-wide Ortholog Identification and Visualization
Our goal with this assignment is to identify (coding) gene orthologs between genome pairs (among 3 ray-finned fish species) using a reciprocal best BLAST hit approach. The biological objective is to learn something about the evolution of genome structure in these fish species. You should complete the BLAST-ing on Talapas and the data formatting/plotting using R. The assignment is worth 75 points total.
# Finding and enumerating putative orthologs via RBBH
1.	On Talapas in /projects/bgmp/shared/Bi623/assign2, find the 3 FASTA files containing protein sequences (i.e. XXX.pep.all.fa.gz) from the spotted gar genome, the green spotted puffer genome, and the threespine stickleback genome. Although it would be ideal to select only one isoform per gene from each sequence set, it’s not essential for this exercise. (5 points)
2.	Construct a BLAST database for each of 3 the protein sets using makeblastdb (part of the blast+ suite of tools.) (4 points) Include -hash_index for efficiency, but otherwise run this program with default settings. To load blast+ on Talapas, in your sbatch script type: ml purge ; ml blast/2.2.29+
3.	Execute BLASTP using gar proteins as the query and puffer as the subject database. Specify the tabular BLAST output. Important: Be sure to run BLASTP as follows, so all of us get consistent results:
```
zcat $src/Lepisosteus_oculatus.LepOcu1.pep.all.fa.gz | \
blastp -db $src/blast_dbs/puffer -evalue 1e-5 -max_target_seqs 1 \
-num_threads 12 -outfmt '6 qseqid qlen sseqid slen pident nident \  
length mismatch qstart qend sstart send evalue bitscore' \
-out $src/garVSpuff.out.tsv

```
-	If you don’t understand what these flags/parameters are doing, check the BLAST documentation!
-	Obviously you’ll have to change the paths and filenames, but you get the idea. Note that -max_target_seqs 1 returns only the top target entry from the database, which is what we want in this case.
4.	You’ll have to perform 6 total blastp searches this way: 1. gar against puffer (example above); 2. puffer against gar; 3. gar against stickleback; 4. stickleback against gar; 5. puffer against stickleback; and 6. stickleback against puffer. Each blast search should take 50 minutes, tops. Make sure you are running blastp with the same arguments. (6 points)
-	Verify that you only have a single hit per query in your output files. If not, why is this the case? If there are multiple hits per query, make sure you write a new file in which only the first hit per query is included (hint: look into sort-ing with keys).
5.	Read the tabular blast output files (you should have 6 total) into R. Use your knowledge of R to produce a data frame of reciprocal best blast hits for each species pair (3 total). Each RBBH data frame should have two columns, corresponding to sequence IDs for the two species. Here is one possible strategy:
-	Work with one pair of blast output data frames at a time. One could create a new column in both data frames that contains an ID for each query-hit combination (hint: look into paste()). Then, using the powerful data frame manipulation function merge(), one could get the intersecting rows of the two data frames, and therefore the RBBHs. Other strategies are, of course, possible. (25 points)
6.	How many rows (RBBH orthologs) are there for each species pair? Using head(), print the first 5 rows of each of your 3 RBBH data frames and include this information with what you turn in. (6 points)
7.	Compare the number of gar-puffer, gar-stickleback, and puffer-stickleback orthologs, and briefly explain why (biologically) they are different. (4 points)
# Informative plots
8.	Plot the distributions of percent sequence identity for the stickleback-puffer and stickleback-gar RBBH orthologs side-by-side in a boxplot. Do the distributions appear to differ? Why might this be? (5 points)
9.	Using Unix commands, extract protein IDs, chromosome IDs, and genomic coordinates from the headers of the gar and stickleback XXX.pep.all.fa.gz files, and combine it with your gar-stickleback RBBH information in R to produce an orthology dot plot for gar chromosome 11.
-	The dot plot x-axis should reflect the starting basepair position (on chromosome 11) for gar orthologs. The y-axis should reflect the stickleback chromosome number (there are 21 stickleback chromosomes) of the stickleback orthologs.
-	It may be easiest for plotting to subset only those orthologs that correspond to a stickleback chromosome (they’re named “groupI”, “groupII”, etc.).
-	The goal here is to visualize where (in the stickleback genome) gar chromosome 11 genes have orthologous partners. Again, you’ll have to aggregate the appropriate information into a data frame in R before plotting. Again, the merge() function should come in handy. (15 points)
10.	How are the gar chromosome 11 genes distributed across the stickleback genome? Why is this expected/unexpected?
-	If you are stumped, try consulting the gar genome paper (Braasch et al. 2016 Nat. Genet.) (5 points)



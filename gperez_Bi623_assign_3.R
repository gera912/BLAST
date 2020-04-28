# Program Header
# Course: Bi623
# Name:   Gera Perez
# Description: R script that will identify (coding) gene orthologs between genome pairs (among 
# 3 ray-finned fish species) using a reciprocal best BLAST hit approach.
#
# Assignment 3: Genome-wide Ortholog Identification and Visualization
#
#
# Development Environment: RStudio Version 1.2.1335
# Version: R version 3.6.0
# Date:  08/27/2019
#################################################

#1
# Using Bash, in talapas, set a directory to store files and did the following command to copy files. 

#scp /projects/bgmp/shared/Bi623/assign2/*pep.all.fa.gz .


#2-3 
# Constructed a BLAST database for each of 3 the protein sets using makeblastdb and did 6 total 
# blastp searches using bash in talpas. Batch scripts can be found in Bi623_Assign_batch_scripts.txt file.



#4 
# Did not have a single hit per query in output files. This can be due to BLAST being a local alignment algorithm.

# Writing new files in which only the first hit per query is included using bash:

#sort -u -k 1,1 gar_VS_puff.out.tsv > gar_VS_puff.out_uniq.tsv
#sort -u -k 1,1 gar_VS_stick.out.tsv > gar_VS_stick.out_uniq.tsv
#sort -u -k 1,1  puff_VS_gar.out.tsv >  puff_VS_gar.out.tsv_uniq.tsv
#sort -u -k 1,1  puff_VS_stick.out.tsv >  puff_VS_stick.out.tsv_uniq.tsv
#sort -u -k 1,1 stcik_VS_gar.out.tsv > stick_VS_gar.out.tsv_uniq.tsv
#sort -u -k 1,1 stick_VS_puff.out.tsv > stick_VS_puff.out.tsv_uniq.tsv


# Copied files to local folder on laptop to open with R-studiio.

#scp hpc:/home/gperez8/bgmp/projects/bgmp/gperez8/Bi623/Ex3/*_uniq.tsv .



#5

##### gar_VS_puff_merge

# sets the directory to save files from this Rscript
setwd("/Users/gerardoperez/Documents/Bio623/Ex3")

# Opens and reads tsv file and stores it to a variable as a dataframe.
gar_VS_puff_txt <- read.table("gar_VS_puff.out_uniq.tsv", sep="\t")

# adds column names to dataframe
colnames(gar_VS_puff_txt)<- c("qseqid", "qlen", "sseqid", "slen", "pident", "nident", "length", "mismatch", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extracts columnn 1 and 3 and paste them together in a variable.
gar_VS_puff_col<-paste(gar_VS_puff_txt[,1],gar_VS_puff_txt[,3] )

# adds columns to dataframe and stores the result to a new variable.
gar_VS_puff_df<-cbind(gar_VS_puff_txt, gar_VS_puff_col)

# Opens and reads tsv file and stores it to a variable as a dataframe.
puff_VS_gar_txt <- read.table("puff_VS_gar.out.tsv_uniq.tsv", sep="\t")

# adds column names to dataframe
colnames(puff_VS_gar_txt)<- c("qseqid", "qlen", "sseqid", "slen", "pident", "nident", "length", "mismatch", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extracts columnn 3 and 1 and paste them together in a variable.
gar_VS_puff_col<-paste(puff_VS_gar_txt[,3], puff_VS_gar_txt[,1] )

# adds columns to dataframe and stores the result to a new variable.
puff_VS_gar_col_df<-cbind(puff_VS_gar_txt, gar_VS_puff_col)

# merges two dataframes by similiar column name and stores the result to a new variable.
gar_VS_puff_merge<-merge(gar_VS_puff_df, puff_VS_gar_col_df, by="gar_VS_puff_col")

#### stick_VS_puff_merge

# Opens and read tsv file and stores it to a variable as a dataframe.
puff_VS_stick_txt <- read.table("puff_VS_stick.out.tsv_uniq.tsv", sep="\t")

# adds column names to dataframe
colnames(puff_VS_stick_txt)<- c("qseqid", "qlen", "sseqid", "slen", "pident", "nident", "length", "mismatch", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extracts columnn 1 and 3 and paste them together in a variable.
puff_VS_stick_col<-paste(puff_VS_stick_txt[,1], puff_VS_stick_txt[,3])

# adds columns to dataframe and stores the result to a new variable.
puff_VS_stick_df<-cbind(puff_VS_stick_txt, puff_VS_stick_col)

# Opens and reads tsv file and stores it to a variable as a dataframe.
stick_VS_puff_txt<- read.table("stick_VS_puff.out.tsv_uniq.tsv", sep="\t")

# adds column names to dataframe
colnames(stick_VS_puff_txt)<- c("qseqid", "qlen", "sseqid", "slen", "pident", "nident", "length", "mismatch", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extracts columnn 3 and 1 and paste them together in a variable.
puff_VS_stick_col<-paste(stick_VS_puff_txt[,3], stick_VS_puff_text[,1])

# adds columns to dataframe and stores the result to a new variable.
stick_VS_puff_df<-cbind(stick_VS_puff_txt, puff_VS_stick_col)

# merges two dataframes by similiar column name and stores the result to a new variable.
stick_VS_puff_merge<-merge(puff_VS_stick_df, stick_VS_puff_df, by="puff_VS_stick_col")

#### Stick_VS_gar_merge

# Opens and reads tsv file and stores it to a variable as a dataframe.
gar_VS_stick_txt <- read.table("gar_VS_stick.out_uniq.tsv", sep="\t")

# adds column names to dataframe
colnames(gar_VS_stick_txt)<- c("qseqid", "qlen", "sseqid", "slen", "pident", "nident", "length", "mismatch", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extracts columnn 1 and 3 and paste them together in a variable.
gar_VS_stick_col <- paste(gar_VS_stick_txt[,1], gar_VS_stick_txt[,3] )

# adds columns to dataframe and stores the result to a new variable.
gar_VS_stick_df<-cbind(gar_VS_stick_txt, gar_VS_stick_col)

# Opens and reads tsv file and stores it to a variable as a dataframe.
Stick_VS_gar_txt <- read.table("stick_VS_gar.out.tsv_uniq.tsv", sep="\t")

# adds column names to dataframe
colnames(Stick_VS_gar_txt)<- c("qseqid", "qlen", "sseqid", "slen", "pident", "nident", "length", "mismatch", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extracts columnn 3 and 1 and paste them together in a variable.
gar_VS_stick_col<-paste(Stick_VS_gar_txt[,3], Stick_VS_gar_txt[,1])

# adds columns to dataframe and stores the result to a new variable.
Stick_VS_gar_df<- cbind(Stick_VS_gar_txt, gar_VS_stick_col )

# merges two dataframes by similiar column name and stores the result to a new variable.
Stick_VS_gar_merge<- merge(gar_VS_stick_df, Stick_VS_gar_df, by="gar_VS_stick_col" )


#6

# computes the number of rows
length((rownames(gar_VS_puff_merge)))
#There are 12288 rows for each species pair in gar_VS_puff_merge.

# prints the first 5 line of the data frame
head(gar_VS_puff_merge, n=5)
"
                            gar_VS_puff_col             qseqid.x qlen.x             sseqid.x slen.x pident.x
1 ENSLOCP00000000001.1 ENSTNIP00000006678.1 ENSLOCP00000000001.1    368 ENSTNIP00000006678.1    369    76.92
2 ENSLOCP00000000005.1 ENSTNIP00000008370.1 ENSLOCP00000000005.1    618 ENSTNIP00000008370.1    852    82.90
3 ENSLOCP00000000006.1 ENSTNIP00000012440.1 ENSLOCP00000000006.1    313 ENSTNIP00000012440.1    338    73.80
4 ENSLOCP00000000008.1 ENSTNIP00000010908.1 ENSLOCP00000000008.1    509 ENSTNIP00000010908.1    511    83.75
5 ENSLOCP00000000010.1 ENSTNIP00000010906.1 ENSLOCP00000000010.1    149 ENSTNIP00000010906.1    199    51.01
  nident.x length.x mismatch.x qstart.x qend.x sstart.x send.x evalue.x bitscore.x             qseqid.y qlen.y
1      280      364         82        1    364        1    362    0e+00        596 ENSTNIP00000006678.1    369
2      320      386         51        1    380        1    377    0e+00        556 ENSTNIP00000008370.1    852
3      231      313         82        1    313       26    338   4e-180        505 ENSTNIP00000012440.1    338
4      433      517         60        1    504        1    506    0e+00        824 ENSTNIP00000010908.1    511
5       76      149         62        3    149       59    198    4e-38        131 ENSTNIP00000010906.1    199
              sseqid.y slen.y pident.y nident.y length.y mismatch.y qstart.y qend.y sstart.y send.y evalue.y
1 ENSLOCP00000000001.1    368    76.92      280      364         82        1    362        1    364    0e+00
2 ENSLOCP00000000005.1    618    94.44      289      306         17        1    306        1    306    0e+00
3 ENSLOCP00000000006.1    313    73.80      231      313         82       26    338        1    313   5e-180
4 ENSLOCP00000000008.1    509    83.75      433      517         60        1    506        1    504    0e+00
5 ENSLOCP00000000010.1    149    51.01       76      149         62       59    198        3    149    6e-38
  bitscore.y
1        596
2        576
3        505
4        822
5        131
"


# computes the number of rows
length(rownames(stick_VS_puff_merge))
#There are 14628 rows for each species pair in stick_VS_puff_merge.

# prints the first 5 line of the data frame
head(stick_VS_puff_merge, n=5)
"
                          puff_VS_stick_col             qseqid.x qlen.x             sseqid.x slen.x pident.x
1 ENSTNIP00000000004.1 ENSGACP00000025030.1 ENSTNIP00000000004.1    705 ENSGACP00000025030.1    708    87.50
2 ENSTNIP00000000005.1 ENSGACP00000002858.1 ENSTNIP00000000005.1    272 ENSGACP00000002858.1    282    63.97
3 ENSTNIP00000000007.1 ENSGACP00000022871.1 ENSTNIP00000000007.1   1057 ENSGACP00000022871.1   1046    85.08
4 ENSTNIP00000000008.1 ENSGACP00000014423.1 ENSTNIP00000000008.1    108 ENSGACP00000014423.1    108    84.11
5 ENSTNIP00000000009.1 ENSGACP00000025969.1 ENSTNIP00000000009.1    843 ENSGACP00000025969.1    846    87.98
  nident.x length.x mismatch.x qstart.x qend.x sstart.x send.x evalue.x bitscore.x             qseqid.y qlen.y
1      623      712         78        1    705        1    708    0e+00       1256 ENSGACP00000025030.1    708
2      174      272         97        1    272        1    271   1e-131        378 ENSGACP00000002858.1    282
3      901     1059        143        1   1057        1   1046    0e+00       1831 ENSGACP00000022871.1   1046
4       90      107         17        1    107        1    107    9e-65        194 ENSGACP00000014423.1    108
5      732      832         95       14    842       14    843    0e+00       1476 ENSGACP00000025969.1    846
              sseqid.y slen.y pident.y nident.y length.y mismatch.y qstart.y qend.y sstart.y send.y evalue.y
1 ENSTNIP00000000004.1    705    87.50      623      712         78        1    708        1    705    0e+00
2 ENSTNIP00000000005.1    272    63.97      174      272         97        1    271        1    272   1e-131
3 ENSTNIP00000000007.1   1057    85.08      901     1059        143        1   1046        1   1057    0e+00
4 ENSTNIP00000000008.1    108    84.11       90      107         17        1    107        1    107    8e-65
5 ENSTNIP00000000009.1    843    87.40      735      841        101        5    843        5    842    0e+00
  bitscore.y
1       1256
2        378
3       1831
4        194
5       1479
"

# computes the number of rows
length(rownames(Stick_VS_gar_merge))
#There are 13091 rows for each species pair in Stick_VS_gar_merge.

# prints the first 5 line of the data frame
head(Stick_VS_gar_merge, n=5)
"
                           gar_VS_stick_col             qseqid.x qlen.x             sseqid.x slen.x pident.x nident.x length.x
1 ENSLOCP00000000001.1 ENSGACP00000013354.1 ENSLOCP00000000001.1    368 ENSGACP00000013354.1    357    78.06      281      360
2 ENSLOCP00000000004.1 ENSGACP00000007080.1 ENSLOCP00000000004.1    215 ENSGACP00000007080.1    218    76.64      164      214
3 ENSLOCP00000000005.1 ENSGACP00000015206.1 ENSLOCP00000000005.1    618 ENSGACP00000015206.1    859    84.74      261      308
4 ENSLOCP00000000006.1 ENSGACP00000022751.1 ENSLOCP00000000006.1    313 ENSGACP00000022751.1    338    78.21      244      312
5 ENSLOCP00000000008.1 ENSGACP00000015192.1 ENSLOCP00000000008.1    509 ENSGACP00000015192.1    508    82.43      427      518
  mismatch.x qstart.x qend.x sstart.x send.x evalue.x bitscore.x             qseqid.y qlen.y             sseqid.y slen.y pident.y
1         75        1    359        1    357    0e+00        589 ENSGACP00000013354.1    357 ENSLOCP00000000001.1    368    78.06
2         50        1    214        5    218   2e-122        350 ENSGACP00000007080.1    218 ENSLOCP00000000004.1    215    76.64
3         40        1    306        1    303   4e-153        466 ENSGACP00000015206.1    859 ENSLOCP00000000005.1    618    84.74
4         68        1    312       26    337    0e+00        518 ENSGACP00000022751.1    338 ENSLOCP00000000006.1    313    78.21
5         66        1    504        1    507    0e+00        808 ENSGACP00000015192.1    508 ENSLOCP00000000008.1    509    82.21
  nident.y length.y mismatch.y qstart.y qend.y sstart.y send.y evalue.y bitscore.y
1      281      360         75        1    357        1    359    0e+00        589
2      164      214         50        5    218        1    214   2e-122        350
3      261      308         40        1    303        1    306   6e-165        497
4      244      312         68       26    337        1    312    0e+00        518
5      425      517         69        1    507        1    504    0e+00        801
"

#7 

# stick_VS_puff had 14628 orthologs, the most compared to Stick_VS_gar at 13091 and gar_VS_puff at 12288 orthologs.
# They are different because this indicates that stickleback and the puffer fish are closer related species from sharing
# the most orthologous which indicates having the most common ancestor. This can be due to the stickleback and the puffer fish being
# teleost fishes which went through the teleost genome duplication (TGD). Where the spotted gar diverged from the teleost fishes before
# the TGD and having less shared orthologs to stickleback and the puffer fish.

#8

# Extracts column 6, which is the sseqid, from the merged dataframe ans stores it to a new variable.
stick_VS_puff_pident<-stick_VS_puff_merge[,6]
Stick_VS_gar_pident<-Stick_VS_gar_merge[,6]

# opens pdf file to write
pdf("stickleback-puffer__stickleback-gar_RBBH_orthologs.pdf")

# creats box plot
boxplot(stick_VS_puff_pident, Stick_VS_gar_pident, main="stickleback-puffer and stickleback-gar RBBH orthologs" ,names=c("stick_VS_puff_pident","Stick_VS_gar_pident"), ylab="percent sequence identity")
dev.off()

# The distributions appears to be different because the puffer fish has a higher sequence homology to the stickleback
# than the gar to stickleback. This can be be explained by the stickleback and the puffer fish being
# teleost fishes which went through the teleost genome duplication (TGD). Where the spotted gar diverged from the teleost fishes before
# the TGD which shows less sequence homology to the stickleback when compared to stickleback to the puffer fish.


#9

# Using Bash, in talapas, extracted protein IDs, chromosome IDs, and genomic coordinates from the headers of the Lepisosteus_oculatus.LepOcu1.pep.all.fa.gz file into a new file.

#zcat Lepisosteus_oculatus.LepOcu1.pep.all.fa.gz | grep "^>" | grep "LG11:"| sed 's/>//g'  | sed 's/chromosome:LepOcu1:/\t/g'| awk '{print $1,$3}'| grep -v "scaffold" | sort -u -k 1,1| sed 's/ /\t/g'| sed 's/:/\t/g'| awk '{print $1,$2,$3}'| sed 's/ /\t/g'> gar_LG11.tsv



# Using Bash, in talapas, extracted protein IDs, chromosome IDs, and genomic coordinates from the headers of the Gasterosteus_aculeatus.BROADS1.pep.all.fa.gz file into a new file.

#zcat Gasterosteus_aculeatus.BROADS1.pep.all.fa.gz| grep "^>" | sed 's/>//g'  | sed 's/chromosome:LepOcu1:/\t/g'| awk '{print $1,$3}'| grep -v "scaffold" | sort -u -k 1,1| sed 's/ /\t/g'| sed 's/group:BROADS1//g'| sed 's/:/\t/g' | awk '{print $1,$2,$3}'| sed 's/ /\t/g'| grep -v "MT" >stick_groups.tsv




# Opens and reads tsv file and stores it to a variable as a dataframe.
gar_header<- read.table("gar_LG11.tsv", sep="\t")

# adds column names to dataframe
colnames(gar_header)<-c("qseqid.x", "chromosome", "region")

# Opens and reads tsv file and stores it to a variable as a dataframe.
stick_header<- read.table("stick_groups.tsv", sep="\t")

# adds column names to dataframe
colnames(stick_header)<-c("qseqid.y", "chromosome", "region")

# merges two dataframes by similiar column name and stores the result to a new variable.
Stick_VS_gar_combine_ch11<-merge(Stick_VS_gar_merge, gar_header, by ="qseqid.x")


# merges two dataframes by similiar column name and stores the result to a new variable.
Stick_VS_gar_combine<-merge(Stick_VS_gar_combine_ch11, stick_header, by ="qseqid.y")

# opens pdf file to write
pdf("gar_chromosome_11_vs_Stickleback's_chromosomes.pdf")

# creates a scatter plot with certain parameters
plot(Stick_VS_gar_combine$region.x/1000000, Stick_VS_gar_combine$chromosome.y ,yaxt="n", ylab="Stickleback chromosome numbers", xlab="Starting basepair position on chromosome 11 \n gar orthologs (MegaBases)", main = "Spotted gar chromosome 11 vs Stickleback's chromosomes ")

# customizes the y-axis scale
axis(side=2, at=1:21, labels=levels(Stick_VS_gar_combine$chromosome.y), las = 1, cex.axis=0.5)
dev.off()

#10 
# From the scatterplot, you can observe the genome duplication on the stickleback. Spotted gar chromosome 11 is distributed 
# mostly in stickleback's chromosome 10 and chromosome 20. This supports that the spotted gar diverged while stickleback did 
# a genome duplication which lead to two orthologous. Also, in both stickleback's chromosome 10 and 20 you can observe gaps that could 
# possibly be a result of chromosome translocation or assembly error. 

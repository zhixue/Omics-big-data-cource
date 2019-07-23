<h2 style="text-align:center">Omics Big data

Final exam 

</h2>

# Q1 Determining sequencing platform and sequencing depth.
The sequencing error of platform A is about 15/100, and the read length can be as long as
100kb. The sequencing error of platform B is about 6/1000, and the read length can be
about 150 bases. We want to use one or both of these two sequencing platforms to
identify a mutation that happens about once in a million replicates. Which sequencing
platform will you choose and what is the minimum sequencing depth we should do in
order to determine with high confidence whether this mutation has happened or not? You
can give your own criterion about the “high confidence”, and you need to show why this
sequencing strategy or sequencing depth is enough. 
### Answer 

Define: 
```math
x_i = 1(mutation), x_i = 0(no~ mutation)   

P(x_i=1) = 10^{-6}

P(x_i=0) + P(x_i=1) =1
```
and
```math
y_i = 1(sequencing~error), y_i = 0(sequencing~no~error)
```

Using platform A:
```math
P(y_{iA} = 1) = 0.15

P(y_{iA} = 1) + P(y_{iA} = 0) = 1
```
Using platform B:
```math
P(y_{iB} = 1) = 0.006

P(y_{iB} = 1) + P(y_{iB} = 0) = 1

```

Suppose:
1. all reads from platform A and B can be mapped to the correct genome location `$i$` 
2. the error of base a,c,g,t probablity is equal 
3. `$x_i$` and `$y_i$` are i.i.d 


#### If use platform A:
```math
P(x_i = 1, y_{iA} = 0) = 10^{-6} * 0.85 = 8.5 * 10^{-7}

P(x_i = 0, y_{iA} = 1) = (1-10^{-6}) * 0.15 = 0.14999985
```
if we know that the location `$i$` is a mutaton(1/million), which is base **g->a**, and we will observe
```math
P(base_{i}=a) = P(x_i = 1, y_{iA} = 0) + P(x_i = 0, y_{iA} = 1)/3 = 8.5 * 10^{-7} + 0.14999985 /3 = 0.0500008

P(base_{i}=g) = P(x_i = 0, y_{iA} = 0) + P(x_i = 1, y_{iA} = 0)/3 = 0.85 * (1-10^{-7}) + 0.15 * 10^{-7} /3= 0.8499992

P(base_{i}=c) =  P(base_{i}=t) = 0.05


```
if we know that the location `$i$` is base **g**, 
**not** a mutaton, and we will observe
```math
P(base_{i}=g) = 0.85

P(base_{i}=a) = P(base_{i}=c) = P(base_{i}=t) = 0.05
``` 
So the true orignal base(such as base g in this example) fraction in the mutation location is lower about **`$8*10^{-7}$`** than in the no mutation location, and the mutated base(such as base a in this example) fraction is higher about  **`$8*10^{-7}$`**.  
In other words, mutated base fraction will be higher about **`$8*10^{-7}$`** than normal error base. 

We define **'high confidence'** that if an error base depth number is **more than** other error bases at least`$\delta(\delta = 5)$` in the location `$i$`, it will be a mutation. 
We need an expect sequencing depth of  **n**.
```math
E(\delta,n) = 8*10^{-7} n >=\delta=5  

n_{min} = 6.25*10^{6}
```
#### If use platform B:  
The mutated base fraction will be higher about **`$10^{-6}*0.994+(1-10^{-6})*(1-0.994)/3 - (1-0.994)/3 = 9.92*10^{-7}$`** than normal error base. We need an expect sequencing depth of  **n**.
```math
E(\delta,n) = 9.92*10^{-7}n >=\delta=5

n_{min} = 5.0403225806452*10^{6}
```
----------
# Q2 Motif finding.
Please find the top 20 most frequent k-mers (DNA fragments with length k bases) in the
human genome (version GRCh38 or hg38) for each k (=3, …, 20). Compared to the
background frequency (each base is treated as independent identical distribution), which
k-mers have the top 20 enrichment scores (observed_frequency/expected_frequency) for
each k (k=3, 4, 5, …, 20)? The human genome can be downloaded from genome.ucsc.edu
or from NCBI website or you can access it (genome.fa) under directory
/share/data/reference/hg38/ in our teaching server. 
### Answer
* top 20 most frequent k-mers (complete table in [/share/home/omics2019/018080910011/final/q2/total_k_top20.txt]())


K-mer seq | K | Observed Frequeny
---|---|---
 TTT| 3| 117589724
...|...|...
GCCTGTAATCCCAGC |13 | 234964
...|...|...

* expected_frequency = (genome size-k+1)/4^k, 
 genome size = 3*10^9
* the top 20 enrichment scores (observed_frequency/expected_frequency) 

K-mer seq |  Observed Frequeny | Enrichment Scores
---|---|---
TTTTTTTTTTTTTTTTTTTT   | 471174|  172687098.328928
AAAAAAAAAAAAAAAAAAAA|    466351  |170919450.124145
GTGTGTGTGTGTGTGTGTGT |   253886|  93050203.632497
ACACACACACACACACACAC|    251092 | 92026191.796676
TGTGTGTGTGTGTGTGTGTG |   249255 | 91352924.172337
CACACACACACACACACACA|    246547  |90360431.670046
CTCCCAAAGTGCTGGGATTA |   178385|  65378794.321006
TAATCCCAGCACTTTGGGAG|    178308 | 65350573.522381
AATCCCAGCACTTTGGGAGG |   175237 | 64225040.112286
CCTCCCAAAGTGCTGGGATT|    175044  |64154304.863785
GCCTCCCAAAGTGCTGGGAT   | 171658  |62913322.732042
ATCCCAGCACTTTGGGAGGC  |  171402  |62819497.739211
TCCCAAAGTGCTGGGATTAC |   169097  |61974706.299853
GTAATCCCAGCACTTTGGGA|    169065  |61962978.175749
CCCAAAGTGCTGGGATTACA    |168661  |61814910.608937
TGTAATCCCAGCACTTTGGG    |168443 | 61735012.763479
CCAAAGTGCTGGGATTACAG   | 167385  |61347251.660295
CTGTAATCCCAGCACTTTGG  |  167164  |61266254.303202
CAAAGTGCTGGGATTACAGG |   164608  |60329470.390404
CCTGTAATCCCAGCACTTTG|    164487  |60285123.421136


### Code
```shell
#shell
filename='/share/data/reference/hg38/genome.fa'
for k in `seq 3 20`
do
    echo $k
    jellyfish count -m $k -s 3G -o human_k${k} $filename 
    jellyfish dump -c -t -L 10000 human_k${k}_0 | sort -rnk 2 | head -20 > k${k}_top20.txt 
done

cat k*_top20.txt > total_k_top20.txt
awk '{b[NR]=$0; for(i=2;i<=NF;i++)a[NR]=$2/(3*10^9-length($1)+1)*(4^length($1));}END{for(i=1;i<=NR;i++) printf(b[i]"\t%f\n",a[i])}' total_k_top20.txt | sort -rnk 3 | head -20   
```

----------------
# Q3 Evaluation of the novel genes found in the rice pan-genome.
Please create a pipeline to evaluate the quality of the novel genes in the rice pan-genome
(but not included in the reference rice genome). You need to answer what is the
percentage of the novel proteins (with names started with “Un_maker_”) in the pangenome protein sequences that can be aligned to the reference rice genome with at least
30%, 50%, or 90% of the protein coding region length? For alignment tools, such as Blast
or any other alignment tools that are comfortable for you to use, please explain why you
choose this tool. Or you can develop your own program to do the alignment.
The rice pan-genome and its annotations are available at
http://cgm.sjtu.edu.cn/3kricedb/data_download.php. The reference rice genome is
included in the pan-genome, i.e., chromosome 1-12. The protein sequences of the rice
pan-genome are also provided in the above link.
### Answer
Sequences from reference and from novel genes are splited, blastp is used in this question (E-value threshold 10). The subject is the reference cds sequence(protein sequence), queries are the novel proteins sequences. 
Novel proteins  | Number | Percent
---|---|---
All | 15362 | 100%
Coverage>30% | 12780 | 83%
Coverage>50% | 9449 | 62%
Coverage>90% |  3855 | 25%

 

#### Pipeline
```
graph TB
subgraph Evaluation of novel proteins
A(Novel proteins not from reference  )
B(Mapping to proteins from reference - blastp)
C(Filtering by subject coverage 30%, 50%, 90% - awk)
D(Filtered novel proteins id and numbers)

end

A-->B
B-->C
C-->D
```

#### Code
```shell
#shell
wget http://cgm.sjtu.edu.cn/3kricedb/data/pep.fa
grep -n '>Un' pep.fa |head -1
awk 'NR>=249470{gsub(/,/," ");print}' pep.fa > pan_pep.fa
head -249469 pep.fa > ref_pep.fa

#create database
#makeblastdb -in ref_pep.fa -dbtype prot -parse_seqids

# mapping
thread=1
blastp -db ref_pep.fa -query pan_pep.fa -num_threads $thread -max_target_seqs 1 -out pep.blastout -outfmt "6 qacc sacc evalue pident length qstart qend sstart send"

# filtering
grep '>' pan_pep.fa | wc -l # 15362
awk '{print $1}' pep.blastout |sort|uniq|wc -l # 14538
awk '{if ($10>90) print $1}' pep.blastout |sort|uniq|wc -l # 3855
awk '{if ($10>50) print $1}' pep.blastout |sort|uniq|wc -l # 9449
awk '{if ($10>30) print $1}' pep.blastout |sort|uniq|wc -l # 12780
```
* notes: subject coverage(not query coverage)
----------------
# Q4 Planning a big omics research project.
Dr. Z is planning a project to investigate the impact of the host genome and gut
microbiomes on a therapeutic food intervention. A cohort of 30,000 patients is going to be
collected. At the first stage of the project, 100 patients will have their genomes sequenced
(with 30x coverage), together with RNA-seq data (~30GB each) before and after a 3-month
therapeutic food intervention. All these 100 patients will also have their gut microbiomes
sequenced (at least 30GB per sample) before and after the 3 months’ therapeutic food
intervention. Including these 100 patients, a total number of 1,000 patients will have 16S
rRNA sequenced (at least 30MB per sample) before and after the 3 months’ therapeutic
food intervention. Please give a strategy to cover the sequencing and data analysis stages.  
You need to show the strategy (with diagrams) and give your estimation about the time and cost for the strategy. For example, you can pick the different sequencing platforms for
genome, transcriptome and metagenome sequencing. For the 16S rRNA sequencing, 100
samples can be merged into one run of illumine sequencing. For the cost estimation,
please give a reasonable estimation. For example, you can pick cloud computing or the Pi
supercomputer from SJTU to do the data analysis. You can find current market price from
computing resource providers such as cn.aliyun.com or www.huaweicloud.com.  
If Dr. Z wants to do the sequencing in a scale 10 times bigger in the second stage of the
project above, what strategy will you propose and what will be the cost in terms of time
and money? 
### Answer

#### stage 1, 100 patients, 2 times
##### Sequencing and cost
* sequencing type   
(**a**) DNA-seq (100x2, 30x)  
(**b**) RNA-seq (100x2, 30G)  
(**c**) metagenome sequencing (100x2,30G)  
(**d**) 16S rRNA-seq (1000x2, 30MB, near 20 samples DNA-seq or RNA-seq)  

* computing in SJTU PI HPC
* * ￥0.1(or 0.05￥) = 1 point (1p) 
* * cpu (16p/h), cpu128 (24p/h), fat (32p/h)  
Supposing this is a project 1 years finished, using average 32 cpu 24 hour/every day    
12 * 365 * 32 * 16 * 0.1 = ￥224 256

* storing in SJTU PI HPC
* * upper than 2T (166p/TB/month)  
raw data: 0.03 * (100 * 2+100 * 2+100 * 2+20) = 18.6TB
temp data/results: 18.6 * 5 = 93TB  
Supposing 1 years finished  
166 * 0.1 * (18.6+93) * 12 = ￥22 230.72

* sequencing platforms
*  * Illumina NextSeq/MiSeq/HiSeq (supposing ￥6000/sample)  
6000 * (100 * 2+100 * 2+100 * 2+20) = ￥3 720 000

Total, nearly ￥ 4 000 000
##### Analysis pipeline

```
graph TB
A[Raw sequences - a,b,c,d] 
B(Quality control)
subgraph human
C1(Mapping - a,b)
C2(De novo assembly - a)
D1(SV, Indel, SNV, CNV ,Gene fusion Dectection - a,b)
D2(Different gene expression - b)
E1(GWAS -a,b)
E2(GSEA, Pathway analysis - b)
end

subgraph microbiome
C3(Mapping - c,d)
C4(De novo assembly - c)
D3(Taxonomic analysis - c,d)
E3(ORF/Gene detection - c)
F1(Protein function prediction - c)
end


Z[Merged results / Experiment / New discovery - impact of the host genome and gut microbiomes on a therapeutic food intervention]

A-->B
B-->C1
B-->C2
B-->C3
B-->C4
C1-->D1
C2-->D1
C1-->D2
C3-->D3
C4-->D3
D1-->E1
D2-->E2
C4-->E3
E3-->F1

E1-->Z
E2-->Z
D3-->Z
F1-->Z

```

#### stage 2, 1000 patients, 2 times
using sequencing platform, computing source, analysis in stage 1
* time: supposing this a project 2 years finished
* money:  
* * sequencing: 6000 * (1000 * 2+1000 * 2+1000 * 2+200) = ￥37 200 000
* * storing: 166 * 0.1 * (18.6+93) * 10 * 12 * 2 = ￥444 614
* * computing: Supposing this is a project 1 years finished, using average 64 cpu 24 hour/every day   
2 * 12 * 365 * 64 * 16 * 0.1 = ￥897 024

Total, nearly ￥38 500 000

* notes: 

1. analysis step - data source(a=DNA-seq,b=RNA-seq,c=metagenome,d=16sRNA-seq)
2. computing price - https://hpc.sjtu.edu.cn/regulation_20181008.pdf

------------------
# Q5 Estimating the detectability of peptides .
Proteomics search engines identify peptides by comparing the theoretical and observed
peaks in the MS/MS spectrum, and measure the fitness between theoretical and observed
distribution of m/z of tryptic peptides.
Using human proteome sequences, generate theoretical tryptically-digested peptides
database with two misscleavages at most. Find out the percentage of detectable peptides,
and evaluate the relationship between the detectability in MS/MS and the expression
level, misscleavages, charge state and lengths of the candidate peptides. Graphs are
encouraged to show the results.
### Answer 
We define detected peptides with two misscleavages at most. Nearly 98% peptides 
are with two misscleavages at most. In this part, we use rTANDEM(a package of R) to deal with the xml ouput of X!TANDEM from public sample files  'TCGA-A6-3807-01A-22_Proteome_VU_20121019_OUT' and 'TCGA-A6-3808-01A-22_Proteome_VU_20121205_OUT' in our server.

The relationship between the detectability in MS/MS and the expression
level, misscleavages, charge state and lengths of the candidate peptides:  
* File 1  
![VvmlYF.png](https://s2.ax1x.com/2019/06/20/VvmlYF.png)
* File 2  
![VvnUjs.png](https://s2.ax1x.com/2019/06/20/VvnUjs.png)

Detectability in MS/MS will be higher if the length of the candidate peptide is higher.
 
#### Code
```R
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("rTANDEM")
library('rTANDEM')
library('corrplot')

read_xml = function(xmlfile){
    result.R=GetResultsFromXML(xmlfile)

    # peptides misscleavages<=2
    detected.idx = which(result.R@peptides[["missed.cleavages"]]<=2)
    detected_num = length(result.R@peptides[["missed.cleavages"]][detected.idx])
    uniq_detected_num = length(unique(result.R@peptides[["spectrum.id"]][detected.idx]))
    total_num = length(result.R@peptides[["spectrum.id"]])
    uniq_total_num = length(unique(result.R@peptides[["spectrum.id"]]))
    print(detected_num/total_num)
    print(uniq_detected_num/uniq_total_num)


    result = data.frame(matrix(NA,1,5))
    for (i in detected.idx){
       # 1 detected score
       # 2 peptides.missed.cleavages: misscleavages
       # 3 peptides.spectrum.z:  charge
       # 4 a@peptides[["end.position"]]-a@peptides[["start.position"]]: length
       # 5 sum(a@peptides[["spectrum.id"]]=="6500"): expression 
    
       result[i,1] = result.R@peptides[["tandem.score"]][i]
       result[i,2] = result.R@peptides[["missed.cleavages"]][i]
       result[i,3] = result.R@peptides[["spectrum.z"]][i]
       result[i,4] = result.R@peptides[["end.position"]][i]-result.R@peptides[["start.position"]][i]
       result[i,5] = sum(result.R@peptides[["spectrum.id"]]==result.R@peptides[["spectrum.id"]][i])
    }
    colnames(result) = c('detected score','misscleavages','charge','length','expression')
    M = cor(na.omit(result), method='spearman')
    corrplot(M, method = "color", type = "upper")
}

xmlfile='output.2017_06_09_10_56_26.t.xml'
read_xml(xmlfile)
read_xml('output.2017_06_09_11_12_31.t.xml')
```


* notes: about mgf file http://www.matrixscience.com/help/data_file_help.html#GEN
-----------------
# Q6 Correlation between gene expression and protein expression.
Of all proteome spectrum or RNA-Seq reads from human samples, only part of them can
be mapped into human sequence. Estimate the percentage of unmatched reads or
spectrums in human RNA-Seq or proteomics sample. How about the average correlation
between gene expression and protein expression? Which genes show the higher
correlation, which not? Is there any rules?
### Answer 
There is no significant correlation between gene expression and protein expression (r = 0.1~0.9, average r = 0.2~0.3).   
Genes **encoding intermediary metabolism functions** showed high mRNA–protein correlations, whereas genes **involved in oxidative phosphorylation, RNA splicing and ribosome components** showed low or negative correlations.  
Many factors such as **protein degradation, mRNA degradation, polyadenylation, codon preference, translation rates, alternative splicing, translation lag** may have a marked impact on it.

![VvkAzt.png](https://s2.ax1x.com/2019/06/20/VvkAzt.png)

![VvAQ0O.jpg](https://s2.ax1x.com/2019/06/20/VvAQ0O.jpg)

![VvE8bT.jpg](https://s2.ax1x.com/2019/06/20/VvE8bT.jpg)
* Source: Gygi et al. Mol Cell Biol, 1999; Zhang et al. Nature, 2014







---------------
# Q7 Where are the unexplained spectrums in shotgun proteomics from?
Many spectrums from real sample still cannot be explained by the known human protein.
Please try to propose a possible reason, and test your hypothesis.   

Some data for questions 5-7 are available at our teaching server under directory
/share/data/proteomics/final-examination/, which include mgf and xml (searching results by
X!tandem) files from two human samples and a fasta file of human proteome sequences. 

### Answer
* Possible reasons to explain unmapped spectrums:  
1. some proteins are not from microbiome instead of human genome(the microbial products have been reported in the tissue far away) 
2. some proteins expression in the special tissue(e.g. brain, malignant solid tumors) or situation(e.g. no oxygen, huge physical attack, other disease)
3. althouth the annotation of transcrptions in human is better than in other species, all possibilities of alternative splicing and noncoding RNA effect(have been reported to code proteins) remain unknown at the level of RNA.



We try to test the first hypothesis/reason, using the strategy as follows:
* **Step 1**: download large amount of protein sequencing data from various databases such as Proteomics DB(https://www.proteomicsdb.org/) , www.peptideatlas(http://www.peptideatlas.org/repository/repository_public_Hs_Plasma2.php), Pride(https://www.ebi.ac.uk/pride/archive/)     
> **Status in Proteomics DB**  
Human Proteome  
Coverage:80%  
Proteins:15721 of 19629  
Isoforms:11353 of 86771  
Unique Peptides (Isoform):113944  
Unique Peptides (Gene):455289  
Spectra:43237800  
**Repository Proteomics DB**  
Registered Users:786  
Projects:80  
Experiments:706  
Files:22906  
Data Volume:8.85 TB  
* **Step 2**: filter the unmapped spectrums from known proteins databases   
* **Step 3**: use microbiome references(proteins manually reported form experiments and predicted by algorithm)  
* **Step 4**: map the sepectrums to microbiome references

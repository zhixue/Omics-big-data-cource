<h2 style="text-align:center">Omics Big data

Project 1  

</h2>

### Which genome annotation set is the best for the rice reference genome?
> The goal of this project is to evaluate the quality of the three genome annotation datasets for the rice genome.  
You are expected to evaluate two rice genome annotation data by these RNA-seq data. Please pick at least 5 RNA-seq data, and align these RNA-seq data to the rice genome and check how many of the genes/transcripts can be validated by at least 2 RNA-seq datasets. You can use any sequencing alignment (or mapping ) tools to align these RNA-seq data back to the rice genome. You are expected to report the description of steps in your analysis pipeline, the reasons why you choose these tools for each step, together with the following results.
> 1. The total number of RNA-seq reads for your analysis;
> 2. The number of reads that can be aligned to the rice genome (with  the criteria you use) for
each RNA-seq dataset;
> 3. The number of reads that can’t be aligned to the rice genome;
> 4. For all genome annotation sets,  
> a) How many of the aligned reads locate inside an exon, i.e., not include exon junctions;  
> b) How many of the aligned reads include exon junction;  
> c) How many genes are validated with at least two exons with their junctions supported
with at least k reads, where k=1, 5, 10 and 50;  
> d) the percentages of exons in the rice genome that can be covered by at least k reads,
where k=1, 5, 10, and 50; 
-----------------
#### 1. Pipeline
Steps in my analysis pipeline are as following
```
graph TB
step0["raw data(RNAseq,fastq)"]

subgraph  
step1("quaility control(fastqc)")
step2("mapping(hisat2)")
step3("filtering aligned reads and sort(samtools)")
step4("running scripts to compare different gff files")
end 

step4o["comparation detail for annotations"]

step2o["aligning summary(output by hisat2)"]


step0-->step1
step1-->step2
step2-->step2o
step2-->step3
step3-->step4
step4-->step4o
```
#### 2. Method description
* 2.1 **quaility control(fastqc)**  

In this step, RNA-seq raw data(fastq format) 's quaility are controled by ***fastqc***[(FastQC v0.10.1)](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), using the following command: 
```shell
fastqc -o fastqc_out/${sample} -t 4 ${fastqgz}
```
The samples with **short read length** and **low quaility** are ***not*** taken into account in next steps. All the results of fastqc reports are in `/share/home/omics2019/018080910011/rice_project/fastqc_out`. Here show an example of a sample with bad quaility.(DRR001965)
![EFmx2T.png](https://s2.ax1x.com/2019/04/21/EFmx2T.png)
* 2.2 **mapping**   

***Hisat2***[(Hisat2 version 2.1.0)](http://ccb.jhu.edu/software/hisat2/index.shtml) is the ideal  alignment program for mapping next-generation sequencing reads based on BWT and GFM. Before mapping, reference index should be built. The reference(fasta format) index is built with command 
```shell
# MSU and RAP-DP
hisat2-build /share/home/ccwei/courses/2019/omics/proj1/reference/IRGSP-1.0_genome.fasta IRGSP1
# BGI
hisat2-build /share/home/ccwei/courses/2019/omics/proj1/reference/IRGSP4_genome.fasta IRGSP4
```
Run the alignment program for samples, using serveral threads. The outputs are SAM format file(`${sample}.sam`) and the alignment summary(`${sample}.summary`).
```shell
hisat2 -p ${thread} -x ${ref_idx} -U ${fastqgz} -S ${sample}.sam 2>${sample}.summary
```
The alignment summary is like this(an example of singe-end sample), which includes **the number of RNA-seq reads, the number of reads that can be aligned(1 times or more), the number of reads that can’t be aligned.**
```
2997625 reads; of these:
  2997625 (100.00%) were unpaired; of these:
    594675 (19.84%) aligned 0 times
    2299688 (76.72%) aligned exactly 1 time
    103262 (3.44%) aligned >1 times
80.16% overall alignment rate
```
* 2.3 **filtering aligned reads and sort**  

Aligned reads are selected and sorted by ***samtools***[(samtools Version: 1.3.1)](samtools.sourceforge.net).
```bash
# single-end
## drop unaligned reads
samtools view -h -F 4 ${sample}.sam > ${sample}.aligned.sam
## sort reads by reference location
samtools sort -@ ${thread} ${sample}.aligned.sam >${sample}.sorted.sam
```

* 2.4 **running scripts to compare different gff files**
##### 2.4.1 storing with less file size
The information of **Query_name, Flag, Reference name, Start position, End position, Cigar** in SAM format file(part of columns) is useful in next step and extracted by `get_location.py`, stored with a smaller size of file.
```bash
python3 get_location.py ${sample}.sorted.sam ${sample}.loc
```
${sample}.loc is like this
```
#query_name     flag    reference       ref_start       ref_end cigar
SRR1005301.28822924     16      chr01   76654   76928   114M72N88M
...
```
##### 2.4.2 comparing anotations
`compare_annotation.py` is a script to output the answer from question 4.a to 4.d.   
In this step, reads with **best aligning location** are considered(Flag do not include 256).   
In question 4.a and 4.b, aligned reads are classified with **3 types(reads completely in one exon, reads over exon junctions, reads completely not in exons)**, accorrding to their mapping locations and exon locations in gff file.  Using counters to count that the number of reads satisty these 3 types.  
In question 4.c and 4.d, I count **the number of reads(either completely in exons or over exon junctions) in all the exons and the exons belong to which genes**, because there are many genes in a chromosome, and there are at least one exon(including UTR,CDS) in a gene. Next, I filter the genes which are validated with at least two exons with their junctions supported
with at least k reads as well as the percentages of exons that can be covered by at least k reads on each chromosome(where k=1, 5, 10 and 50).

```bash
python3 compare_annotation.py ${sample}.loc $gff_file '' >${sample}.all.ans
```
${sample}.all.ans is like this
```
#DRR001963.loc
#aligned reads locate inside an exon
chr01   120848
chr02   97546
...
#aligned reads include exon junction
chr01   101212
chr02   74581
...
#genes are validated with at least two exons,the percentages of exons(k=1,5,10,50)
chr01   2903    1178    695     118
chr02   2389    1083    623     108
...
chr01   0.5018446968596578      0.15768526857043255     0.07858118118699707     0.011852540452604364
chr02   0.5151441725224134      0.17719547232510643     0.09145349440963689     0.014296098861158226
...
```

#### 3. Results
75 samples have been run. All the results is in `merge_result.xlsx`.  
##### For example ,in chr2:
![EFJSte.png](https://s2.ax1x.com/2019/04/21/EFJSte.png)
the value of differnece has been produced with `mus-rapdp`
* **column 1**: `mus-rapdp`.aligned reads locate inside an exon
* **column 2**: `mus-rapdp`.aligned reads include exon junction
* **column 3~6**: `mus-rapdp`.genes are validated with at least two exons with their junctions supported with at least k reads, where k=1, 5, 10 and 50
* **column 7~10**: `mus-rapdp`.the percentages of exons that can be covered by at least k reads,where k=1, 5, 10, and 50  


 1. It seems that there are more **reads aligned completely in exons and exon junctions** in `msu` annotations, which suggests that the ranges of exon in `msu` is wider than them in `rapdp`.
 2. The numbers of `mus-rapdp`.**reads aligned completely in exons** in some samples are low because of low number in both `msu` and `rapdp` with a lower aligning rate at the step of mapping. So these samples are not good enough for us to decide which annotation is better.
 3. Fewer **genes detected in `msu` by parameter `k=1,5,10` while more genes detected in `msu` by parameter `k=50`**, which indicates that `rapdp` catches more genes than `msu` but some genes with high confidence are not catched by `rapdp`.However, `msu` catches more genes with high confidence. 
 4. Larger **percentages of exons in `rapdp`**  by parameter `k=1,5,10,50` suggest that the sizes of exons in `rapdp` are more fragmental than them in `msu`. In other words, the lengths of exons in `msu` are longer but less detailed.  
 5. There are no obvious evidence that either of two annotations has a lack of the RNA in the specific status or tissue of there samples.
 
 #### 4. Conclusion
 `rapdp` may be a little better than `msu`.
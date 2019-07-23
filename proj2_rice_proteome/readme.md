<h2 style="text-align:center">Omics Big data

Project 2  

</h2>

### Is the proteome helpful to evaluate genome annotation?
> In the project #1, you have compared two or three commonly used genome annotations for the reference rice genome using RNA-seq datasets. The goal of project #2 is to evaluate your results further using proteome datasets.  
We have collected 4 proteome datasets for rice in different experiments. You may find the data and the brief summary at /share/data/project2.  
Please pick up two shotgun proteomics data of rice at least to evaluate the rice annotation data. You are expected to report the following results:  
>1. The total number of proteins annotated in the rice reference genome.
>2. The total number of spectra mapped to sequence(FDR<=0.05) and the
unmapped spectra in your analysis.
>3. The distribution of coverage of detected peptides in a protein.
>4. How many gene products are validated with both RNA-seq and proteome data.
-------------------------------------
#### 1. The total number of protein


Anotation | Protien numbers
---|---
MSU |  66338
RAPDP | 42229
BGI | 40745


#### 2. the total number of spectra mapped to sequence(FDR<=0.05) and the unmapped spectra
2 samples(PXD001058, PXD003156) are used in this analysis.
Anotation | Spectra number(mapped) | Spectra number(unmapped) | Spectra number(total)
---| --- | --- | ---
MSU | 217888 | 215424 | 433312
RAPDP | 212944 | 220368 | 433312
BGI|  221883 | 211429 | 433312

#### 3. the distribution of coverage of detected peptides in a protein
X: coverage(%), Y: peptides count

![Vn4oJf.png](https://s2.ax1x.com/2019/05/29/Vn4oJf.png)
![Vn4TW8.png](https://s2.ax1x.com/2019/05/29/Vn4TW8.png)

![Vn5pWT.png](https://s2.ax1x.com/2019/05/29/Vn5pWT.png)
![Vn5EwR.png](https://s2.ax1x.com/2019/05/29/Vn5EwR.png)

![Vn53md.png](https://s2.ax1x.com/2019/05/29/Vn53md.png)
![Vn580A.png](https://s2.ax1x.com/2019/05/29/Vn580A.png)


#### 4. how many gene products are validated with both RNA-seq and proteome data

Anotation | Sample | Detected protein number(unique)
---|--- | ---
MSU | PXD001058 | 2334
MSU | PXD003156 | 2203
RAPDP | PXD001058 | 2300
RAPDP | PXD003156 | 2254
BGI(IRGSP4)| PXD001058 | 2317
BGI(IRGSP4)| PXD003156 | 2285

#### Conclusion 
* It seems that `msu` anotated more proteins but had no significent advantage in  discovering more proteins or peptides.    
* Although the distributions of coverage of detected peptides in a protein from different anotation are similar, `rapdp` is better in sample PXD001058 by support of more high coverage peptides.
* As is summaried in project 1 and project 2, `rapdp` may be a little better than `msu`.
#### Methods
1. Download protein sequences(fasta format)
```shell
### database version time ###
# BGI 2017.05.17
wget ftp://public.genomics.org.cn/BGI/rice/rise2/9311_glean_gene_pep.fa.gz

# RAP 2018.05.30
wget http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_protein_2018-03-29.fasta.gz

# MSU 2018.05.30
wget http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.pep
```
2. Count anotated proteins
```shell
grep '>' msu_pep.fa | wc -l 
#66338
grep '>' IRGSP-1.0_protein_2018-03-29.fasta | wc -l 
#42229
grep '>' bgi_pep.fa | wc -l 
#40745
```

3. Run Maxquast to get the information of proteins and peptides  
Load the '*.RAW' data from different samples and set the experiment id and fractions. Remember to load the fasta file of proteins produced by different anotations.

![VnICNt.md.png](https://s2.ax1x.com/2019/05/29/VnICNt.md.png)

4. Filter number of spectra  
Total number of spectra can be found in 'summary.txt', and number of mapped spectra can be found in 'proteinGroups.txt' in every anotation.
5. Draw the distribution of coverage of peptides in a protein  
The value of distrubtion of coverage can be selected from a column of 'proteinGroups.txt'.
6. Count gene products  
Proteins(gene products) are counted from a column of 'proteinGroups.txt'.
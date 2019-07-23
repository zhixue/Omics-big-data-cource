<h2 style="text-align:center">Omics Big data

Homework 7 

</h2>

> 1. Write down your comments for three problems that you picked from the ‘Big Challenges’ in the lecture and four problems given by the teacher in the class.

* Is the best hit the right choice?  
No, we think top-k should be selected because best hit may be randomly generated. But the choice of k remains a challenge.
 
* Is it reasonable to pick minimal protein list to explain peptides?  
We don’t think so. Because long sequences are more likely to be selected. We think that adding a parameter about length may reduce the bias and improve this problem.
 
* Can identify novel peptide using data of colorectal cancer?  
Yes. And If novel peptide is found only in the patient sample rather than control sample, it is likely to be related to the occurrence of the disease.

--------

> 2. Go to the database UniProt www.uniprot.org; let us know the number of human proteins (UniProtKB), as well as the amount of well-annotated human proteins (Swiss-Prot) in the newest release.  

* 171,063 human proteins(20,421 reviewed + 150,642 unreviewed) human proteins(UniProtKB)  
* 20,421 reviewed human proteins(Swiss-Prot)  
(Updated time: 2019.5.13)

------

> 3. Using Swiss-Prot sequence database, try to analyze the proteomics data human proteomics samples (5 cancer and 5 normal) using MaxQuant software (http://www.biochem.mpg.de/5111795/maxquant).  
Describe the number of peptides and protein you identified (FDR<0.01), and report and annotate the significantly expressed proteins between the case and the control group.

* 4655(raw), 4584(filtered) proteins are detected by Maxquant(FDR<0.01).  
* 141859 petides(30470 unique) are identified.  

* Significantly expressed proteins：

Proteins | 
---|
ENSP00000462357;ENSP00000419497;ENSP00000462830;ENSP00000011292;ENSP00000462257;ENSP00000419408;ENSP00000475021;ENSP00000474175 |
ENSP00000381386;ENSP00000215773;ENSP00000399768;ENSP00000215770;ENSP00000384866;ENSP00000385714;ENSP00000409440 |
ENSP00000446231;ENSP00000437351;ENSP00000229314 |
ENSP00000230588 |
ENSP00000241052 |
ENSP00000244137;ENSP00000391890;ENSP00000380226;ENSP00000476869 |
ENSP00000258201 |
ENSP00000259396 |
ENSP00000444935;ENSP00000262890;ENSP00000473065;ENSP00000471933;ENSP00000469352;ENSP00000472950;ENSP00000470753 |
ENSP00000412936;ENSP00000397738;ENSP00000392920;ENSP00000407643;ENSP00000401992;ENSP00000394376;ENSP00000386868;ENSP00000342035;ENSP00000387084;ENSP00000380868;ENSP00000380867;ENSP00000265810 |
ENSP00000270142;ENSP00000374645 |
ENSP00000422448;ENSP00000274458;ENSP00000427319 |
ENSP00000280333 |
ENSP00000290349;ENSP00000434613;ENSP00000395132;ENSP00000382143 |
ENSP00000299299 |
ENSP00000300060;ENSP00000452934;ENSP00000453405;ENSP00000453545 |
ENSP00000300900;ENSP00000465837;ENSP00000464757;ENSP00000466964 |
ENSP00000456319;ENSP00000301729;ENSP00000457900;ENSP00000456565 |
ENSP00000374665;ENSP00000305692;ENSP00000460543;ENSP00000459972;ENSP00000458306 |
ENSP00000417774;ENSP00000373192;ENSP00000419642;ENSP00000327880;ENSP00000346560;ENSP00000419874;ENSP00000417617 |
ENSP00000442094;ENSP00000435113;ENSP00000346403;ENSP00000331209;ENSP00000433723;ENSP00000434382;ENSP00000433721;ENSP00000434627;ENSP00000434291;ENSP00000434065;ENSP00000435204;ENSP00000435231;ENSP00000435729 |
ENSP00000356015;ENSP00000437850 |
ENSP00000357861;ENSP00000397261;ENSP00000413960;ENSP00000408263;ENSP00000406222;ENSP00000396209;ENSP00000412816;ENSP00000402531;ENSP00000395637;ENSP00000416206;ENSP00000400841;ENSP00000390433 |
ENSP00000358106 |
ENSP00000361845 |
ENSP00000387649;ENSP00000363021;ENSP00000363017;ENSP00000363015 |
ENSP00000363603 |
ENSP00000363641;ENSP00000363639 |
ENSP00000445116;ENSP00000419628;ENSP00000368219;ENSP00000368209;ENSP00000368222;ENSP00000436844;ENSP00000407856 |
ENSP00000414376;ENSP00000370010;ENSP00000370009;ENSP00000370007 |
ENSP00000376028;ENSP00000356986;ENSP00000356988;ENSP00000356987 |
ENSP00000378394;ENSP00000378392 |
ENSP00000381607;ENSP00000381604 |
ENSP00000436419;ENSP00000437184;ENSP00000403925 |
ENSP00000454537 |

![VSgcjI.th.jpg](https://s2.ax1x.com/2019/05/21/VSgcjI.th.jpg) ![VSg4US.th.jpg](https://s2.ax1x.com/2019/05/21/VSg4US.th.jpg) ![VSgHvn.th.jpg](https://s2.ax1x.com/2019/05/21/VSgHvn.th.jpg)

![VSgLD0.jpg](https://s2.ax1x.com/2019/05/21/VSgLD0.jpg)
--------------

> 4. The expressions of mini-peptides have been reported in human proteome studies. In order to investigate the expression of possible mini-peptides in human proteome, please perform the following steps:  
(1) Use 3-frame translation to predict the ORF in human ref genome, of which length is between 6 and 100 AA [6,100].  
(2) Remove the entries that same with any known human proteins, as well as the
redundancy.  
(3) Construct a fasta database using the putative mini-peptides, and detect their
expression in ten human proteomics samples mentioned in #3.  
(4) Summarize your finding.  
You may find the proteomics data at /share/data/proteomics/colorectal_cancer, the reference genome at /share/data/proteomics/. 


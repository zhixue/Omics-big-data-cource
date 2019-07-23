# Homework 2, week 2
1.
A sequence assembly engine? (theory)
You’re given two sequences x and y that represent overlapping sequencing reads from the same region of a chromosome. You know (from other data) that x and y overlap, that they’re reads from the same strand, and that x is to the left of y: that is, the overlap is like this:  
```
...1---i-M......  
x: -------......  
y: ...----------  
......1--j-----N  
```
so that the optimal alignment of the two reads involves x_i..M and Y_1..j. This problem is a hybrid of local and global alignment. We want the overhangs to be free, as in local alignment (the gap penalties for x_1..x_i-1 and y_j+1..y_N are 0). Within the aligned section, though, we use normal gap penalties.

Using a scoring system of +1 for a match, -3 for a mismatch and a linear gap penalty of -10 per gap symbol, you are asked to give  
1). the dynamic programming algorithm that finds an optimal alignment of this type;   
2). your analysis about the time and space complexity of your algorithm.
-----
2. 
Paring Smith/Waterman output (Perl).
You are the new person in a bioinformatics group that’s just decided that everything in their informatics pipeline will be converted from WU-BLAST searches to Smith/Waterman searches. Fortunately you notice that the program ssearch, part of Bill Pearson’s FASTA package, a robust implementation of the Smith/Waterman local alignment algorithm, is installed on our course server for you to use. (Do a man ssearch to see the man page.)
You also have a legacy Perl script that takes a WU-BLAST output file and parses it to find the name of the query seuqence, and the name, socre and P-value of the top scoring hit. The source code for the script is here (/share/home/ccwei/courses/2019/omics/hw2/blastparser.pl). An example of WU-BLAST output is here (/share/home/ccwei/courses/2019/omics/hw2/blast.out) . Make a copy of this script and the output in files called blastparser.pl and blast.out, and make the parser executable as a program (chmod +x blastparser.pl). When you run the script on the sample output file, it produces a single summary line of output as follows:
```
(Gpu)./blastparser.pl blast.out
Best hit to RU1_HUMAN is: K08D10.3, with score 378, P-value
3.2e-53.
```
Your task is to modify the script to do the same job on an ssearch output file. An example of ssearch output with the same query is here(/share/home/ccwei/courses/2019/omics/hw2/
ssearch.out) . Parse the file to get the query name, and the name, score and E-value of the top hit. Have the script print out a summary line similar to what the BLAST parser script did, showing this information.

Getting ahead...  
p.s. Want to run the searches yourself? The query sequence is the human U1A RNA binding
protein, and its sequence is here (/share/home/ccwei/courses/2018/omics/hw2/u1_human.fa) . Get the sequence and save it to a file called u1_human.fa. The database is Wormbase C.elegans 215, the complete genome of C.elegans (/share/home/ccwei/courses/2019/omics/C.elegans/ Proteome/ws_215.protein) . The command line I used to run the ssearch program was:
```
/share/home/ccwei/tools/fasta-36.2.6/bin/ssearch36 -q u1_human.fa /share/home/ccwei/courses/2019/omics/C.elegans/Proteome/ws_215.p rotein.fa > ssearch.out
the command line I used to run the BLAST search was
/share/home/ccwei/tools/wu-blast/blastp /share/home/ccwei/courses/2019/omics/C.elegans/Proteome/ws_215. protein u1_human.fa filter=seg+xnu > blast.out
```
-----


* Q1A sequence assembly engine? (theory) 
1). the dynamic programming algorithm that finds an optimal alignment of this type;
```
GetOverlap(x,y)
	// init
	MisMatch = -3
	Match = 1
	Gap = -10
	ScoreMatrix[1…M+1][1… N+1] = 0
	PointerMatrix[1…M+1][1…N+1] = 0
# PointerMatrix Status - 1: Up trace,2: Left trace,3: Diag trace
	FOR i = 1 TO M
		DO ScoreMatrix[i][1] = 0
			PointerMatrix[i][1] = 1
	FOR j = 1 TO N
		DO ScoreMatrix[1][j] = 0
			PointerMatrix[1][i] = 2

	// fill score matrix
	FOR i = 1 TO M
		FOR j = 1 TO N
			DO UpScore = ScoreMatrix[i+1][j] + Gap
				LeftScore = ScoreMatrix[i][j+1] + Gap
				IF x[i] = y[j]
					DO DiagScore = ScoreMatrix[i][j] + Match
				ELSE
					DO DiagScore = ScoreMatrix[i][j] +MisMatch
ScoreMatrix[i+1][j+1] = max(DiagScore, UpScore, LeftScore)
// mark trace
IF ScoreMatrix[i+1][j+1] = DiagScore
	DO PointerMatrix[i+1][j+1] = 3
ELSE IF ScoreMatrix[i+1][j+1] = LeftScore
	DO PointerMatrix[i+1][j+1] = 2
ELSE IF ScoreMatrix[i+1][j+1] = UpScore
	DO PointerMatrix[i+1][j+1] = 1

	// search best score and mark the aligned end position 
	MaxScore = -100000000
	StarPosj = N
	StarPosi = M
	FOR j = 1 TO N
		IF ScoreMatrix[M+1][j] >= MaxScore
			DO MaxScore = ScoreMatrix[M+1][j]
				StarPosj = j
	FOR i = 1 TO M
		IF ScoreMatrix[i][N+1] >= MaxScore
			DO MaxScore = ScoreMatrix[i][N+1]
				StarPosi = i
				StarPosj = N

	// get overlap sequences
	Overlapx = “”
	Overlapy =””
	i = StarPosi
	j = StarPosj 
	WHILE i!=1 and j!=1
		IF PointerMatrix[i][j] = 3 
			DO Overlapx = x[i] + Overlapx
				Overlapy = y[j] + Overlapy
				i -= 1
				j -= 1
		ELSE IF PointerMatrix[i][j] = 2
			DO	Overlapx = “-” + Overlapx
				Overlapy = y[j] + Overlapy
				j -=1
		ELSE IF PointerMatrix[i][j] = 1
			DO	Overlapx = x[i] + Overlapx
				Overlapy = “-” + Overlapy
				i -=1
	RETURN MaxScore, Overlapx, Overlapy
```					
* Q1 2). your analysis about the time and space complexity of your algorithm. 
time complexity: O(M*N)
space complexity: O(M*N)




* Q2
根据ssearch输出结果结构，修改代码正则表达式提取规则，结果如下，
```
$./ssearchparser.pl ssearch.out
Best hit to u1_human is: K08D10.3, with score 633, P-value 7.5e-18.
```
perl code
```
#!/usr/bin/perl -w
use strict;

my ($ssearch_out) = @ARGV;
my $usage = "This script is to get the best hit from ssearch output file with 1 input sequence.
usage: $0 <ssearch_output_file>
";
die $usage if @ARGV<1;

open(SSEARCHOUT,$ssearch_out)||die("open $ssearch_out error!\n");

my $query = "";
my $score = "";
my $p_value = "";
my $hit = "";
my $flag = 0;
while(<SSEARCHOUT>){
    chomp;
    if($flag == 0){
	if(/^Query:\s*(\w+)/){
	    $query = $1;
	}
	elsif(/^The best scores are:/){
	    $flag = 1;
	}
	else{
	    next;
	}
    }
    else{
	if(/^([\w\.\-]+)\s+.+\s+([0-9]+)\s+[0-9]+.[0-9]+\s+([0-9e\-\.]+)$/){
	    $hit = $1;
	    $score = $2;
	    $p_value = $3;
	    last;
	}
	else{
	    next;
	}
    }
}

close SSEARCHOUT;

print "Best hit to $query is: $hit, with score $score, P-value $p_value.\n";

exit;

```
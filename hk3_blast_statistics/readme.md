# Omics Big Data Spring 2019
Homework 3, week 3

1. Problem 1 (theory)  
Given the single-letter scoring system and sequence shown below:
```
A = +2, C = -1, G = -4, T = +2
TTACTGCGCCTTATAGCTATACGCTGTCGATCTGCGCAATTCCCCCCAATATCCCTCGGTTGATATTAC
```
A. What is the maximum segment score?  
B. What are the start and end points of the maximum-scoring segment (MSS)?  
C. Can the A+T composition of the maximum-scoring segment be deduced from its score
alone, or from its score and length?  
D. Suggest scoring systems for which the answer to 1C would be different.  

--------

2. Problem 2 (theory)  
Using the crude “hydrophobicity” scoring system shown below, if the 20 common amino acid codes are equiprobable in the sequences being searched, are Karlin-Altschul statistics applicable? Why or why not?
```
        AGILVFHPW = + 3
        CSTMY = 0
        NQ = -1
        RKDE = -2
```

----

3. Problem 3 (theory)  
A. What P-values (probabilities) are associated with the E-values (expected frequencies of
chance occurrence): 100., 10., 1., 0.1, 0.01?  
B. What trend is observed as E approaches infinity?  
C. What trend is observed a

----

4. Problem 4 (theory)  
A protein sequence of length 600 with typical amino acid composition was compared against the “NR” protein sequence database (70.028 billion amino acids in 192.341 million sequences on March 1, 2019). The scoring matrix was BLOSUM62 and the gap penalties were infinite. BLASTP computed precise values for λ, K and H of 0.304, 0.171, and 0.587, respectively. Clearly state any assumptions you make in answering the following questions. [Note: values forλand H are expressed here in units of “nats”, not bits. ]  
A. What is the highest alignment score to be expected between our query sequence and
the unrelated sequences in the database? [Hint: the expected high score corresponds to
the situation where E=1 ].  
B. How many bits of information are associated with expected high score?  
C. What is the the expected length of the highest scoring alignment between our query and
unrelated sequences?  
D. Answer the above questions ABC again, but this time using values forλ, K and H
appropriate for a second database search with the same query sequence, but using affine gap penalties u=10 and v=2. [Note: as used here, the penalty for a gap of length n is u+(n-1)v]. There are at least two sources to which you might refer for these values ofλ, k and H, including:  
The output from an appropriately parameterized BLASTP run, be it WU-BLASTP or NCBI-BLASTP.  
Table V. in Altchul and Gish, 1996.

Be sure to indicate how you obtained the values you used.

----

5. Problem 5 (theory)  
The first search in problem 4 uncovered a few high-scoring segment pairs (HSPs) having only marginal significance against the N-terminal portion of our query sequences. The search was repeated using just the N-terminal 100 residues (1/6th of the length of our original query sequence). This time, the same HSPs were reported, but their statistical significance was 12- to 40-fold higher than before (P-value were 12 to 40-fold lower), making the alignments appear to be somewhat significant.  

From the Karlin-Altschul equation, we expected the significance of the HSPs to increase (their P-values to decrease) when a shorter query sequence was used, but why did the statistical significance increase 12- to 40-fold, when the length of our query sequence was reduced by only 6-fold?

-----

6. Problem 6 (Perl/Python or any programming language)  
This computational experiment is intended to illustrate how the score and length expected for the maximal scoring segment(MSS) vary with the length and residue composition of the input model sequences. We’ll also see how the experimentally determined values for the expected score and length affect the values of K and relative entropy H.
For these experiments, you’ll need to write a Perl script that:  
 Generates random model sequences;  
 Finds the MSS in each model sequence;  
 Gathers statistics about the score and length of the MSSes;  
 Gathers statistics on the frequency of occurrence of the letter A in the MSS;  
 Outputs just three numbers: the expected score, expected length and expected
frequency of A(qA) in the MSS.  

The script must be capable of generating model sequences of arbitrary length and residue composition (although for any given experiment these will remain fixed). The process for selecting residues in the sequences should be i.i.d. Use Perl’s rand function as the random number generator. Using the model sequences so generated, the same script will need to search each sequence for the MSS, noting its score, length and fraction that is A(%A).

The alphabet for the model sequences should consist of the four letters A, C, G and T. You might find it convenient to use instead the digits 0, 1, 2 and 3 to represent these letters in your Perl script. (Maybe a numerical alphabet will even speed up your Perl script. You can test it.) Use any residue composition you wish for the nonuniform distribution, subject to the constraints that  
Pi >= 0.05 for all i  
All Pi differ from the uniform distribution by at least 0.1 The conditions for using Karlin-Altschul statistics still hold.  

The scoring system to use in all experiments is A=+1, C=G=T=-1. For all experimental results, have your script perform 104 iterations (i.e., collect statistics for MSSes from 104 model sequences). After iterating over the 104 sequences, have the script report just the expected (arithmetic mean) values for the score, length and %A.

Your fully contained Perl script (random sequence generator, MSS finder, statistics gather and reporter) should utilize 3 command line arguments in this order:  
a. The model sequence length;  
b. Whether to use a nonuniform distribution (0=> no, non-zero => yes);  
c. The number of model sequences to generate and examine (usually 10000).  
Note: full credit for this homework will require implementation of an O(N) (i.e., linear in the length of the model sequences) algorithm for finding the MSS.

You are forewarned that, even using an O(N) algorithm, if implemented inefficiently the lengthiest simulations performed here may take several hours to complete. To speed things up a good deal, it is acceptable to capitalize on the specialized conditions of the simulations and the limited output your program is required to produce, to create an efficient script that produces correct results far faster than a more general script. (For instance, no model sequences are to be reported by the program – only the results of computation on these sequences are to be reported – so it is unnecessary to expend time actually storing a randomly generated sequence, only to have to read it again to find the MSS. Why not scan a sequence for the MSS simultaneously with its generation?) In fact, if you program this problem in the C programming language instead of Perl, you might create a much faster program. However, there are no extra points for faster code or code written in C versus Perl, so don’t necessarily put your time into optimizing code if you haven’t already finished the requirements for answering all the questions. On the other hand, if you wait until the last day to start work on this problem, you may be forced to optimize your code in order to finish all simulations by the due date.

Independently of the experiments, you should populate the table below with values in units of nats computed to 2 decimal digits of accuracy. (Hint: identify an equation thatλsolves. ) Compute as well the target frequencies for the letter A in the MSS, reporting the value of qA for each distribution.

With the results of your experiments, fill in the table with the expected score, E(S), expected length, E(L), and expected fraction of A residues in the MSS, E(%A), when model sequence length of 102, 103, and 104, letters are used. Report values with 3 digits of precision for E(S), E(L) and E(%A). Using the experimental results, compute values for H (in units of nats) and K to 2 digits of precision. (Hint: solve E=Knme-λSfor K in the case of the expected high score. )
```
Distribution    Uniform -      -       Nonuniform  -   -
λ   
qA
Model Length, N 10^2    10^3    10^4    10^2      10^3 10^4
E(S)
E(L)
H
K
```

Questions:   
(1) What nonuniform distribution did you use?
A= C= G= T=    
(2) Which model sequence length yielded values for E(%A) closest to the theoretical target frequency qA? Why?  
(3) Which pair of model sequence lengths yielded results for H and K that are most similar? Why?   
(4) For both the uniform and nonuniform distributions and a single model sequence of length 10^4, what is the expected frequency of occurrence of a MSS having a score twice the expected high score?

------

Homework 3 

1.
A. 
14
B.
 start:38, end:68
C.
No.
can not be deduced from score alone, e.g. AAAAACCCCAAAAA(score=6), AAAAAGAAAAA(score=6)
can not be deduced from score and length e.g. TTACCTT(socre=8, length=7,AT%=5/7),  TTGAATT(socre=8,length=7,AT%=6/7)  
D. 
If A=T=+k, C=G=-m
let the length of A,T,C,G be la, lt ,lc ,lg, socre be s, length be l
	la + lt + lc + lg = l
	k*(la + lt) – m*(lc + lg) = s
	la, lt, lc, lg are positive integer
AT% = (la + lt)/l can be calculated from score and length


2.
No. It requires that the expected score sum is negative.(3*9+5*0-1*2-2*4>0)


3.
A.
```
P	E
1	100
0.99995460007024	10
0.63212055882856	1
0.09516258196404	0.1
0.0099501662508319	0.01
```
B. P is near 1, as E approaches infinity.
C. P is near E, as E approaches 0.
 


4. λ, K and H of 0.304, 0.171, and 0.587
(70.028 billion amino acids in 192.341 million sequences on March 1, 2019)
```
A. 

S = -1/lambda*ln(E/(KMN)) = -1/ 0.304*ln(1/0.171/600/(70.028*10^9)) = 97.37828595971635
B.
  
S’ = 31.369090654233243
C.
E(L) = log(KMN)/H=log(0.171*600*(70.028*10^9))/0.587= 50.43100329089228
D.
\lambda=0.266, K=0.04,H(nats)=0.24 (source: Table V. in Altchul and Gish, 1996.)
N=70.028*10^9
S = 105.82787529836482
S’ = 31.369090654233243
E(L) = 117.29256178902102
```

5.
A high-scoring alignment must have some length, and therefore can not begin near to the end of either of two sequences being compared. This "edge effect" may be corrected for by calculating an "effective length" for sequences. The BLAST programs implement such a correction(M’=M-E(L)). M’ may be small. 

6.
```
  the equation to solve \lambda
r = 4 , for which stands {A,C,G,T}
s = {+1,-1,-1,-1}
if uniform, p = {0.25,0.25,0.25},\lambda = ln(3)
if not uniform, I use p = {0.1, 0.4, 0.4, 0.1} which means P(A)=0.1,P(C)=0.4,P(G)=0.4,P(T)=0.1, \lambda=2ln(3) 
```
```
Distribution	Uniform	Nonuniform
\lambda	1.10	1.10	1.10	2.20	2.20	2.20
qA	4.000	6.971	10.146	1.778	2.928	4.103
Model Length, N	10^2	10^3	10^4	10^2	10^3	10^4
E(S)	3.277	5.311	7.410	1.704	2.740	3.791
E(L)	4.722	8.631	12.882	1.852	3.116	4.415
E(%A)	0.909	0.858	0.829	0.984	0.965	0.954
H	3.37	2.29	1.87	9.15	6.71	5.71
K	0.08	0.04	0.03	0.23	0.12	0.09
```
Questions:
(1) A=0.1,C=0.4,G=0.4,T=0.1
(2) Nonuniform model is close, after I compare qA/E(L) and E(%A) in each group. 

if qAj/Lj is equal(or similar) to each other for j=1,…,n then
 .
which means that qAj/Lj variance is small in nonuniform model.(E(%A) is higher, E(S),E(L) is smaller.)
(3) neither. 
(4) uniform: E(frequence(S>2E(S))) = 1.9
nonuniform: E(frequence(S>2E(S))) = 1.46


-------
* python code
```
# python3 om_hk3.py seqlength nonuniform seqnumber
import sys
import random

def strtolistint(string):
    # char to int 
    trans_dic = {'A':0,
    'C':1,
    'G':2,
    'T':3}
    return tuple([trans_dic[string[i]] for i in range(len(string))])


def freq(intlist,element=0):
    element_freq = sum([1 for i in intlist if i==element])
    return element_freq


def mss(seq,score):
    # transform to Kadane's algorithm, f[i]: the score of position i of seq
    n = len(seq)
    f = tuple([score[seq[i]] for i in range(n)])
    start = 0
    end = 0
    temp_sumscore = 0
    max_score = 0
    j = 0
    for i in range(n):
        temp_sumscore += f[i]
        if temp_sumscore>max_score:
            max_score = temp_sumscore
            start = j
            end = i 
        if temp_sumscore<0:
            temp_sumscore = 0
            j = i + 1      
    ele_freq = freq(seq[start:end+1])
    return max_score, end-start+1, ele_freq, ele_freq/(end-start+1)


def randomseq(length,uniform=1):
    if uniform:
        elements = (0,1,2,3)
        return tuple([random.choice(elements) for i in range(length)])
    else:
        p = (0.1,0.4,0.4,0.1) # pa, pc, pg ,pt
        elements = (0,1,1,1,1,2,2,2,2,3)
        return tuple([random.choice(elements) for i in range(length)])


if __name__ == '__main__':
    # get input 
    try:
        seqlength = int(sys.argv[1])
        nonuniform =  sys.argv[2]
        seqnumber = int(sys.argv[3])
        if nonuniform == '0':
            uniform = 1
        else:
            uniform = 0
    except:
        print('''python3 om_hk3.py seqlength nonuniform seqnumber
          seqlength <int> : the length of each sequence
          nonuniform <1 or 0> : use nonunifrom generation or not(0:uniform, 1:nonuniform)
          seqnumber <int> : the number of sequences to generate''')
        exit()

    #A,C,G,T score
    score = (1,-1,-1,-1)
    result = tuple([mss(randomseq(seqlength,uniform),score) for i in range(seqnumber)])
    print("E(S),E(L),E(qA),E(%A):")
    for outidx in range(4):
       print(sum([each[outidx] for each in result])/seqnumber,end=' ')


    #zzz = [0]*10
    #for i in range(10):
    #    result = tuple([mss(randomseq(seqlength,uniform),score) for i in range(seqnumber)])
        #for outidx in range(4):
        #print(sum([each[outidx] for each in result])/seqnumber)

    #    avgscore = sum([each[0] for each in result])/seqnumber
    #    zzz[i] = sum([1 for j in range(seqnumber) if result[j][0]>=avgscore])
    ##print(zzz)
    #print(sum(zzz)/10)  

    #print(max([result[j][0] for j in range(seqnumber)]))
    #s = 'TTACTGCGCCTTATAGCTATACGCTGTCGATCTGCGCAATTCCCCCCAATATCCCTCGGTTGATATTAC'
    #s_list_dig = strtolistint(s)
    #score,start,end = mss(s_list_dig,score)
    #print(score,start,end,s[start:end+1]) 
```

    

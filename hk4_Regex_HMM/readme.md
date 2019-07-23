# Omics Big Data Spring 2019
## Homework 4, week 4 Regex&HMM
1. Roll your own PROSITE pattern scanner in Perl or any programming language of your choice. Write a Perl script or a program with your favorite programming language that searches a protein sequence database with a PROSITE pattern.

Input: Your script should take two arguments: a PROSITE pattern string, and a filename for a FASTA format sequence database. For example:
% ./protsite.pl “[RK]-G-{EDRKHPCG}-[AGSCI]-[FY]-[LIVA]-x-[FYLM].” /nfs/wa
Be sure that your program can accept any PROSITE pattern gracefully. A full definition of PROSITE syntax is listed in the a file located in our course server at /share/home/ccwei/courses/2019/omics/hw4/pattern_syntax_rules. You will want to deal with the following pattern elements:  
[]: sets of allowed amino acids  
{}: sets of disallowed amino acids  
-: separator between pattern elements (x) and (x,y): length ranges  
<,>: anchor to N- or C-terminus  
.: end of pattern.  
(Hint: By far the easiest strategy will be to convert the PROSITE pattern to a Perl regular expression.)  
Output: for every match, your script should print out the name and description of the target sequence and the start/end point of the pattern match. Your script should be able to deal with multiple hits per target sequence. For example:
```
52 59 RU1A_HUMAN    U1 SMALL NUCLEAR RIBONUCLEOPROTEIN A (U1 SNRNP A PROTEIN).
138 145 PAB1_HUMAN  POLYADENYLATE-BINDING PROTEIN 1 (POLY(A) BINDING PROTEIN 1 
231 238 PAB1_HUMAN  POLYADENYLATE-BINDING PROTEIN 1 (POLY(A) BINDING PROTEIN 
233 340 PAB1_HUMAN  1 POLYADENYLATE-BINDING PROTEIN 1 (POLY(A) BINDING PROTEIN 1
```
2. The Rosencrantz and Guildenstern inference problem  
For some people, this may be a nontrivial problem. Leave time to do it properly. I think it’ll be worth your time... you should learn how HMM algorithms work.

In the Tom Stoppard play Rosencrantz and Guildenstern Are Dead, the play opens with Guildenstern flipping coins that always come up heads, against all odds, leading to a discussion of probabilities that sets the tone for the play. In the play, Rosencrantz has no doubt about what the next flip will produce, and he complains that the flipping game is boring. Consider the following more complicated problem, in which Guildenstern’s flips are a hidden Markov process, and Rosencrantz must use HMM theory to figure out what’s going on. Stoppard, I’m sure, would approve.

Here are the rules to a game that our more interesting Guildenstern plays every day. Guildenstern has two coins. The first coin is fair, and the second coin is biased. He starts with the first coin. After each flip, there is a small probability that he will secretly switch to the second coin. Then after each flip of the second coin, there is a small probability that he will stop flipping coin for the day. (He never switches back to the fair coin.) He asks Rosencrantz to guess when he switched to the second coin, based solely on the observed sequence of heads and tails.

Clearly Rosencrantz can model this as a hidden Markov process. The two coins behave as two HMM states, emitting H and T with some probability; the small switching probability behaves as a transition probability from the first to the second state, and the small ending probability behaves as a transition probability to the usual HMM special end state.

An example of a string that Rosencrantz might see:   
THHHHTTHTTHTTHHTTTTHHTHTHTHTTTTTTHHHHHHTTHTHTHTHHHH

Rosencrantz’s problem: he can’t tell for sure just from looking at a string like this when Guildenstern switched to the biased coin. But, if he learns HMM theory, he can make a really good guess (or, as we say to make it sounds less like guesswork, a statistical inference).

Rosencrantz sees a different sequence of heads and tails every day that they play the game. A file of the sequences for 100 days of the game is the example file /share/home/ccwei/courses/2019/omics/hw4/example; the format is one sequence per line, and some of the lines can be fairly long. Have a look at the file... do you think you could guess where the switch happens in each sequence? It’s pretty hard to do by eye (I can’t do it).

Part a. The HMM architecture.  
Draw the HMM architecture (states, emissions and state transitions) that corresponds to Guildenstern’s game. You don’t need to turn this picture in; but you need it for your own use for the next parts of the question, and we’ll see its structure anyway in the two Perl scripts you write. (You might include a description of the architecture in the comments of your Perl script, though – it’ll help us understand what you’re doing, if you do something we don’t expect!).

Part b. Viterbi HMM alignment  
Let’s make it a bit easy on Rosencrantz: let’s assume he’s played the game so much that he knows Guildenstern’s parameters. The first (fair) coin has emission probabilities p(H)=p(T)=0.5. The second (biased) coin has emission probabilities p(H) = 0.8, p(T) = 0.2. The probability of switching from the first to the second coin is 0.01; the probability of ending the sequence after flipping the second coin is 0.05.  
Implement a Viterbi dynamic programming alignment procedure in Perl for the HMM you drew above, using these parameters. Have your script read a file such as the example data file /share/home/ccwei/courses/2019/omics/hw4/example and, for each sequence, find the

maximum likelihood state path. For each sequence, have your Perl script print out the maximum likelihood guess for which position was the first flip of the biased coin (e.g. a number from 2..N for a sequence of length N).  
Part c. Forward-Backward HMM alignment and posterior decoding  
Of course, the maximum likelihood position is still just a guess; its count can easily be wrong. We can use Forward-Backward and posterior decoding to get an even more detailed inference for Rosencrantz.  
Implement Forward and Backward dynamic programming procedures in Perl for the HMM, using the same parameters as above. Have your Perl script read a file such as the example data file /share/home/ccwei/courses/2019/omics/hw4/example), and then, for each sequence in the file, use the Forward and Backward variables to calculate the posterior probability distribution over positions of the first flip of the biased coin... e.g. for a given position, what is the probability that this position is the correct guess? For each sequence, have your Perl script print out the best (most probable) five positions, followed by the summed posterior probability of all five (e.g. is the probability that one of the five is right).  

Hints:  
In this case, the position with maximum a posterior probability will agree with the position indicated by the Viterbi trace. This is not always the case for more complicated HMMs.
For debugging, it might be helpful to write an auxiliary Perl script that generates strings from Guildenstern’s “HMM”. The script I used, for example, is   /share/home/ccwei/courses/2019/omics/hw4/generate.pl. This way you can generate data and keep track of the correct answer (when the model switched to the biased coin), and see how well your inference engines guess.  
Another debugging tip: in this case, there’s only N-1 possible state paths for a sequence of length N. You can actually enumerate them all, and calculate likelihoods explicitly! This isn’t terribly efficient but it is a good way to double check the answers from your DP algorithms. Normally, the number of possible state paths for an HMM is exponentially large, and you can’t do this.

-----

1. Regex
```python
import re
import sys

def openfasta(fastafile):
    seqdict = dict()
    for line in open(fastafile):
        if line[0] == '>':
            seqname = line[1:].split(' ')[0]
            seqdict[seqname] = [line.rstrip(),'']
        else:
            seqdict[seqname][1] += line.rstrip()
    return seqdict


def patterntorestring(pattern):
    # check point(.)
    if pattern[-1] == '.':
        blocks = pattern.split('.')[0]
    else:
        print('Long long ago, you forgot a point(.).')
        return
    # split -
    blocks = blocks.split('-')
    restring = ''
    for block in blocks:
        # check x
        if 'x' in block:
            block = block.replace('x', '.')
            if len(block) == 1:
            	continue
            # check ( and )
            if block[1] == '(' and block[-1] == ')':
                templength = block[2:-1].split(',')
                if len(templength) == 2:
                    minl = templength[0]
                    maxl = templength[1]
                    tempstring = '{'+minl+','+maxl+'}'
                    restring += block[0] + tempstring
                elif len(templength) == 2:
                    length = templength[0]
                    tempstring = '{'+length+'}'
                    restring += block[0] + tempstring
        else:
            if ')' in block or '(' in block:
                print('Not valid pattern!')
                return

        # check > and <
        if '<' in block:
            block = '^' + block.replace('<','')
        if '>' in block:
            block = block.replace('>','') + '$'

        # check single character
        if len(block) == 1:
            restring += block
            continue

        # check [,],{ and }
        if block[0]  == '[' and block[-1] == ']':
            tempstring = '[' + block[1:-1] + ']'
            restring += tempstring
        elif block[0]  == '{' and block[-1] == '}':
            tempstring = '[^' + block[1:-1] + ']'
            restring += tempstring
    return restring


def findpattern(patrern,fastafile):
    printnum = 0
    restring = patterntorestring(patrern)
    if not restring:
        exit()
    rereg =re.compile(r''+str(restring),re.S)
    #print(rereg)
    seqdict = openfasta(fastafile)
    for seq in seqdict:
        results = rereg.finditer(seqdict[seq][1])
        if results:
            for match in results:
                idxs = match.span()
                content = match.group()
                seqanotation = seqdict[seq][0][1:]
                seqname = ''
                seqdiscribe = ''
                if len(seqanotation.split('\t')) > 1:
                    seqname = seqanotation.split('\t')[0]
                    seqdiscribe = ' '.join(seqanotation.split('\t')[1:])
                elif len(seqanotation.split(' ')) > 1:
                    seqname = seqanotation.split(' ')[0]
                    seqdiscribe = ' '.join(seqanotation.split('\t')[1:])
                printnum += 1
                print(str(idxs[0]+1)+'\t'+str(idxs[1]+1)+'\t'+seqname+'\t'+content+'\t'+seqdiscribe)
    print('Total {n} Hits!'.format(n=str(printnum)))


#if __name__ == '__main__':
try:
    patt = sys.argv[1]
    file = sys.argv[2]
    findpattern(patt,file)
except Exception as e:
    print(e)
    print('''python3 hk4_1.py [pattern] [fastafile]
    [pattern] <string>: the Re pattern(e.g. '[RK]-G-{EDRKHPCG}-[AGSCI]-[FY]-[LIVA]-x-[FYLM].')
    [fastafile] <string>: the location of fasta file
    ''')

```

------

2. HMM
```python
import sys

def startposition(stateseq,state='b'):
    return stateseq.find(state)+1


def viterbi(observations,transiton,emission):
    observations += '$'

    maxp_m = [] #possiblity matrix
    maxt_m = [] #trace matrix
    for k in range(len(observations)):
        if k == 0:
            f = emission['f' + observations[k]]
            maxp_m.append([f,0,0])
            maxt_m.append(['f','b','e'])
        else:
            # init
            tof = 0
            tob = 0
            toe = 0
            temptracef = ''
            temptraceb = ''
            temptracee = ''
            if observations[k] == '$':
                # this state is e
                be = maxp_m[k - 1][1] * transiton['be'] * emission['e' + observations[k]]
                toe = be
                temptracee = maxt_m[k - 1][1] + 'e'
            else:
                # this state is f
                ff = maxp_m[k - 1][0] * transiton['ff'] * emission['f' + observations[k]]
                tof = ff
                temptracef = maxt_m[k-1][0] + 'f'

                # this state is b
                fb = maxp_m[k - 1][0] * transiton['fb'] * emission['f' + observations[k]]
                bb = maxp_m[k - 1][1] * transiton['bb'] * emission['b' + observations[k]]
                if fb > bb:
                    tob = fb
                    temptraceb = maxt_m[k - 1][0] + 'b'
                else:
                    tob = bb
                    temptraceb = maxt_m[k - 1][1] + 'b'

            maxp_m.append([tof,tob,toe])
            maxt_m.append([temptracef,temptraceb,temptracee])
    u =len(maxp_m)-1
    if maxp_m[u][0] >= maxp_m[u][1]:
        if maxp_m[u][2] > maxp_m[u][0]:
            return(startposition(maxt_m[u][2]))
        else:
            return(startposition(maxt_m[u][0]))
    else:
        return(startposition(maxt_m[u][1]))


def gettop(longlist,n=5):
    n = min(n,len(longlist)) # exists length < 5
    sortedlist = sorted(longlist,reverse=True)
    # top n index
    topidx = []
    for i in range(len(longlist)):
        for k in range(n):
            if longlist[i] == sortedlist[k]:
                topidx.append(i+1)
    return topidx,sum(sortedlist[:n])


def f (observations,transiton,emission):
    observations += '$'
    # forward
    forward = []
    fpre = {}
    for idx in range(len(observations)-1):
        fcur = {}
        if idx == 0:
            fpre_sumf = 1
            fpre_sumb = 0
        else:
            fpre_sumf = fpre['f'] * transiton['ff']
            fpre_sumb = fpre['f'] * transiton['fb'] + fpre['b'] * transiton['bb']
        fcur['f'] = fpre_sumf * emission['f'+observations[idx]]
        fcur['b'] = fpre_sumb * emission['b'+observations[idx]]
        forward.append(fcur)
        fpre = fcur
    p_forward = fcur['b'] * transiton['be']

    # backward
    backward = []
    bpre = {}
    for idx in range(len(observations)-2,-1,-1):
        bcur = {}
        if idx == len(observations)-2:
            bpre_sumf = 0
            bpre_sumb = 1 * transiton['be']
        else:
            bpre_sumb = bpre['b'] * transiton['bb']
            bpre_sumf = bpre['b'] * transiton['fb'] + bpre['f'] * transiton['ff']
        bcur['f'] = bpre_sumf * emission['f'+observations[idx]]
        bcur['b'] = bpre_sumb * emission['b'+observations[idx]]
        backward = [bcur] + backward
        bpre = bcur
    p_backward = bcur['f'] # is equal with p_forward

    # get posterior
    posterior = []
    for i in range(len(forward)):
        posterior.append({'f':forward[i]['f']*backward[i]['f']/p_forward,'b':forward[i]['b']*backward[i]['b']/p_forward})
    p_startpos = [0]
    for i in range(1,len(posterior)):
        p_startpos.append(posterior[i]['b']/emission['b'+observations[i]] - posterior[i-1]['b']/emission['b'+observations[i-1]])
    return(gettop(p_startpos))


### parameters
# alphbet H, T, $(end)
alphabet = ['H','T','$']
# state Fair ,Bias, End
state = ['f','b','e']
# transition
s = {}
s['fb'] = 0.01
s['ff'] = 0.99
s['be'] = 0.05
s['bb'] = 0.95
# emission
p = {}
p['fH'] = 0.5
p['fT'] = 0.5
p['bH'] = 0.8
p['bT'] = 0.2
p['e$'] = 1
#if __name__ == '__main__':

try:
    vresult = []
    fbresult = []
    standard = []
    #Viterbi:
    for line in open(sys.argv[1]):
        vresult.append(viterbi(line.rstrip(),s,p))

    #Forward-Backward:(postion & posterior probability)
    for line in open(sys.argv[1]):
        fbresult.append(forwardbackward(line.rstrip(),s,p))

    #print('Debug#########')
    #for line in open(sys.argv[1]+'.sta'):
    #    standard.append(startposition(line.strip(),'2'))

    print('viterbi result, forward-backward result(top 5,sum_P)')#, standard, whether standard is in forward-backward result')
    for i in range(len(vresult)):
        print(vresult[i],fbresult[i])#,standard[i],int(standard[i]) in fbresult[i][0],sep='\t')
except:
    print('''
    python3 hk4_2.py [file]
    file <string>: observations file''')
```
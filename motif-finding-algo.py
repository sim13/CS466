import numpy as np
def motif_algorithm(path):
    file_seq = open(str(path) + '/sequences.fa', 'r')
    file_len = open(str(path) + '/motiflength.txt', 'r')    
    motiflength = int(file_len.readline())
    sequences = []
    for line in file_seq:
        if(line[0]!='>'):
            line.rstrip('\n')
            sequences.append(line)
            
    file_seq.close()
    file_len.close()
    bestMotif = [0] * len(sequences)
    s = [0] * len(sequences)
    s1 = 0
    s2 = 0
    s3 = 0
    i1 = [0,0]
    i2 = [0,0]
    i3 = [0,0]
    lmer1 = ''
    lmer2 = ''
    lmer3 = ''
    for p in arange(0, len(sequences[0])-motiflength):
        s[0] = p
        for q in arange(0, len(sequences[0])-motiflength):
            s[1] = q
            seq1_s = sequences[0][s[0]:s[0]+motiflength]
            seq2_s = sequences[1][s[1]:s[1]+motiflength]
            
            seq1_b = sequences[0][bestMotif[0]:bestMotif[0]+motiflength]
            seq2_b = sequences[1][bestMotif[1]:bestMotif[1]+motiflength]
            
            score_s = 0
            score_b = 0
            for i in arange(0, len(seq1_s)):
                if(seq1_s[i]==seq2_s[i]):
                    score_s += 1
                if(seq1_b[i]==seq2_b[i]):
                    score_b += 1
                
            if(score_s >= score_b):
                if(score_s > s1):
                    s3 = s2
                    lmer3 = lmer2
                    i3 = i2
                    s2 = s1
                    lmer2 = lmer1
                    i2 = i1
                    s1 = score_s
                    lmer1 = seq1_s
                    i1 = [s[0], s[1]]
                elif(score_s > s2):
                    s3 = s2
                    lmer3 = lmer2
                    i3 = i2
                    s2 = score_s
                    lmer2 = seq1_s
                    i2 = [s[0], s[1]]
                elif(score_s > s3):
                    s3 = score_s
                    lmer3 = seq1_s
                    i3 = [s[0], s[1]]
                else:
                    continue
   # print lmer1, lmer2, lmer3
    elmer = ''                
    for x in arange(0, len(sequences[0])-motiflength):
        s[2] = x
        seq_s = sequences[2][x:x+motiflength]
        seq_b = sequences[2][bestMotif[2]:bestMotif[2]+motiflength]
        score_s1 = 0
        score_s2 = 0
        score_s3 = 0
        score_b1 = 0
        score_b2 = 0
        score_b3 = 0
        for n in arange(0, len(seq_s)):
            if(seq_s[n]==lmer1[n]):
                score_s1 += 1
            if(seq_s[n]==lmer2[n]):
                score_s2 += 1
            if(seq_s[n]==lmer3[n]):
                score_s3 += 1
            if(seq_b[n]==lmer1[n]):
                score_b1 += 1
            if(seq_b[n]==lmer2[n]):
                score_b2 += 1
            if(seq_b[n]==lmer3[n]):
                score_b3 += 1
        score_s = max(score_s1, score_s2, score_s3)
        '''if(score_s==score_s1):
            score_b = score_b1
        elif(score_s==score_s2):
            score_b = score_b2
        elif(score_s==score_s3):
            score_b = score_b3
        '''
        score_b = max(score_b1, score_b2, score_b3)
        
        if(score_s > score_b):
            bestMotif[2] = s[2]
            #if(s[2]==4):
                #print score_s, score_b, score_s1, score_s2, score_s3, s[2]
            #if(s[2]==452):
                #print score_s, score_b, score_s1, score_s2, score_s3, s[2]
            if(score_s==score_s1):
                elmer = lmer1
                bestMotif[0] = i1[0]
                bestMotif[1] = i1[1]
               
            elif(score_s==score_s2):
                elmer = lmer2
                bestMotif[0] = i2[0]
                bestMotif[1] = i2[1]
               
            elif(score_s==score_s3):
                elmer = lmer3
                bestMotif[0] = i3[0]
                bestMotif[1] = i3[1]

    #print bestMotif[0:2]
    #print i1, i2, i3
    s[0] = bestMotif[0]
    s[1] = bestMotif[1]
    s[2] = bestMotif[2]
    
    for i in arange(3, len(sequences)):
        for p in arange(0, len(sequences[0])-motiflength):
            s[i] = p
            seq_s = sequences[i][p:p+motiflength]
            seq_b = sequences[i][bestMotif[i]:bestMotif[i]+motiflength]
            score_s = 0
            score_b = 0
            for n in arange(0, len(seq_s)):
                if(seq_s[n]==elmer[n]):
                    score_s += 1
             
                if(seq_b[n]==elmer[n]):
                    score_b += 1
        
            
            if(score_s > score_b):
                bestMotif[i] = s[i]
        s[i] = bestMotif[i]
    aln_mat = []
    for i in arange(0, len(sequences)):
        aln_mat.append(sequences[i][bestMotif[i]:bestMotif[i]+motiflength])
    
    tran = np.asarray(aln_mat)
    tran = zip(*tran)
    pwm = []
    for i in arange(0, motiflength):
        sh = []
        sh.append(tran[i].count('A'))
        sh.append(tran[i].count('C'))
        sh.append(tran[i].count('G'))
        sh.append(tran[i].count('T'))
        pwm.append(sh)
    file_prsites = open(str(path) + '/predictedsites.txt', 'w')
    file_prmotif = open(str(path) + '/predictedmotif.txt', 'w')
    for x in arange(0, len(sequences)):
        file_prsites.write(str(bestMotif[x]) + '\n')
    file_prmotif.write('> PMOTIF    ' + str(motiflength) + '\n')
    for i in arange(0, motiflength):
        for x in arange(0, 4):
            file_prmotif.write(str(pwm[i][x]))
            if(x!=3):
                file_prmotif.write('    ')
            else:
                file_prmotif.write('\n')
    file_prmotif.write('<')
    file_prmotif.close()
    file_prsites.close()
    #print bestMotif
    #print pwm
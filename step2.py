import numpy as np
import os
import math


def predict_motif(PWM, DIR = './sample_data/', file = 'predictedmotif.txt'):
    ML = PWM.shape[0]
    np.savetxt(DIR+file, PWM, fmt="%f", header=f">PREDICTED MOTIF {ML}", comments = "") 


def predict_site(sites, DIR = './sample_data/', file = 'predictedsites.txt'):
    with open(DIR+file, 'w+') as f:
        f.write("\n".join(str(idx) for idx in sites))
        f.close()


def read_data(DIR='./sample_data/', ML_file = "motiflength.txt", seq_file = 'sequences.fa'):

    with open(DIR+seq_file) as f:
        next(f)
        seqs = [line.strip() for line in f]
        f.close()
    
    with open(DIR+ML_file) as f:
        ML = f.readline()
        f.close()
    ML = int(ML)
    return seqs, ML


def get_PWM(sites):
    ML = len(sites[0])
    PWM = np.zeros((ML, 4), dtype=np.float)
    for i in range(ML):
        for seq in sites:
            PWM[i, dic[seq[i]]] += 1
    PWM = PWM / len(sites)
    return PWM


def get_scores(PWM, sub_seqs):
    scores = [sum(np.log(PWM[i][dic[c]]) for i, c in enumerate(sub_seq)) for sub_seq in sub_seqs]
    return scores


def gibbs(seqs, ML):
    INDEX = [np.random.randint(0, len(seq) - ML + 1) for seq in seqs]
    prev = None
    cnt = 0
    while INDEX != prev: 
        cnt += 1
        if cnt >= 100: break # break if over 100 rounds.
        prev = INDEX[:]
        for i, seq in enumerate(seqs):
            # compute PWM using all sequence except the current sequence
            PWM = get_PWM([s[idx : idx + ML] for j, (s, idx) in enumerate(zip(seqs, INDEX)) if j != i])
            # find the optimal binding sites in the current sequence
            sub_seqs = [seq[idx : idx + ML] for idx in range(len(seq) - ML + 1)]
            scores = get_scores(PWM, sub_seqs)
            optimal_idx = np.argmax(scores)
            INDEX[i] = optimal_idx
    # print(cnt)
    return INDEX, get_PWM([x[j : j + ML] for x, j in zip(seqs, INDEX)])


if __name__ == "__main__":

    dic = {"A":0, "C":1, "G":2, "T":3}
    for i in range(1,71):
        DIR = f"./data/dataset{i}/"
        seqs, ML = read_data(DIR=DIR)
        sites, PWM = gibbs(seqs, ML)
        predict_motif(PWM, DIR=DIR)
        predict_site(sites, DIR=DIR)

    # seqs, ML = read_data(DIR='./sample/')
    # sites, PWM = gibbs(seqs, ML)
    # print(sites)
    # print(PWM)
    # predict_motif(PWM, DIR='./sample/')
    # predict_site(sites, DIR='./sample/')

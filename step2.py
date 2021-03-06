from cmath import inf
import socketserver
import numpy as np
import os
import math
import argparse


dic = {"A":0, "C":1, "G":2, "T":3}

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



def get_seq_ICs(PWM, sub_seqs):
    m = PWM/0.25
    log_m = np.log2(m, out=np.zeros_like(m), where=(m!=0))
    IC_matrix = np.multiply(PWM, log_m)
    # get the IC for each sub_seq, given the current PWM (theta(t) in the lecture)
    ICs = [sum([IC_matrix[i][dic[c]] for i, c in enumerate(sub_seq)]) for sub_seq in sub_seqs]
    # print(IC_matrix)
    return np.array(ICs)


def compute_total_IC(PWM):
    m = PWM/0.25
    log_m = np.log2(m, out=np.zeros_like(m), where=(m!=0))
    IC_matrix = np.multiply(PWM, log_m)
    ICPCs = np.sum(IC_matrix)
    return ICPCs


def gibbs_version2(seqs, ML):
    # generate the inital INDEX for site of each sequence
    INDEX = [np.random.randint(0, len(seq) - ML + 1) for seq in seqs]
    prev = None
    cnt = 0
    best_ic = 0
    Best_index = None
    while INDEX != prev: 
        cnt += 1
        if cnt >= 20: break # break if over loop over 100 iteration 
        prev = INDEX[:]
        for i, seq in enumerate(seqs):
        # lst = np.random.permutation(len(seqs))
        # for i in lst:
            seq = seqs[i]
            # compute PWM using all sequence except the current sequence
            PWM = get_PWM([s[idx : idx + ML] for j, (s, idx) in enumerate(zip(seqs, INDEX)) if j != i])
            # find the optimal binding sites in the current sequence
            # get all sub-sequence of length ML of the current sequence
            sub_seqs = [seq[idx : idx + ML] for idx in range(len(seq) - ML + 1)]
            # calculate information content for each sub-sequecne
            ICs = get_seq_ICs(PWM, sub_seqs) # shape: (493, )
            # take the one with highest IC as the new optimal index
            optimal_idx = np.argmax(ICs)
            INDEX[i] = optimal_idx
            temp_PWN = get_PWM([x[j : j + ML] for x, j in zip(seqs, INDEX)])
            cur_ic = compute_total_IC(temp_PWN)
            if best_ic < cur_ic:
                best_ic = cur_ic
                Best_index = INDEX
    return Best_index, get_PWM([x[j : j + ML] for x, j in zip(seqs, Best_index)])


def get_QxPx(PWM, sub_seqs):
    m = PWM/0.25
    # log_m = np.log2(PWM, out=np.zeros_like(PWM), where=(PWM!=0))
    # log_m[log_m == 0] = -inf
    prods = [np.prod([m[i][dic[c]] for i, c in enumerate(sub_seq)])  for sub_seq in sub_seqs]
    # pdd = [np.sum([log_m[i][dic[c]] for i, c in enumerate(sub_seq)]) for sub_seq in sub_seqs]
    # print(np.array(pdd).shape)
    return np.array(prods)


def gibbs(seqs, ML):
    # generate the inital INDEX for site of each sequence
    INDEX = [np.random.randint(0, len(seq) - ML + 1) for seq in seqs]
    prev = None
    cnt = 0
    best_ic = 0
    Best_index = None
    # while INDEX != prev: 
    while True:
        cnt += 1
        if cnt >= 20: break # break if over loop over 100 iteration 
        prev = INDEX[:]
        # for i, seq in enumerate(seqs):
        lst = np.random.permutation(len(seqs))
        for i in lst:
            seq = seqs[i]
            # compute PWM using all sequence except the current sequence
            PWM = get_PWM([s[idx : idx + ML] for j, (s, idx) in enumerate(zip(seqs, INDEX)) if j != i])
            # find the optimal binding sites in the current sequence
            # get all sub-sequence of length ML of the current sequence
            sub_seqs = [seq[idx : idx + ML] for idx in range(len(seq) - ML + 1)]
            # calculate information content for each sub-sequecne
            QxPx = get_QxPx(PWM, sub_seqs) # shape: (493, )
            if all(v == 0 for v in QxPx): continue
            # take the one with highest IC as the new optimal index
            # Select a idx according to the value of Qx/Px
            chosen_idx = np.random.choice([idx for idx in range(len(seq)-ML+1)], p=QxPx/QxPx.sum(), size=1)
            INDEX[i] = chosen_idx[0]
            temp_PWN = get_PWM([x[j : j + ML] for x, j in zip(seqs, INDEX)])
            cur_ic = compute_total_IC(temp_PWN)
            if best_ic < cur_ic:
                best_ic = cur_ic
                Best_index = INDEX
    return Best_index, get_PWM([x[j : j + ML] for x, j in zip(seqs, Best_index)])



def main(args):
    dic = {"A":0, "C":1, "G":2, "T":3}


    method = gibbs
    if args.method == "mod_gibbs":
        method = gibbs_version2

    n_data = args.n_data
    folder = "data"
    for i in range(1, n_data * 7 + 1):
        print(f"Computing the {i}-th data")
        DIR = f"./{folder}/dataset{i}/"
        seqs, ML = read_data(DIR=DIR)
        sites, PWM = method(seqs, ML)
        # total_IC = compute_ICPC(PWM)
        # print(f"ICPC:{total_IC/ML}\n")
        predict_motif(PWM, DIR=DIR)
        predict_site(sites, DIR=DIR)

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description='Motif-Finding')
    arg_parser.add_argument('--method', type=str, default="gibbs", help="method to use")
    arg_parser.add_argument('--n_data', type=int, default=10, help="number of data for each parameter combination")
    args = arg_parser.parse_args()
    main(args)


    # seqs, ML = read_data(DIR='./sample/')
    # sites, PWM = gibbs_version2(seqs, ML)
    # print(sites)
    # total_IC = compute_total_IC(PWM)
    # print(total_IC)
    # predict_motif(PWM, DIR='./sample/')
    # predict_site(sites, DIR='./sample/')

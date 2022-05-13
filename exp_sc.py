from audioop import avg
import argparse
import os
from evaluator import Evaluator, relative_entropy
from step2 import read_data, gibbs, predict_motif, predict_site, gibbs_version2
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Regular expression for sorting file names.
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


def perform_motif_finding(DIR, method):
    seqs, ML = read_data(DIR=DIR)
    sites, PWM = method(seqs, ML)
    predict_motif(PWM, DIR=DIR)
    predict_site(sites, DIR=DIR)


def plot_a_graph(keys, r_dic, e_dic, o_dic, labels, fn):
    plt.figure(figsize=(16, 7), dpi=200)
    plt.subplot(1, 3, 1)
    plt.boxplot([e_dic[keys[0]], e_dic[keys[1]], e_dic[keys[2]]], labels=labels)
    plt.title("relative entropy")
    plt.subplot(1, 3, 2)
    plt.boxplot([o_dic[keys[0]], o_dic[keys[1]], o_dic[keys[2]]], labels=labels)
    plt.title("overlapped sites")
    plt.subplot(1, 3, 3)
    plt.boxplot([r_dic[keys[0]], r_dic[keys[1]], r_dic[keys[2]]], labels=labels)
    plt.title("run times")
    plt.savefig(f"./graphs/{fn}", bbox_inches='tight')


def plot_graphs(run_time_dict, entropy_dict, overlapped_dict, algor):
    combinations = ["sc_20", "sc_50", "sc_100"]
    label_sc = ['SC-20', 'SC-50', 'SC-100']
    keys_sc = ['sc_20', 'sc_50', 'sc_100']
    # plot_a_graph(keys_icpc, run_time_dict, entropy_dict, overlapped_dict, label_icpc, fn= f"ICPC_{algor}.png")
    # plot_a_graph(keys_ml, run_time_dict, entropy_dict, overlapped_dict, label_ml, fn= f"ML_{algor}.png")
    plot_a_graph(keys_sc, run_time_dict, entropy_dict, overlapped_dict, label_sc, fn= f"SC_test_{algor}.png")
    
    
 


def main(args):
    # Disable this term for debugging
    eval_runtime = True
    evaluator = Evaluator()

    dic = {"A":0, "C":1, "G":2, "T":3}
    method = gibbs
    algor = "gibbs"
    if args.method == "mod_gibbs":
        method = gibbs_version2
        algor = "mod_gibbs"

    folder = "data_sc"
    n_data = args.n_data

    # To get running time, evaluation should be simultaneous with predition
    if eval_runtime:
        run_time_dict = defaultdict(list)
        combinations = ["sc_20", "sc_50", "sc_100"]
        for cnt, comb in enumerate(combinations):
            for i in range(cnt*n_data+1, (cnt+1)*n_data+1):
                DIR = f"./{folder}/dataset{i}/"
                evaluator.set_start_time()
                perform_motif_finding(DIR, method)
                run_time_dict[comb].append(evaluator.get_running_time())

    overlapped_dict = defaultdict(list)
    entropy_dict = defaultdict(list)

    combinations = combinations = ["sc_20", "sc_50", "sc_100"]
    for cnt, comb in enumerate(combinations):
        for i in range(cnt*n_data+1, (cnt+1)*n_data+1):
            DIR = f"./{folder}/dataset{i}/"
            motif_loss, sites_loss, _, _ = evaluator.evaluate(DIR)
            entropy_dict[comb].append(motif_loss)
            overlapped_dict[comb].append(sites_loss)
    
    plot_graphs(run_time_dict, entropy_dict, overlapped_dict, algor)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description='Evaluation')
    arg_parser.add_argument('--method', type=str, default="gibbs")
    arg_parser.add_argument('--n_data', type=int, default=10, help="number of data for each parameter combination")
    args = arg_parser.parse_args()
    main(args)
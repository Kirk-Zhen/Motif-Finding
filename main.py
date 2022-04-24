from audioop import avg
import os
from evaluator import Evaluator, relative_entropy
from step2 import read_data, gibbs, predict_motif, predict_site
import re
import numpy as np
import matplotlib.pyplot as plt

# Regular expression for sorting file names.
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]






if __name__ == "__main__":
    # Disable this term for debugging
    eval_runtime = True
    evaluator = Evaluator()
    # To get running time, evaluation should be simultaneous with predition
    if eval_runtime:
        run_times_1 = []
        run_times_15 = []
        run_times_2 = []
        dic = {"A":0, "C":1, "G":2, "T":3}
        for i in range(1,51):
            DIR = f"./data_icpc/dataset{i}/"
            evaluator.set_start_time()
            seqs, ML = read_data(DIR=DIR)
            sites, PWM = gibbs(seqs, ML)
            predict_motif(PWM, DIR=DIR)
            predict_site(sites, DIR=DIR)
            run_times_1.append(evaluator.get_running_time())
        for i in range(51,101):
            DIR = f"./data_icpc/dataset{i}/"
            evaluator.set_start_time()
            seqs, ML = read_data(DIR=DIR)
            sites, PWM = gibbs(seqs, ML)
            predict_motif(PWM, DIR=DIR)
            predict_site(sites, DIR=DIR)
            run_times_15.append(evaluator.get_running_time())
        for i in range(101,151):
            DIR = f"./data_icpc/dataset{i}/"
            evaluator.set_start_time()
            seqs, ML = read_data(DIR=DIR)
            sites, PWM = gibbs(seqs, ML)
            predict_motif(PWM, DIR=DIR)
            predict_site(sites, DIR=DIR)
            run_times_2.append(evaluator.get_running_time())
        run_times = [run_times_1, run_times_15, run_times_2]
    
    motif_losses = []
    sites_losses = []
    ICPC_list = []
    
    num_overlapped_sites_1 = []
    num_overlapped_sites_15 = []
    num_overlapped_sites_2 = []
    
    relative_entropy_1 = []
    relative_entropy_15 = []
    relative_entropy_2 = []
    
    for i in range(1,51):
        DIR = f"./data_icpc/dataset{i}/"
        motif_loss, sites_loss, _, _ = evaluator.evaluate(DIR)
        relative_entropy_1.append(motif_loss)
        num_overlapped_sites_1.append(sites_loss)
    for i in range(51, 101):
        DIR = f"./data_icpc/dataset{i}/"
        motif_loss, sites_loss, _, _ = evaluator.evaluate(DIR)
        relative_entropy_15.append(motif_loss)
        num_overlapped_sites_15.append(sites_loss)
    for i in range(101,151):
        DIR = f"./data_icpc/dataset{i}/"
        motif_loss, sites_loss, _, _ = evaluator.evaluate(DIR)
        relative_entropy_2.append(motif_loss)
        num_overlapped_sites_2.append(sites_loss)
    
    relative_entropys = [relative_entropy_1, relative_entropy_15, relative_entropy_2]
    num_overlapped_sites = [num_overlapped_sites_1, num_overlapped_sites_15, num_overlapped_sites_2]
        
    # draw boxplots on all three datasets on all relative entropynum overlapped sites, and run times
    plt.subplot(1, 3, 1)
    plt.boxplot(relative_entropys, labels=['ICPC-1.0', 'ICPC-1.5', 'ICPC-2.0'])
    plt.title("relative entropy")
    plt.subplot(1, 3, 2)
    plt.boxplot(num_overlapped_sites, labels=['ICPC-1.0', 'ICPC-1.5', 'ICPC-2.0'])
    plt.title("overlapped sites")
    plt.subplot(1, 3, 3)
    plt.boxplot(run_times, labels=['ICPC-1.0', 'ICPC-1.5', 'ICPC-2.0'])
    plt.title("run times")
    plt.savefig('ICPC_impact.png')
    
    
    
    # for dataset in sorted(os.listdir(ROOTDIR), key=natural_keys):
    #     if dataset.startswith("dataset"):
    #         print("----------------------------")
    #         print(dataset)
    #         motif_loss, sites_loss, icpc_motif, icpc_predictedmotif \
    #             = evaluator.evaluate(os.path.join("./data/", dataset))

    #         print("motif_loss:", motif_loss)
    #         print("sites_loss:", sites_loss)
    #         print("ICPC(original): ", round(icpc_motif,1))
    #         print("ICPC(predicted):", icpc_predictedmotif)
    #         motif_losses.append(motif_loss)
    #         sites_losses.append(sites_loss)
    #         ICPC_list.append(icpc_predictedmotif)
    # print("============================")
    # print("Average motif_loss:", np.mean(motif_losses))
    # print("Average sites_loss:", np.mean(sites_losses))
    # print("Average ICPC(predicted):", np.mean(ICPC_list))
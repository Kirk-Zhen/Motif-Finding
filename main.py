from audioop import avg
import os
from evaluator import Evaluator
from step2 import read_data, gibbs, predict_motif, predict_site
import re
import numpy as np

# Regular expression for sorting file names.
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]






if __name__ == "__main__":
    # Disable this term for debugging
    eval_runtime = False
    evaluator = Evaluator()
    # To get running time, evaluation should be simultaneous with predition
    if eval_runtime:
        run_times = []
        dic = {"A":0, "C":1, "G":2, "T":3}
        for i in range(1,71):
            DIR = f"./data/dataset{i}/"
            evaluator.set_start_time()
            seqs, ML = read_data(DIR=DIR)
            sites, PWM = gibbs(seqs, ML)
            predict_motif(PWM, DIR=DIR)
            predict_site(sites, DIR=DIR)
            run_times.append(evaluator.get_running_time())
    
    motif_losses = []
    sites_losses = []
    ICPC_list = []
    
    for dataset in sorted(os.listdir("./data/"), key=natural_keys):
        if dataset.startswith("dataset"):
            print("----------------------------")
            print(dataset)
            motif_loss, sites_loss, icpc_motif, icpc_predictedmotif \
                = evaluator.evaluate(os.path.join("./data/", dataset))

            print("motif_loss:", motif_loss)
            print("sites_loss:", sites_loss)
            print("ICPC(original): ", round(icpc_motif,1))
            print("ICPC(predicted):", icpc_predictedmotif)
            motif_losses.append(motif_loss)
            sites_losses.append(sites_loss)
            ICPC_list.append(icpc_predictedmotif)
    print("============================")
    print("Average motif_loss:", np.mean(motif_losses))
    print("Average sites_loss:", np.mean(sites_losses))
    print("Average ICPC(predicted):", np.mean(ICPC_list))
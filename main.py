import os
from evaluator import Evaluator
from step2 import read_data, gibbs, predict_motif, predict_site








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
    
    for dataset in os.listdir("./data/"):
        if dataset.startswith("dataset"):
            print("----------------------------")
            print(dataset)
            motif_loss, sites_loss = evaluator.evaluate(os.path.join("./data/", dataset))
            print("motif_loss:", motif_loss)
            print("sites_loss:", sites_loss)
            motif_losses.append(motif_loss)
            sites_losses.append(sites_loss)
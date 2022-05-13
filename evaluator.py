import os
import time
import numpy as np
from step2 import compute_total_IC

def relative_entropy(p, q, esp=0.00001):
    # Scale p and q to [0, 1]
    p = p / np.linalg.norm(p) + esp
    q = q / np.linalg.norm(q) + esp
    rn = 0
    for i in range(len(p)):
        if p[i] != 0 and q[i] != 0:
            rn += p[i] * np.log(p[i] / q[i])
    return rn

class Evaluator(object):
    def __init__(self) -> None:
        self._start_time = 0
        self.motif_length = 0
        self.path = ""
        
    def _get_motif_length(self):
        motiflength_file = os.path.join(self.path, "motiflength.txt")
        self.motif_length = int(np.loadtxt(motiflength_file))
        
        return self.motif_length
    
    def _motif_loss(self):
        """
        Relative Entropy between “motif.txt” and “predictedmotif.txt”
        """
        motif_file = os.path.join(self.path, "motif.txt")
        predictedmotif_file = os.path.join(self.path, "predictedmotif.txt")
        motif = np.loadtxt(motif_file, dtype=float, skiprows=1)
        predictedmotif = np.loadtxt(predictedmotif_file, dtype=float, skiprows=1)
        entropy_list = []
        for l in range(len(motif)):
            entropy_list.append(relative_entropy(motif[l], predictedmotif[l]))
        return np.mean(entropy_list)
            
    def _sites_loss(self):
        """
        Number of overlapping positions or overlapping sites
        between “sites.txt” and “predictedsites.txt”
        """
        sites_file = os.path.join(self.path, "sites.txt")
        predictedsites_file = os.path.join(self.path, "predictedsites.txt")
        sites = np.loadtxt(sites_file, dtype=int)
        predictedsites = np.loadtxt(predictedsites_file, dtype=int)

        sites_list = []
        for site in sites:
            cur_site = [site+i for i in range(self.motif_length)]
            sites_list.append(cur_site)

        predictedsites_list = []
        for psite in predictedsites:
            cur_site = [psite+i for i in range(self.motif_length)]
            predictedsites_list.append(cur_site)

        overlap = 0
        for ori_site, pred_site in zip(sites_list, predictedsites_list):
            # print(ori_site, pred_site)
            for s in pred_site:
                if s in ori_site:
                    overlap += 1
        return overlap/len(sites)
        
        
    def _get_ICPC_difference(self):
        motif_file = os.path.join(self.path, "motif.txt")
        predictedmotif_file = os.path.join(self.path, "predictedmotif.txt")
        motif = np.loadtxt(motif_file, dtype=float, skiprows=1)
        predictedmotif = np.loadtxt(predictedmotif_file, dtype=float, skiprows=1)
        icpc_motif = compute_total_IC(motif)/self.motif_length
        icpc_predictedmotif = compute_total_IC(predictedmotif)/self.motif_length
        return icpc_motif, icpc_predictedmotif

    def set_start_time(self):
        self._start_time = time.time()
        
    def get_running_time(self):
        """
        Running time
        """
        running_time = time.time() - self._start_time
        self._start_time = 0

        return running_time
    
    def evaluate(self, dataset_path):
        self.path = dataset_path
        self._get_motif_length()
        # relative entropy
        motif_loss = self._motif_loss()
        # overlap
        site_loss = self._sites_loss()
        icpc_motif, icpc_predictedmotif = self._get_ICPC_difference()
        
        return motif_loss, site_loss, icpc_motif, icpc_predictedmotif
        
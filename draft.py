from urllib.request import proxy_bypass
import numpy as np


def single_random_seq(SL):
    """
    Generate one random sequence with length SL.
    """
    return "".join([np.random.choice(["A","C","G","T"]) for _ in range(SL)])


# 2
def multi_random_seq(SL, SC):
    """
    Generate SC random sequences. 
    Each random sequence has length SL.
    """
    return [single_random_seq(SL) for _ in range(SC)]


# 3, 8
def generate_PWM(ML, ICPC):
    """
    Generate a random motif (PWM) of length ML, 
    with total information content being ICPC * ML.
    """
    P = {2.0: 1, 1.5: 0.9245, 1.0:0.8105}
    prob = P[ICPC]
    res_p = (1 - prob)/3 
    PWM = np.full((ML, 4), res_p, dtype=np.float)
    PWM[np.arange(ML), np.random.randint(0, 4, ML)] = prob
    return PWM

def write_motif(PWM, motif_idx = 1, DIR = './sample_data/', file = 'motif.txt'):
    ML = PWM.shape[0]
    np.savetxt(DIR+file, PWM, fmt="%f", header=f">MOTIF{motif_idx} {ML}", comments = "") 


def PWM_to_motif(PWM):
    """ Get motif e.g., "CCGT" from PWM """
    idx = np.argmax(PWM, axis=1)
    dic = {0:"A", 1:"C", 2:"G", 3:"T"}
    return ''.join([dic[i] for i in idx])





# 4
def generate_binding_sites(PWM, SC):
    """Generate SC binding sites"""
    return [random_binding(PWM) for _ in range(SC)]


def random_binding(PWM):
    ML = PWM.shape[0]
    return "".join([np.random.choice(["A","C","G","T"], p=PWM[i]/np.sum(PWM[i])) for i in range(ML)])


# 5, 7
def plant_site_single(seq, site):
    site_length = len(site)
    idx = np.random.randint(0, len(seq)-site_length+1)
    return idx, seq[:idx]+site+seq[idx+site_length:]

def plant_site(seqs, sites, DIR = './sample_data/', file = 'sites.txt'):
    """
    “Plant” one sampled site at a random location.
    And write the locations in the file.
    """
    location = []
    new_seqs = []

    for i,seq in enumerate(seqs):
        idx, new_seq = plant_site_single(seq, sites[i])
        location.append(idx)
        new_seqs.append(new_seq)

    # Write locations in the file.
    with open(DIR+file, 'w+') as f:
        f.write("\n".join(str(idx) for idx in location))
        f.close()

    return new_seqs


# 6
def write_to_FASTA(seqs, DIR = './sample_data/', file = 'sequences.fa'):
    """
    Write sequences into a FASTA format file “sequences.fa”
    """
    with open(DIR+file, 'w+') as f:
        f.write(f">{DIR+file}\n")
        f.write("\n".join(str(seq) for seq in seqs))
        f.close()


# 9
def write_ML(ML, DIR = './sample_data/', file = 'motiflength.txt'):
    """
    write down the motif length.
    """
    with open(DIR+file, 'w+') as f:
        f.write(str(ML))
        f.close()
        




if __name__ == "__main__":
    ICPC = 1.5
    ML = 8
    SL = 20
    SC = 10
    motif_idx = 0


    seqs = multi_random_seq(SL, SC)
    PWM = generate_PWM(ML, ICPC)
    mtfs = generate_binding_sites(PWM, SC)
    new_seqs = plant_site(seqs, mtfs)
    write_to_FASTA(new_seqs)
    write_ML(ML)
    write_motif(PWM, motif_idx)


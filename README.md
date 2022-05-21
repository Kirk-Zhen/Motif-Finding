# Gibbs Sampling for Motif-Finding


##  Environment and Dependencies
|Package|Version|
|:---|:---|
|||
|numpy|1.21.6
|matplotlib|3.2.2
|scipy|1.4.1

To install all requirements:
```bash
pip install -r requirements.txt
```


<!-- ## Framework:  
<img src="fig/MasEnc.png" width="200"/><img src="/fig/BiteNet.png" width="400"/>  -->

## Step 1: Generate Data

Execute the following command to generate data:
```bash
python step1.py --n_data 10
```
`n_data`: the number of data per parameter combination. 

As is stated in the project requirement, `n_data = 10`.  

Data will be generated in the `data` folder.


## Step 2: Run Motif Finding Algorithm

Execute the following command to run the Gibbs Sampling Algorithm for Motif Finding:
```bash
python step2.py --n_data 10 --method gibbs
```
`n_data`: must be consistent with Step 1. 

`method`: parameter for the algorithm. `gibbs` for Gibbs Sampling Algorithm.


## Step 3: Evaluation

Execute the following command to generate the graphs for performance evaluation, graphs will be generated in the `graphs` folder.
```bash
python step3.py --n_data 10 --method gibbs
```
`n_data`:  must be consistent with Step 1.

`method`: parameter for the algorithm. `gibbs` for Gibbs Sampling Algorithm. 


## Reproduce the Results
Executing Step 3 will reproduce the results of `Figure 1`, `Figure 2` and `Figure 3` in `Result.pdf`. 
##### To reproduce Figure 4:
```bash
python exp_sc.py --n_data 10 --method gibbs
```
##### To reproduce Figure 5:
```bash
python exp_ml.py --n_data 10 --method gibbs
```
##### To reproduce Figure 6:
```bash
python exp_icpc.py --n_data 10 --method gibbs
```

All results will be generated in `graphs` folder.
### Deterministic version of Gibbs

We also implemented a modified, and more deterministic, and more time-efficient version of Gibbs Sampling. But the performance of this algorithm is not ideal as expected. 

If you want to run experiments on this algorithm, you can edit the parameter `--method gibbs` as `--method mod_gibbs`. The results for such algorithm will be also generated in `graphs` folder.

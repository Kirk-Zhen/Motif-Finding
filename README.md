# Motif-Finding



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
`n_data`: <font color="red"> must be consistent with Step 1</font>. 
`method`: parameter for the algorithm. `gibbs` for Gibbs Sampling Algorithm, `gibbs2` for a modified version of Gibbs Sampling Algorithm.


## Step 3: Evaluation

Execute the following command to generate the graphs for performance evaluation.
```bash
python step3.py --n_data 10 --method gibbs
```
`n_data`: <font color="red"> must be consistent with Step 1</font>. 
`method`: parameter for the algorithm. `gibbs` for Gibbs Sampling Algorithm, `gibbs2` for a modified version of Gibbs Sampling Algorithm.
#Evolutionary Bilevel Optimization based on Covariance Matrix Adaptation

This is the code to reproduce the experimental results in:

[1]X. He, Y. Zhou, and Z. Chen, “Evolutionary Bilevel Optimization based on Covariance Matrix Adaptation,” *IEEE Transactions on Evolutionary Computation*, pp. 1–1, 2018.

## File structure

In total, three files are contained in the root directory:

1. BLCMAES.m: The source code of the proposed BL-CMA-ES algorithm. 
2. BLEA_runExperiment.m: The code to launch the experiment.
3. BLEA_D: A folder containing all necessary files. Please note that it contains some files provided in BLEAQ2. So please cite the following paper if you decide to use the code formally: [2]A. Sinha, P. Malo, and K. Deb, “Evolutionary algorithm for bilevel optimization using approximations of the lower level optimal solution mapping,” *European Journal of Operational Research*, vol. 257, no. 2, pp. 395–411, Mar. 2017.



## Run an experiment



### Create an test instance

The first thing you need to know is to get a well-defined test instance.

Currently, this code contains two benchmark suites, namely SMD and TP. SMD contains 12 problems while TP contains 10.

Two problems from the real world including "GoldMining" and "DecisionMaking" are also available.

To create a predefined problem, use the code:

```matlab
BI_list = getBLOPinfo('SMD',1:4,5)
```

This will return (2+3)-dimensional SMD1 to SMD4 in an array. Each array element is a BI structure which is the valid input parameter of the  function BLCMAES.m. This function accepts three parameter: benchmark name, function No. list, and the dimensionality.

### Specify the stopping criteria 

At both levels, you should specify two stopping criteria.

1. maximum number of function evaluations allowed (denoted by UmaxFEs and LmaxFEs)
2. maximum number of function evaluations consumed in a stagnation (denoted by UmaxImprFEs and LmaxImprFEs)

### Run 

Almost all experimental settings have been defined in the file BLEA_runExperiment.m

Please launch this file and observe the console output.

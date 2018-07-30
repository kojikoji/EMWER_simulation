# EMWER_simulation
Simulation pipeline for evaluating EMWER performance. This pipeline simulate Evolve and Resequence experiments and provide sync format file (see detail in https://sourceforge.net/p/popoolation2/wiki/Tutorial/). This pipeline basically assume that each locus in the experiments is independent of each other.

## Installation
This pipeline depends on [MimicrEE](https://sourceforge.net/projects/mimicree/) and [Poisson sampler script](https://github.com/handetopa/BBGP/blob/master/scripts/simulation/poisson-3fold-sample.py). Hence, you firstly need to install them using shell scripts as below:


```shell-session
$ bash download.sh
```

## Usage
You can execute this program by typing as below:


```shell-session
$ python script/conduct_simulation.py -o <OUTPUT SYNC FILE> -v <SELECTION INFO FILE>
```
OUTPUT SYNC FILE: Sync format file including the simulation results.
SELECTION INFO FILE: File including chromosome (first column), genomic position (2nd column), Variant allele (3rd column), selection coefficients (4th column) and dominance (5th column.


You can see how to specify other experimental conditions by typing as below:


```shell-session
$ python script/conduct_simulation.py -h
```

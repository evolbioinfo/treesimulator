# treesimulator

Simulation of rooted phylogenetic trees under a given Multitype Birth–Death (MTBD) model, with or without contact tracing (CT).

[//]: # ([![DOI:10.1093/sysbio/syad059]&#40;https://zenodo.org/badge/DOI/10.1093/sysbio/syad059.svg&#41;]&#40;https://doi.org/10.1093/sysbio/syad059&#41;)
[![GitHub release](https://img.shields.io/github/v/release/evolbioinfo/treesimulator.svg)](https://github.com/evolbioinfo/treesimulator/releases)
[![PyPI version](https://badge.fury.io/py/treesimulator.svg)](https://pypi.org/project/treesimulator/)
[![PyPI downloads](https://shields.io/pypi/dm/treesimulator)](https://pypi.org/project/treesimulator/)
[![Docker pulls](https://img.shields.io/docker/pulls/evolbioinfo/treesimulator)](https://hub.docker.com/r/evolbioinfo/treesimulator/tags)

## MTBD
The MTBD models were introduced by Stadler & Bonhoeffer [[_Philos. Trans. R. Soc. B_ 2013]](https://royalsocietypublishing.org/doi/10.1098/rstb.2012.0198).

An MTBD model with m states has 

m(m-1) state transition rate parameters:
* μ<sub>ij</sub> -- transition rate from state i to state j (1 ≤ i, j ≤ m; i ≠ j), where μ<sub>ij</sub> ≥ 0

m<sup>2</sup> transmission rate parameters:
* λ<sub>ij</sub> -- transmission rate from state i (donor) to state j (recipient) (1 ≤ i, j ≤ m), where λ<sub>ij</sub> ≥ 0

m removal (becoming non-infectious) rate parameters:
* ψ<sub>i</sub> -- removal rate of state i (1 ≤ i ≤ m), where ψ<sub>i</sub> ≥ 0

m sampling probability upon removal parameters:
* p<sub>i</sub> -- probability to sample the pathogen of an individual in state i upon removal (1 ≤ i ≤ m), where 0 < p<sub>i</sub> ≤ 1

## Contact Tracing (CT)
Contact tracing adds two parameters to the initial MTBD model:

* υ -- probability to notify contacts upon sampling
* φ -- notified contact removal and sampling rate: φ >> ψ<sub>i</sub> ∀i (1 ≤ i ≤ m). The pathogen of a notified contact is sampled automatically (with a probability of 1) upon removal. 


We pay particular interest to the classical BD model, the BD Exposed-Infectious (BDEI) model, 
and BD with superspreading (BDSS), 
as they are described in [[Voznica _et al._ 2021]]((https://www.biorxiv.org/content/10.1101/2021.03.11.435006v1)), and to their -CT versions.


## BD
3 parameters:
* λ -- transmission rate
* ψ -- removal rate
* p -- sampling probability upon removal

Epidemiological parameters:
* R<sub>0</sub>=λ/ψ -- reproduction number
* 1/ψ -- infectious time


## BD-CT
5 parameters:
* λ -- transmission rate
* ψ -- removal rate
* p -- sampling probability upon removal
* υ -- probability to notify contacts upon sampling
* φ -- notified contact removal and sampling rate: φ >> ψ 


Epidemiological parameters:
* R<sub>0</sub>=λ/ψ -- reproduction number
* 1/ψ -- infectious time
* 1/φ -- notified contact removal time

## BDEI
2 states: 
* E, exposed, i.e. infected but not yet infectious
* I, infectious

4 parameters:
* μ -- transition rate from E to I (becoming infectious)
* λ -- transmission rate from I to E
* ψ -- removal rate of I
* p -- sampling probability upon removal

Epidemiological parameters:
* R<sub>0</sub>=λ/ψ -- reproduction number
* 1/ψ -- infectious time
* 1/μ -- incubation period


## BDEI-CT
2 states: 
* E, exposed, i.e. infected but not yet infectious
* I, infectious

6 parameters:
* μ -- transition rate from E to I (becoming infectious)
* λ -- transmission rate from I to E
* ψ -- removal rate of I
* p -- sampling probability upon removal
* υ -- probability to notify contacts upon sampling
* φ -- notified contact removal and sampling rate: φ >> ψ 

Epidemiological parameters:
* R<sub>0</sub>=λ/ψ -- reproduction number
* 1/ψ -- infectious time
* 1/μ -- incubation period
* 1/φ -- notified contact removal time

## BDSS
2 compartments: 
* N, standard infectious individual
* S, superspreader

6 parameters:
* λ<sub>nn</sub> -- transmission rate from N to N
* λ<sub>ns</sub> -- transmission rate from N to S
* λ<sub>sn</sub> -- transmission rate from S to N
* λ<sub>ss</sub> -- transmission rate from S to S

    (with a constraint that λ<sub>ss</sub>/λ<sub>ns</sub>=λ<sub>sn</sub>/λ<sub>nn</sub>)
* ψ -- removal rate of S and of N (the same)
* p -- sampling probability upon removal (the same for N and S)

Epidemiological parameters:
* R<sub>0</sub>=(λ<sub>nn</sub> + λ<sub>ss</sub>)/ψ -- reproduction number
* 1/ψ -- infectious time
* X=λ<sub>ss</sub>/λ<sub>ns</sub>=λ<sub>sn</sub>/λ<sub>nn</sub> -- super-spreading transmission ratio
* f=λ<sub>ss</sub>/(λ<sub>sn</sub> + λ<sub>ss</sub>) -- super-spreading fraction

## BDSS-CT
2 states: 
* N, standard infectious individual
* S, superspreader

8 parameters:
* λ<sub>nn</sub> -- transmission rate from N to N
* λ<sub>ns</sub> -- transmission rate from N to S
* λ<sub>sn</sub> -- transmission rate from S to N
* λ<sub>ss</sub> -- transmission rate from S to S

    (with a constraint that λ<sub>ss</sub>/λ<sub>ns</sub>=λ<sub>sn</sub>/λ<sub>nn</sub>)
* ψ -- removal rate of S and of N (the same)
* p -- sampling probability upon removal (the same for N and S)
* υ -- probability to notify contacts upon sampling
* φ -- notified contact removal and sampling rate: φ >> ψ 

Epidemiological parameters:
* R<sub>0</sub>=(λ<sub>nn</sub> + λ<sub>ss</sub>)/ψ -- reproduction number
* 1/ψ -- infectious time
* X=λ<sub>ss</sub>/λ<sub>ns</sub>=λ<sub>sn</sub>/λ<sub>nn</sub> -- super-spreading transmission ratio
* f=λ<sub>ss</sub>/(λ<sub>sn</sub> + λ<sub>ss</sub>) -- super-spreading fraction
* 1/φ -- notified contact removal time


## Installation

There are 4 alternative ways to run __treesimulator__ on your computer: 
with [docker](https://www.docker.com/community-edition), 
[apptainer](https://apptainer.org/),
in Python3, or via command line (requires installation with Python3).


### Installation in python3 or command-line

You could either install python (version 3.6 or higher) system-wide and then install treesimulator via pip:
```bash
sudo apt install -y python3 python3-pip python3-setuptools python3-distutils
pip3 install treesimulator
```

or alternatively, you could install python (version 3.6 or higher) and treesimulator via [conda](https://conda.io/docs/) (make sure that conda is installed first). 

(Optional) to install treesimulator in a new conda environment (e.g., called _phyloenv_ below), first create and activate the environment:
```bash
conda create --name phyloenv python=3.6
conda activate phyloenv
```

Install treesimulator with conda
```bash
conda install treesimulator
```


### Basic usage in a command line
If you installed __treesimulator__ in a conda environment (here named _phyloenv_), do not forget to first activate it, e.g.

```bash
conda activate phyloenv
```


#### BD and BD-CT
The following command simulates a tree with 200-500 tips under the BD model, with λ=0.5, ψ=0.25, p=0.5,
and saves it to the file tree.nwk, while saving the parameters to the comma-separated file params.csv:
```bash
generate_bd --min_tips 200 --max_tips 500 \
--la 0.5 --psi 0.25 --p 0.5 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BD-CT model, with λ=0.5, ψ=0.25, p=0.5, φ=2.5, υ=0.2, 
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bd --min_tips 200 --max_tips 500 \
--la 0.5 --psi 0.25 --p 0.5 \
--phi 2.5 --upsilon 0.2 --max_notified_contacts 1 \
--nwk tree.nwk --log params.csv
```
To see detailed options, run:
```bash
generate_bd --help
```

#### BDEI and BDEI-CT
The following command simulates a tree with 200-500 tips under the BDEI model, with μ=1, λ=0.5, ψ=0.25, p=0.5, 
and saves it to the file tree.nwk, while saving the parameters to the comma-separated file params.csv:
```bash
generate_bdei --min_tips 200 --max_tips 500 \
--mu 1 --la 0.5 --psi 0.25 --p 0.5 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDEI-CT model, with μ=1, λ=0.5, ψ=0.25, p=0.5, φ=2.5, υ=0.2, 
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdei --min_tips 200 --max_tips 500 \
--mu 1 --la 0.5 --psi 0.25 --p 0.5 \
--phi 2.5 --upsilon 0.2 --max_notified_contacts 1 \
--nwk tree.nwk --log params.csv
```
To see detailed options, run:
```bash
generate_bdei --help
```


#### BDSS and BDSS-CT
The following command simulates a tree with 200-500 tips under the BDSS model, 
with λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=0.5, λ<sub>ss</sub>=1.5, ψ=0.25, p=0.5, 
and saves it to the file tree.nwk, while saving the parameters to the comma-separated file params.csv:
```bash
generate_bdss --min_tips 200 --max_tips 500 \
--la_nn 0.1 --la_ns 0.3 --la_sn 0.5 --la_ss 1.5 --psi 0.25 --p 0.5 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDSS-CT model, 
with λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=0.5, λ<sub>ss</sub>=1.5, ψ=0.25, p=0.5, φ=2.5, υ=0.2, 
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdss --min_tips 200 --max_tips 500 \
--la_nn 0.1 --la_ns 0.3 --la_sn 0.5 --la_ss 1.5 --psi 0.25 --p 0.5 \
--phi 2.5 --upsilon 0.2 --max_notified_contacts 1 \
--nwk tree.nwk --log params.csv
```
To see detailed options, run:
```bash
generate_bdss --help
```

#### User-defined MTBD and MTBD-CT models
The following command simulates a tree with 200-500 tips under a generic MTBD model, with two states A and B, 
with μ<sub>aa</sub>=0.5, μ<sub>ab</sub>=0.6, μ<sub>ba</sub>=0.7, μ<sub>bb</sub>=0.8, 
λ<sub>aa</sub>=0.1, λ<sub>ab</sub>=0.2, λ<sub>ba</sub>=0.3, λ<sub>bb</sub>=0.4, 
ψ<sub>a</sub>=0.05, ψ<sub>b</sub>=0.08,
p=<sub>a</sub>0.15, p=<sub>b</sub>0.65,
and saves it to the file tree.nwk, while saving the parameters to the comma-separated file params.csv:
```bash
generate_mtbd --min_tips 200 --max_tips 500 \
--states A B \
--transition_rates 0.5 0.6 0.7 0.8 \
--transmission_rates 0.1 0.2 0.3 0.4 \
--removal_rates 0.05 0.08 \
--sampling_probabilities 0.15 0.65 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under a generic MTBD model, with two states A and B, 
with μ<sub>aa</sub>=0.5, μ<sub>ab</sub>=0.6, μ<sub>ba</sub>=0.7, μ<sub>bb</sub>=0.8, 
λ<sub>aa</sub>=0.1, λ<sub>ab</sub>=0.2, λ<sub>ba</sub>=0.3, λ<sub>bb</sub>=0.4, 
ψ<sub>a</sub>=0.05, ψ<sub>b</sub>=0.08,
p=<sub>a</sub>0.15, p=<sub>b</sub>0.65,
φ=2.5, υ=0.2, 
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_mtbd --min_tips 200 --max_tips 500 \
--states A B \
--transition_rates 0.5 0.6 0.7 0.8 \
--transmission_rates 0.1 0.2 0.3 0.4 \
--removal_rates 0.05 0.08 \
--sampling_probabilities 0.15 0.65 \
--phi 2.5 --upsilon 0.2 --max_notified_contacts 1 \
--nwk tree.nwk --log params.csv
```
To see detailed options, run:
```bash
generate_mtbd --help
```


#### Basic usage in Python3
To simulate trees with 200-500 tips under the above models and settings:

```python3
from treesimulator.generator import generate
from treesimulator import save_forest
from treesimulator.mtbd_models import Model, BirthDeathModel, BirthDeathExposedInfectiousModel, BirthDeathWithSuperSpreadingModel, CTModel

# BD and BD-CT
bd_model = BirthDeathModel(p=0.5, la=0.5, psi=0.25)
print(bd_model.get_epidemiological_parameters())
[bd_tree], _, _ = generate(bd_model, min_tips=200, max_tips=500)
save_forest([bd_tree], 'BD_tree.nwk')
bdpn_model = CTModel(model=bd_model, upsilon=0.2, notified_removal_rate=2.5)
[bdpn_tree], _, _ = generate(bdpn_model, min_tips=200, max_tips=500)
save_forest([bdpn_tree], 'BDCT_tree.nwk')

# BDEI and BDEI-CT
bdei_model = BirthDeathExposedInfectiousModel(p=0.5, mu=1, la=0.5, psi=0.25)
print(bdei_model.get_epidemiological_parameters())
[bdei_tree], _, _ = generate(bdei_model, min_tips=200, max_tips=500)
save_forest([bdei_tree], 'BDEI_tree.nwk')
bdeipn_model = CTModel(model=bdei_model, upsilon=0.2, notified_removal_rate=2.5)
[bdeipn_tree], _, _ = generate(bdeipn_model, min_tips=200, max_tips=500)
save_forest([bdeipn_tree], 'BDEICT_tree.nwk')

# BDSS and BDSS-CT
bdss_model = BirthDeathWithSuperSpreadingModel(p=0.5, la_nn=0.1, la_ns=0.3, la_sn=0.5, la_ss=1.5, psi=0.25)
print(bdss_model.get_epidemiological_parameters())
[bdss_tree], _, _ = generate(bdss_model, min_tips=200, max_tips=500)
save_forest([bdss_tree], 'BDSS_tree.nwk')
bdsspn_model = CTModel(model=bdss_model, upsilon=0.2, notified_removal_rate=2.5)
[bdsspn_tree], _, _ = generate(bdsspn_model, min_tips=200, max_tips=500)
save_forest([bdsspn_tree], 'BDSSCT_tree.nwk')

# MTBD and MTBD-CT
mtbd_model = Model(states=['A', 'B'], transition_rates=[[0.5, 0.6], [0.7, 0.8]],
                   transmission_rates=[[0.1, 0.2], [0.3, 0.4]],
                   removal_rates=[0.05, 0.08], ps=[0.15, 0.65])
[mtbd_tree], _, _ = generate(mtbd_model, min_tips=200, max_tips=500)
save_forest([mtbd_tree], 'MTBD_tree.nwk')
mtbdpn_model = CTModel(model=mtbd_model, upsilon=0.2, notified_removal_rate=2.5)
[mtbdpn_tree], _, _ = generate(mtbdpn_model, min_tips=200, max_tips=500)
save_forest([mtbdpn_tree], 'MTBDCT_tree.nwk')
```

### Run with apptainer

Once [apptainer](https://apptainer.org/docs/user/latest/quick_start.html#installation) is installed, 
run the following command:

```bash
apptainer run docker://evolbioinfo/treesimlator
```

This will launch a terminal session within the container, 
in which you can run treesimulator following the instructions for the command line ("Basic usage in a command line") above.





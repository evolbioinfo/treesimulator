# treesimulator

Simulation of rooted phylogenetic trees under a given Multi-Type Birth–Death (MTBD) model, 
with or without contact tracing (CT), 
and with or without Skyline. 

## Preprint
treesimulator is described in the S1 Appendix of the following article:

Anna Zhukova, Olivier Gascuel. Accounting for contact tracing in epidemiological birth-death models. medRxiv 2024.09.09.24313296; doi:[10.1101/2024.09.09.24313296](https://www.medrxiv.org/content/10.1101/2024.09.09.24313296v2.supplementary-material)

[![DOI:10.1101/2024.09.09.24313296](https://zenodo.org/badge/DOI/10.1101/2024.09.09.24313296.svg)](https://doi.org/10.1101/2024.09.09.24313296)
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

and a meta-parameter κ, which defines how many most recent contacts can be notified by each index case. Each is notified independently, with the probability υ.

## Skyline
Skyline was introduced by Stadler _et al._ [[PNAS 2013]](https://doi.org/10.1073/pnas.1207965110) 
and extended to MTBD by Kühnert _et al._ [[MBE 2016]](https://doi.org/10.1093/molbev/msw064). 
It enables piece-wise constant parameter value changes. 
To use a skyline with k models, one needs to specify k-1 model change times t<sub>1</sub>, ..., <sub>k-1</sub>, 
and k sets of model parameters (see above).
At time 0 the simulation starts with model 1, it switches to models 2 at time t<sub>1</sub>, etc.
All the models in the Skyline must have the same states. 
CT-related parameters can also change at skyline changing times, 
in that case if some skyline intervals do not have CT, υ=0 and any value for φ must be specified for them. 
The same κ value is shared among all the skyline -CT models.


We pay particular interest to the classical BD model, the BD Exposed-Infectious (BDEI) model, 
and BD with superspreading (BDSS), 
as they are described in [[Voznica _et al._ 2021]]((https://www.biorxiv.org/content/10.1101/2021.03.11.435006v1)), and to their -CT(κ) versions.


## BD
3 parameters:
* λ -- transmission rate
* ψ -- removal rate
* p -- sampling probability upon removal

Epidemiological parameters:
* R<sub>0</sub>=λ/ψ -- reproduction number
* 1/ψ -- infectious time


## BD-CT(κ)
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


## BDEI-CT(κ)
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

## BDSS-CT(κ)
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


#### BD, BD-CT(κ) and BD-CT(κ)-Skyline
The following command simulates a tree with 200-500 tips under the BD model, with λ=0.5, ψ=0.25, p=0.5,
and saves it to the file tree.nwk, while saving the parameters to the comma-separated file params.csv:
```bash
generate_bd --min_tips 200 --max_tips 500 \
--la 0.5 --psi 0.25 --p 0.5 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BD-CT(1) model, with λ=0.5, ψ=0.25, p=0.5, φ=2.5, υ=0.2, 
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bd --min_tips 200 --max_tips 500 \
--la 0.5 --psi 0.25 --p 0.5 \
--phi 2.5 --upsilon 0.2 --max_notified_contacts 1 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BD-CT(1)-Skyline model with two time intervals, 
with λ=0.5, ψ=0.25, p=0.5, φ=2.5, υ=0 between t=0 and t=3,
and λ=1, ψ=0.25, p=0.75, φ=2.5, υ=0.2 starting at t=3,
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bd --min_tips 200 --max_tips 500 \
--la 0.5 1 --psi 0.25 0.25 --p 0.5 0.75 \
--phi 2.5 2.5 --upsilon 0 0.2 --max_notified_contacts 1 \
--skyline_times 3 \
--nwk tree.nwk --log params.csv
```
To see detailed options, run:
```bash
generate_bd --help
```

#### BDEI, BDEI-CT(κ) and BDEI-CT(κ)-Skyline
The following command simulates a tree with 200-500 tips under the BDEI model, 
with μ=1, λ=0.5, ψ=0.25, p=0.5, 
and saves it to the file tree.nwk, while saving the parameters to the comma-separated file params.csv:
```bash
generate_bdei --min_tips 200 --max_tips 500 \
--mu 1 --la 0.5 --psi 0.25 --p 0.5 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDEI-CT(2) model, 
with μ=1, λ=0.5, ψ=0.25, p=0.5, φ=2.5, υ=0.2, 
and allowing to notify last two contacts of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdei --min_tips 200 --max_tips 500 \
--mu 1 --la 0.5 --psi 0.25 --p 0.5 \
--phi 2.5 --upsilon 0.2 --max_notified_contacts 2 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDEI-CT(2)-Skyline model with three time intervals, 
with μ=1, λ=0.5, ψ=0.25, p=0.2, φ=2.5, υ=0.2, between t=0 and t=2,
with μ=1, λ=0.5, ψ=0.3, p=0.3, φ=2.5, υ=0.3, between t=2 and t=3,
and μ=1, λ=0.5, ψ=0.5, p=0.5, φ=5, υ=0.3 starting at t=3,
and allowing to notify last two contacts of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdei --min_tips 200 --max_tips 500 \
--mu 1 1 1 --la 0.5 0.5 0.5 --psi 0.25 0.3 0.5 --p 0.2 0.3 0.5 \
--phi 2.5 2.5 5 --upsilon 0.2 0.3 0.3 --max_notified_contacts 2 \
--skyline_times 2 3 \
--nwk tree.nwk --log params.csv
```

To see detailed options, run:
```bash
generate_bdei --help
```


#### BDSS, BDSS-CT(κ) and BDSS-CT(κ)-Skyline
The following command simulates a tree with 200-500 tips under the BDSS model, 
with λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=0.5, λ<sub>ss</sub>=1.5, ψ=0.25, p=0.5, 
and saves it to the file tree.nwk, while saving the parameters to the comma-separated file params.csv:
```bash
generate_bdss --min_tips 200 --max_tips 500 \
--la_nn 0.1 --la_ns 0.3 --la_sn 0.5 --la_ss 1.5 --psi 0.25 --p 0.5 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDSS-CT(3) model, 
with λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=0.5, λ<sub>ss</sub>=1.5, ψ=0.25, p=0.5, φ=2.5, υ=0.2, 
and allowing to notify last three contacts of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdss --min_tips 200 --max_tips 500 \
--la_nn 0.1 --la_ns 0.3 --la_sn 0.5 --la_ss 1.5 --psi 0.25 --p 0.5 \
--phi 2.5 --upsilon 0.2 --max_notified_contacts 3 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDSS-CT(3)-Skyline model with two time intervals, 
with λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=0.5, λ<sub>ss</sub>=1.5, 
ψ=0.25, p=0.5, φ=2.5, υ=0.2 between t=0 and t=2,
and λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=1, λ<sub>ss</sub>=3, 
ψ=0.25, p=0.5, φ=5, υ=0.5 starting at t=2, 
and allowing to notify last three contacts of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdss --min_tips 200 --max_tips 500 \
--la_nn 0.1 0.1 --la_ns 0.3 0.3 --la_sn 0.5 1 --la_ss 1.5 3 --psi 0.25 0.25 --p 0.5 0.5 \
--phi 2.5 5 --upsilon 0.2 0.5 --max_notified_contacts 3 \
--skyline_times 2 \
--nwk tree.nwk --log params.csv
```
To see detailed options, run:
```bash
generate_bdss --help
```

#### User-defined MTBD, MTBD-CT(κ) and MTBD-CT(κ)-Skyline models
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
The following command simulates a tree with 200-500 tips under a generic MTBD-CT(1) model, with two states A and B, 
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
The following command simulates a tree with 200-500 tips under a generic MTBD-CT(1)-Skyline model, 
with two states A and B, 
with μ<sub>aa</sub>=0.5, μ<sub>ab</sub>=0.6, μ<sub>ba</sub>=0.7, μ<sub>bb</sub>=0.8, 
λ<sub>aa</sub>=0.1, λ<sub>ab</sub>=0.2, λ<sub>ba</sub>=0.3, λ<sub>bb</sub>=0.4, 
ψ<sub>a</sub>=0.05, ψ<sub>b</sub>=0.08,
p=<sub>a</sub>0.15, p=<sub>b</sub>0.65,
φ=2.5, υ=0.2  between t=0 and t=8,
and μ<sub>aa</sub>=1.5, μ<sub>ab</sub>=1.6, μ<sub>ba</sub>=1.7, μ<sub>bb</sub>=1.8, 
λ<sub>aa</sub>=1.1, λ<sub>ab</sub>=1.2, λ<sub>ba</sub>=1.3, λ<sub>bb</sub>=1.4, 
ψ<sub>a</sub>=1.05, ψ<sub>b</sub>=1.08,
p=<sub>a</sub>0.1, p=<sub>b</sub>0.6,
φ=3.5, υ=0.4 starting at t=8,
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_mtbd --min_tips 200 --max_tips 500 \
--states A B \
--transition_rates 0.5 0.6 0.7 0.8 1.5 1.6 1.7 1.8\
--transmission_rates 0.1 0.2 0.3 0.4 1.1 1.2 1.3 1.4 \
--removal_rates 0.05 0.08 1.05 1.08 \
--sampling_probabilities 0.15 0.65 0.1 0.6 \
--phi 2.5 3.5 --upsilon 0.2 0.4 --max_notified_contacts 1 \
--skyline_times 8 \
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

# 1. BD, BD-CT(1) and BD-CT(1)-Skyline
## BD model
bd_model = BirthDeathModel(p=0.5, la=0.5, psi=0.25)
print(bd_model.get_epidemiological_parameters())
[bd_tree], _, _ = generate([bd_model], min_tips=200, max_tips=500)
save_forest([bd_tree], 'BD_tree.nwk')
## Adding -CT to the model above
bdct_model = CTModel(model=bd_model, upsilon=0.2, notified_removal_rate=2.5, max_notified_contacts=1)
[bdct_tree], _, _ = generate([bdct_model], min_tips=200, max_tips=500)
save_forest([bdct_tree], 'BDCT_tree.nwk')
## BD-CT(1)-Skyline models
bdct_model_1 = CTModel(BirthDeathModel(p=0.5, la=0.5, psi=0.25), 
                       upsilon=0, notified_removal_rate=2.5, max_notified_contacts=1)
bdct_model_2 = CTModel(BirthDeathModel(p=0.75, la=1, psi=0.25), 
                       upsilon=0.2, notified_removal_rate=2.5, max_notified_contacts=1)
[bdct_skyline_tree], _, _ = generate([bdct_model_1, bdct_model_2], skyline_times=[3], 
                             min_tips=200, max_tips=500)
save_forest([bdct_skyline_tree], 'BDCTSkyline_tree.nwk')

# BDEI, BDEI-CT(2) and BDEI-CT(2)-Skyline
## BDEI model
bdei_model = BirthDeathExposedInfectiousModel(p=0.5, mu=1, la=0.5, psi=0.25)
print(bdei_model.get_epidemiological_parameters())
[bdei_tree], _, _ = generate([bdei_model], min_tips=200, max_tips=500)
save_forest([bdei_tree], 'BDEI_tree.nwk')
## Adding -CT to the model above
bdeict_model = CTModel(model=[bdei_model], upsilon=0.2, notified_removal_rate=2.5, max_notified_contacts=2)
[bdeict_tree], _, _ = generate(bdeict_model, min_tips=200, max_tips=500)
save_forest([bdeict_tree], 'BDEICT_tree.nwk')
## BDEI-CT(2)-Skyline with three time intervals
bdeict_model_1 = CTModel(model=BirthDeathExposedInfectiousModel(p=0.2, mu=1, la=0.5, psi=0.25), upsilon=0.2, notified_removal_rate=2.5, max_notified_contacts=2)
bdeict_model_2 = CTModel(model=BirthDeathExposedInfectiousModel(p=0.3, mu=1, la=0.5, psi=0.3), upsilon=0.3, notified_removal_rate=2.5, max_notified_contacts=2)
bdeict_model_3 = CTModel(model=BirthDeathExposedInfectiousModel(p=0.5, mu=1, la=0.5, psi=0.5), upsilon=0.3, notified_removal_rate=5, max_notified_contacts=2)
[bdeict_skyline_tree], _, _ = generate([bdeict_model_1, bdeict_model_2, bdeict_model_3], skyline_times=[2, 3], min_tips=200, max_tips=500)
save_forest([bdeict_skyline_tree], 'BDEICTSkyline_tree.nwk')

# BDSS, BDSS-CT(3) and BDSS-CT(3)-Skyline
## BDSS model
bdss_model = BirthDeathWithSuperSpreadingModel(p=0.5, la_nn=0.1, la_ns=0.3, la_sn=0.5, la_ss=1.5, psi=0.25)
print(bdss_model.get_epidemiological_parameters())
[bdss_tree], _, _ = generate(bdss_model, min_tips=200, max_tips=500)
save_forest([bdss_tree], 'BDSS_tree.nwk')
## Adding -CT to the model above
bdssct_model = CTModel(model=bdss_model, upsilon=0.2, notified_removal_rate=2.5, max_notified_contacts=3)
[bdssct_tree], _, _ = generate(bdssct_model, min_tips=200, max_tips=500)
save_forest([bdssct_tree], 'BDSSCT_tree.nwk')
## BDSS-CT(3)-Skyline with two time intervals, using the model above for the first interval
bdssct_model_2 = CTModel(model=BirthDeathWithSuperSpreadingModel(p=0.5, la_nn=0.1, la_ns=0.3, la_sn=1, la_ss=3, psi=0.25), 
                         upsilon=0.5, notified_removal_rate=5, max_notified_contacts=2)
[bdssct_skyline_tree], _, _ = generate([bdssct_model, bdssct_model_2], skyline_times=[2], min_tips=200, max_tips=500)
save_forest([bdssct_skyline_tree], 'BDSSCTSkyline_tree.nwk')

# MTBD, MTBD-CT(1) and MTBD-CT(1)-Skyline
## MTBD model with two states: A and B
mtbd_model = Model(states=['A', 'B'], transition_rates=[[0.5, 0.6], [0.7, 0.8]],
                   transmission_rates=[[0.1, 0.2], [0.3, 0.4]],
                   removal_rates=[0.05, 0.08], ps=[0.15, 0.65])
[mtbd_tree], _, _ = generate(mtbd_model, min_tips=200, max_tips=500)
save_forest([mtbd_tree], 'MTBD_tree.nwk')
## Adding -CT to the model above
mtbdct_model = CTModel(model=mtbd_model, upsilon=0.2, notified_removal_rate=2.5, max_notified_contacts=1)
[mtbdct_tree], _, _ = generate(mtbdct_model, min_tips=200, max_tips=500)
save_forest([mtbdct_tree], 'MTBDCT_tree.nwk')
## MTBD-CT(1)-Skyline with two time intervals, using the model above for the first interval
mtbdct_model_2 = CTModel(model=Model(states=['A', 'B'], transition_rates=[[1.5, 1.6], [1.7, 1.8]],
                                     transmission_rates=[[1.1, 1.2], [1.3, 1.4]],
                                     removal_rates=[1.05, 1.08], ps=[0.1, 0.6]), 
                         upsilon=0.4, notified_removal_rate=3.5, max_notified_contacts=1)
[mtbdct_skyline_tree], _, _ = generate([mtbdct_model, mtbdct_model_2], skyline_times=[8], min_tips=200, max_tips=500)
save_forest([mtbdct_skyline_tree], 'MTBDCTSkyline_tree.nwk')
```

### Run with apptainer

Once [apptainer](https://apptainer.org/docs/user/latest/quick_start.html#installation) is installed, 
run the following command:

```bash
apptainer run docker://evolbioinfo/treesimlator
```

This will launch a terminal session within the container, 
in which you can run treesimulator following the instructions for the command line ("Basic usage in a command line") above.





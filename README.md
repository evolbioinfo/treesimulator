# treesimulator

Simulation of rooted phylogenetic trees under a given Multi-Type Birth–Death (MTBD) model, 
with or without contact tracing (CT), 
and with or without Skyline. 

## Article
treesimulator is described in the [S1 Appendix](https://doi.org/10.1371/journal.pcbi.1012461.s001) of the following article:

Anna Zhukova, Olivier Gascuel (2025) Accounting for contact tracing in epidemiological birth-death models. _PLOS Comput Biol_ 21(5): e1012461; doi:[10.1371/journal.pcbi.1012461](https://doi.org/10.1371/journal.pcbi.1012461)

[![DOI:10.1371/journal.pcbi.1012461](https://zenodo.org/badge/DOI/10.1371/journal.pcbi.1012461.svg)](https://doi.org/10.1371/journal.pcbi.1012461)
[![GitHub release](https://img.shields.io/github/v/release/evolbioinfo/treesimulator.svg)](https://github.com/evolbioinfo/treesimulator/releases)
[![PyPI version](https://badge.fury.io/py/treesimulator.svg)](https://pypi.org/project/treesimulator/)
[![PyPI downloads](https://shields.io/pypi/dm/treesimulator)](https://pypi.org/project/treesimulator/)
[![Docker pulls](https://img.shields.io/docker/pulls/evolbioinfo/treesimulator)](https://hub.docker.com/r/evolbioinfo/treesimulator/tags)

## MTBD
The MTBD models were introduced by Stadler & Bonhoeffer [[_Philos. Trans. R. Soc. B_ 2013]](https://royalsocietypublishing.org/doi/10.1098/rstb.2012.0198).

An MTBD model with m states has 

m(m-1) state transition rate parameters:
* μ<sub>ij</sub> -- transition rate from state i to state j (1 ≤ i, j ≤ m; i ≠ j), where μ<sub>ij</sub> ≥ 0
  (In practice, we ask the user to provide an m x m matrix for MTBD transition rates. μ<sub>ii</sub> must be 0 for all i)

m<sup>2</sup> transmission rate parameters:
* λ<sub>ij</sub> -- transmission rate from state i (donor) to state j (recipient) (1 ≤ i, j ≤ m), where λ<sub>ij</sub> ≥ 0

m removal (becoming non-infectious) rate parameters:
* ψ<sub>i</sub> -- removal rate of state i (1 ≤ i ≤ m), where ψ<sub>i</sub> ≥ 0

m sampling probability upon removal parameters:
* p<sub>i</sub> -- probability to sample the pathogen of an individual in state i upon removal (1 ≤ i ≤ m), where 0 < p<sub>i</sub> ≤ 1

Using these probabilities, one can calculate the equilibrium frequencies π<sub>i</sub> of the model's states, 
where 0 ≤ π<sub>i</sub> ≤ 1 and π<sub>1</sub> + ... + π<sub>m</sub> = 1.

One can also calculate the following rates, which are useful in epidemiological parameter calculations:
* λ<sub>i</sub> = Σ<sub>1≤j≤m</sub>λ<sub>ij</sub> -- total transmission rate from state i
* μ<sub>i</sub> = Σ<sub>1≤j≤m;i≠j</sub>μ<sub>ij</sub> -- total state-change rate from state i
* ε<sub>i</sub> = μ<sub>i</sub> + ψ<sub>i</sub> -- state i exit rate

The MTBD model has the following epidemiological parameters:

* R<sub>i</sub> = λ<sub>i</sub>/ε<sub>i</sub> + Σ<sub>1≤j≤m</sub>(μ<sub>ij</sub>/ε<sub>i</sub>)R<sub>j</sub>-- reproduction number of state i
* t<sub>i</sub> = 1/ε<sub>i</sub> -- exit time for state i
* d<sub>i</sub> = 1/ε<sub>i</sub> + Σ<sub>1≤j≤m</sub>(μ<sub>ij</sub>/ε<sub>i</sub>)d<sub>j</sub> -- time till the end of infection for state i
* R = Σ<sub>1≤i≤m</sub>π<sub>i</sub>R<sub>i</sub> -- average reproduction number
* d = Σ<sub>i ∈ recipient states</sub>(π<sub>i</sub>d<sub>i</sub>)/Σ<sub>j ∈ recipient states</sub>π<sub>j</sub> -- average infection time


## Contact Tracing (CT)
Contact tracing extension was introduced by Zhukova & Gascuel [[_PLOS Comput Biol_ 2025]](https://doi.org/10.1371/journal.pcbi.1012461). 
Here we present a slightly modified version of it. CT adds three parameters to the initial MTBD model:

* υ -- probability to notify contacts upon sampling (default) or becoming non-infectious 
    (a.k.a. removal, if the option _--notify_at_removal_ is set).
* X<sub>C</sub> -- notified removal speeed-up: X<sub>C</sub> >> 1. 
Once notified, the contact's removal rate increases X<sub>C</sub> times. 
* X<sub>p</sub> -- notified sampling probability speeed-up: X<sub>p</sub> >> 1. 
Once notified and removed, the contact's pathogen's sampling probability increases up to X<sub>p</sub> times (capped at 1). 
X<sub>p</sub> = 1 corresponds to no increase, 
while X<sub>p</sub> ≥ 1/p (where p is the original sampling probability) corresponds to automatic sampling upon removal.

CT also adds a meta-parameter κ, which defines how many most recent contacts can be notified by each index case. Each is notified independently, with the probability υ.

CT extension adds m notified contact states (1<sub>C</sub>, ..., m<sub>C</sub>) to its original MTBD model, 
where i<sub>C</sub> is a notified version of the state i.

Transition rates for a notified state i<sub>C</sub> are analogous to those of i (μ<sub>i<sub>C</sub>j<sub>C</sub></sub>=μ<sub>ij</sub>), 
while transitions from non-notified states to notified ones and vice versa are not allowed (μ<sub>i<sub>C</sub>j</sub>=μ<sub>ji<sub>C</sub></sub>=0).

Transmission rates for a notified state i<sub>C</sub> are the same as those of i (λ<sub>i<sub>C</sub>j</sub>=λ<sub>ij</sub>), 
where the recipients are always in a non-notified state (λ<sub>i<sub>C</sub>j<sub>C</sub></sub>=λ<sub>ij<sub>C</sub></sub>=0).

The removal rates for notified states are X<sub>C</sub> times higher than those of the corresponding non-notified states:
ψ<sub>i<sub>C</sub></sub> = X<sub>C</sub>ψ<sub>i</sub>.

The sampling probabilities for notified states are X<sub>p</sub> times higher (but are ≤ 1) than those of the corresponding non-notified states:
p<sub>i<sub>C</sub></sub> = X<sub>p</sub>p<sub>i</sub>.

## Skyline
Skyline was introduced by Stadler _et al._ [[PNAS 2013]](https://doi.org/10.1073/pnas.1207965110) 
and extended to MTBD by Kühnert _et al._ [[MBE 2016]](https://doi.org/10.1093/molbev/msw064). 
It enables piece-wise constant parameter value changes. 
To use a skyline with k models, one needs to specify k-1 model change times T<sub>1</sub>, ..., T<sub>k-1</sub>, 
and k sets of model parameters (see above).
At time 0 the simulation starts with model 1, it switches to models 2 at time T<sub>1</sub>, etc.
All the models in the Skyline must have the same states. 
CT-related parameters can also change at skyline changing times, 
in that case if some skyline intervals do not have CT, υ=0 and any values for X<sub>C</sub> and X<sub>p</sub> must be specified for them. 
The same κ value is shared among all the skyline -CT models.


We pay particular interest to the classical BD model, the BD Exposed-Infectious (BDEI) model, 
and BD with super-spreading (BDSS), 
as they are described in [[Voznica _et al._ 2021]](https://www.biorxiv.org/content/10.1101/2021.03.11.435006v1), 
and to their -CT(κ) versions.


## BD
1 state: 
* I (infectious)

3 parameters:
* λ = λ<sub>I</sub>  -- transmission rate
* ψ = ψ<sub>I</sub> -- removal rate
* p = p<sub>I</sub> -- sampling probability upon removal

Epidemiological parameters:
* R = R<sub>I</sub> = λ/ψ -- reproduction number
* d<sub>I</sub> = 1/ψ -- infectious time 


## BD-CT(κ)
2 states: 
* I, infectious
* I<sub>C</sub>, notified infectious contact

6 parameters:
* λ = λ<sub>I</sub> = λ<sub>I<sub>C</sub></sub> -- transmission rate
* ψ = ψ<sub>I</sub> -- removal rate of I
* p = p<sub>I</sub> -- sampling probability upon removal of I
* υ -- probability to notify contacts upon sampling
* X<sub>C</sub> -- notified removal speed-up: X<sub>C</sub> >> 1
* X<sub>p</sub> -- notified sampling probability speed-up: X<sub>p</sub> ≥ 1


## BDEI
2 states: 
* E, exposed, i.e. infected but not yet infectious
* I, infectious

4 parameters:
* μ = μ<sub>EI</sub> -- transition rate from E to I (becoming infectious)
* λ = λ<sub>IE</sub> -- transmission rate from I to E
* ψ = ψ<sub>I</sub> -- removal rate of I
* p = p<sub>I</sub> -- sampling probability upon removal of I

Epidemiological parameters:
* R = R<sub>E</sub> = R<sub>I</sub> = λ/ψ -- reproduction number
* t<sub>I</sub> = 1/ψ -- infectious time 
* t<sub>E</sub> = 1/μ -- incubation period
* d = d<sub>E</sub> = 1/μ + 1/ψ -- total infection time
* f<sub>E</sub> = t<sub>E</sub>/d -- incubation fraction


## BDEI-CT(κ)
4 states: 
* E, exposed, i.e. infected but not yet infectious
* I, infectious
* E<sub>C</sub>, notified exposed contact
* I<sub>C</sub>, notified infectious contact

7 parameters:
* μ = μ<sub>EI</sub> = μ<sub>E<sub>C</sub>I<sub>C</sub></sub>-- transition rate from an exposed state (notified or not) 
to the corresponding infectious state (becoming infectious)
* λ = λ<sub>IE</sub> = λ<sub>I<sub>C</sub>E</sub> -- transmission rate from an infectious state (notified or not) to E
* ψ = ψ<sub>I</sub> -- removal rate of I
* p = p<sub>I</sub> -- sampling probability upon removal
* υ -- probability to notify contacts upon sampling
* X<sub>C</sub> -- notified removal speed-up: X<sub>C</sub> >> 1
* X<sub>p</sub> -- notified sampling probability speed-up: X<sub>p</sub> ≥ 1

## BDSS
2 states: 
* I, standard infectious individual (a.k.a. normal spreader)
* S, superspreader

5(+1) parameters:
* λ<sub>nn</sub> = λ<sub>II</sub> -- transmission rate from I to I
* λ<sub>ns</sub> = λ<sub>IS</sub> -- transmission rate from I to S
* λ<sub>sn</sub> = λ<sub>SI</sub> -- transmission rate from S to I
* λ<sub>ss</sub> = λ<sub>SS</sub> -- transmission rate from S to S

    (with a constraint that λ<sub>ss</sub>/λ<sub>ns</sub>=λ<sub>sn</sub>/λ<sub>nn</sub>)
* ψ = ψ<sub>I</sub> = ψ<sub>S</sub> -- removal rate
* p = p<sub>I</sub> = p<sub>S</sub> -- sampling probability upon removal

Epidemiological parameters:
* X<sub>S</sub>=λ<sub>ss</sub>/λ<sub>ns</sub>=λ<sub>sn</sub>/λ<sub>nn</sub> -- super-spreading transmission ratio
* f<sub>S</sub>=λ<sub>ss</sub>/(λ<sub>sn</sub> + λ<sub>ss</sub>) -- super-spreading fraction
* R = (1 - f<sub>S</sub>) λ/ψ + f<sub>S</sub> X<sub>S</sub>λ/ψ -- reproduction number
* d = 1/ψ -- total infection time

## BDSS-CT(κ)
4 states: 
* I, standard infectious individual (a.k.a. normal spreader)
* S, superspreader
* I<sub>C</sub>, notified normal spreader
* S<sub>C</sub>, notified superspreader

8(+1) parameters:
* λ<sub>nn</sub> = λ<sub>II</sub> = λ<sub>I<sub>C</sub>I</sub> -- transmission rate from a normal spreader (notified or not) to I
* λ<sub>ns</sub> = λ<sub>IS</sub> = λ<sub>I<sub>C</sub>S</sub> -- transmission rate from a normal spreader (notified or not) to S
* λ<sub>sn</sub> = λ<sub>SI</sub> = λ<sub>S<sub>C</sub>I</sub> -- transmission rate from a superspreader (notified or not) to I
* λ<sub>ss</sub> = λ<sub>SS</sub> = λ<sub>S<sub>C</sub>S</sub> -- transmission rate from a superspreader (notified or not) to S

    (with a constraint that λ<sub>ss</sub>/λ<sub>ns</sub>=λ<sub>sn</sub>/λ<sub>nn</sub>)
* ψ = ψ<sub>I</sub> = ψ<sub>S</sub> -- removal rate of a non-notified state
* p = p<sub>I</sub> = p<sub>S</sub> -- sampling probability upon removal of a non-notified state
* υ -- probability to notify contacts upon sampling
* X<sub>C</sub> -- notified removal speed-up: X<sub>C</sub> >> 1
* X<sub>p</sub> -- notified sampling probability speed-up: X<sub>p</sub> ≥ 1


## BDEISS
3 states: 
* E, exposed, i.e. infected but not yet infectious
* I, standard infectious individual (a.k.a. normal spreader)
* S, infectious superspreader 

6 parameters:
* μ<sub>n</sub> = μ<sub>EI</sub> -- transition rate from E to I (becoming infectious for normal spreaders)
* μ<sub>s</sub> = μ<sub>ES</sub> -- transition rate from E to S (becoming infectious for superspreaders)
* λ<sub>n</sub> = λ<sub>IE</sub> -- transmission rate from I to E
* λ<sub>s</sub> = λ<sub>SE</sub> -- transmission rate from S to E
* ψ = ψ<sub>I</sub> = ψ<sub>S</sub> -- removal rate of I and of S (the same)
* p = p<sub>I</sub> = p<sub>S</sub> -- sampling probability upon removal (the same for I and S)

Epidemiological parameters:
* X<sub>S</sub> = λ<sub>s</sub>/λ<sub>n</sub> -- super-spreading transmission ratio
* f<sub>S</sub> = μ<sub>s</sub>/(μ<sub>n</sub> + μ<sub>s</sub>) -- super-spreading fraction (among infectious individuals)
* t<sub>E</sub> = 1/(μ<sub>n</sub> + μ<sub>s</sub>) -- incubation period
* t<sub>I</sub> = t<sub>S</sub> = 1/ψ -- infectious time
* R = (1 - f<sub>S</sub>) λ/ψ + f<sub>S</sub> X<sub>S</sub>λ/ψ -- reproduction number
* d = d<sub>E</sub> = 1/(μ<sub>n</sub> + μ<sub>s</sub>) + 1/ψ -- total infection time
* f<sub>E</sub> = t<sub>E</sub>/d -- incubation fraction


## BDEISS-CT(κ)
6 states: 
* E, exposed, i.e. infected but not yet infectious
* N, standard infectious individual (a.k.a. normal spreader)
* S, superspreader
* E<sub>C</sub>, notified exposed individual
* I<sub>C</sub>, notified normal spreader
* S<sub>C</sub>, notified superspreader

9 parameters:
* μ<sub>n</sub> = μ<sub>EI</sub> = μ<sub>E<sub>C</sub>I<sub>C</sub></sub> -- transition rate from an exposed state (notified or not) to the corresponding normal infectious state (becoming infectious for normal spreaders)
* μ<sub>s</sub> = μ<sub>ES</sub> = μ<sub>E<sub>C</sub>S<sub>C</sub></sub> -- transition rate from an exposed state (notified or not) to the corresponding superspreader infectious state (becoming infectious for superspreaders)
* λ<sub>n</sub> = λ<sub>IE</sub> = λ<sub>I<sub>C</sub>E</sub> -- transmission rate from I or I<sub>C</sub> to E
* λ<sub>s</sub> = λ<sub>SE</sub> = λ<sub>S<sub>C</sub>E</sub> -- transmission rate from S or S<sub>C</sub> to E
* ψ = ψ<sub>I</sub> = ψ<sub>S</sub> -- removal rate of I and of S (the same)
* p = p<sub>I</sub> = p<sub>S</sub> -- sampling probability upon removal of I or S
* υ -- probability to notify contacts upon sampling
* X<sub>C</sub> -- notified removal speed-up: X<sub>C</sub> >> 1
* X<sub>p</sub> -- notified sampling probability speed-up: X<sub>p</sub> ≥ 1


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
The following command simulates a tree with 200-500 tips under the BD-CT(1) model, with λ=0.5, ψ=0.25, p=0.5, 
υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1, 
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bd --min_tips 200 --max_tips 500 \
--la 0.5 --psi 0.25 --p 0.5 \
--X_C 10 --upsilon 0.2 --X_p 1 --max_notified_contacts 1 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BD-CT(1)-Skyline model with two time intervals, 
with λ=0.5, ψ=0.25, p=0.5, υ=0, X<sub>C</sub>=10, X<sub>p</sub>=1 between t=0 and t=3,
and λ=1, ψ=0.25, p=0.75, υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1 starting at t=3,
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bd --min_tips 200 --max_tips 500 \
--la 0.5 1 --psi 0.25 0.25 --p 0.5 0.75 \
--X_C 10 10 --upsilon 0 0.2 --X_p 1 1 --max_notified_contacts 1 \
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
with μ=1, λ=0.5, ψ=0.25, p=0.5, υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1, 
and allowing to notify last two contacts of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdei --min_tips 200 --max_tips 500 \
--mu 1 --la 0.5 --psi 0.25 --p 0.5 \
--X_C 10 --upsilon 0.2 --X_p 1 --max_notified_contacts 2 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDEI-CT(2)-Skyline model with three time intervals, 
with μ=1, λ=0.5, ψ=0.25, p=0.2, υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1, between t=0 and t=2,
with μ=1, λ=0.5, ψ=0.3, p=0.3, υ=0.3, X<sub>C</sub>=10, X<sub>p</sub>=1, between t=2 and t=3,
and μ=1, λ=0.5, ψ=0.5, p=0.5, υ=0.3, X<sub>C</sub>=20, X<sub>p</sub>=1 starting at t=3,
and allowing to notify last two contacts of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdei --min_tips 200 --max_tips 500 \
--mu 1 1 1 --la 0.5 0.5 0.5 --psi 0.25 0.3 0.5 --p 0.2 0.3 0.5 \
--X_C 10 10 20 --upsilon 0.2 0.3 0.3 --X_p 1 1 1 --max_notified_contacts 2 \
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
with λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=0.5, λ<sub>ss</sub>=1.5, ψ=0.25, p=0.5, 
υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1, 
and allowing to notify last three contacts of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdss --min_tips 200 --max_tips 500 \
--la_nn 0.1 --la_ns 0.3 --la_sn 0.5 --la_ss 1.5 --psi 0.25 --p 0.5 \
--X_C 10 --upsilon 0.2 --X_p 1 --max_notified_contacts 3 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDSS-CT(3)-Skyline model with two time intervals, 
with λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=0.5, λ<sub>ss</sub>=1.5, 
ψ=0.25, p=0.5, υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1 between t=0 and t=2,
and λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=1, λ<sub>ss</sub>=3, 
ψ=0.25, p=0.5, υ=0.5, X<sub>C</sub>=20, X<sub>p</sub>=1 starting at t=2, 
and allowing to notify last three contacts of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdss --min_tips 200 --max_tips 500 \
--la_nn 0.1 0.1 --la_ns 0.3 0.3 --la_sn 0.5 1 --la_ss 1.5 3 --psi 0.25 0.25 --p 0.5 0.5 \
--X_C 10 20 --upsilon 0.2 0.5 --X_p 1 1 --max_notified_contacts 3 \
--skyline_times 2 \
--nwk tree.nwk --log params.csv
```
To see detailed options, run:
```bash
generate_bdss --help
```


#### BDEISS, BDEISS-CT(κ) and BDEISS-CT(κ)-Skyline
The following command simulates a tree with 200-500 tips under the BDEISS model, 
with μ<sub>n</sub>=0.1, μ<sub>s</sub>=0.3, λ<sub>n</sub>=0.5, λ<sub>s</sub>=1.5, ψ=0.25, p=0.5, 
and saves it to the file tree.nwk, while saving the parameters to the comma-separated file params.csv:
```bash
generate_bdeiss --min_tips 200 --max_tips 500 \
--mu_n 0.1 --mu_s 0.3 --la_n 0.5 --la_s 1.5 --psi 0.25 --p 0.5 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDEISS-CT(1) model, 
with μ<sub>n</sub>=0.1, μ<sub>s</sub>=0.3, λ<sub>n</sub>=0.5, λ<sub>s</sub>=1.5, ψ=0.25, p=0.5, 
υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1, 
and allowing to notify the last contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdeiss --min_tips 200 --max_tips 500 \
--mu_n 0.1 --mu_s 0.3 --la_n 0.5 --la_s 1.5 --psi 0.25 --p 0.5 \
--X_C 10 --upsilon 0.2 --X_p 1 --max_notified_contacts 1 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under the BDSS-CT(1)-Skyline model with two time intervals, 
with μ<sub>n</sub>=0.1, μ<sub>s</sub>=0.3, λ<sub>n</sub>=0.5, λ<sub>s</sub>=1.5, ψ=0.25, p=0.5, 
υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1
between t=0 and t=2,
and μ<sub>n</sub>=0.1, μ<sub>s</sub>=0.3, λ<sub>n</sub>=0.5, λ<sub>s</sub>=1.5, ψ=0.25, p=0.2, 
υ=0.5, X<sub>C</sub>=10, X<sub>p</sub>=1
starting at t=2, 
and allowing to notify the last contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_bdeiss --min_tips 200 --max_tips 500 \
--mu_n 0.1 0.1 --mu_s 0.3 0.3 --la_n 0.5 0.5 --la_s 1.5 1.5 --psi 0.25 0.25 --p 0.5 0.2 \
--X_C 10 10 --upsilon 0.2 0.5 --X_p 1 1 --max_notified_contacts 1 \
--skyline_times 2 \
--nwk tree.nwk --log params.csv
```
To see detailed options, run:
```bash
generate_bdss --help
```

#### User-defined MTBD, MTBD-CT(κ) and MTBD-CT(κ)-Skyline models
The following command simulates a tree with 200-500 tips under a generic MTBD model, with two states A and B, 
with μ<sub>ab</sub>=0.6, μ<sub>ba</sub>=0.7 
(note that μ<sub>aa</sub>=μ<sub>bb</sub>=0 as only transitions between different states are possible!), 
λ<sub>aa</sub>=0.1, λ<sub>ab</sub>=0.2, λ<sub>ba</sub>=0.3, λ<sub>bb</sub>=0.4, 
ψ<sub>a</sub>=0.05, ψ<sub>b</sub>=0.08,
p=<sub>a</sub>0.15, p=<sub>b</sub>0.65,
and saves it to the file tree.nwk, while saving the parameters to the comma-separated file params.csv:
```bash
generate_mtbd --min_tips 200 --max_tips 500 \
--states A B \
--transition_rates 0 0.6 0.7 0 \
--transmission_rates 0.1 0.2 0.3 0.4 \
--removal_rates 0.05 0.08 \
--sampling_probabilities 0.15 0.65 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under a generic MTBD-CT(1) model, with two states A and B, 
with μ<sub>ab</sub>=0.6, μ<sub>ba</sub>=0.7, 
λ<sub>aa</sub>=0.1, λ<sub>ab</sub>=0.2, λ<sub>ba</sub>=0.3, λ<sub>bb</sub>=0.4, 
ψ<sub>a</sub>=0.05, ψ<sub>b</sub>=0.08,
p=<sub>a</sub>0.15, p=<sub>b</sub>0.65,
υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1, 
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_mtbd --min_tips 200 --max_tips 500 \
--states A B \
--transition_rates 0 0.6 0.7 0 \
--transmission_rates 0.1 0.2 0.3 0.4 \
--removal_rates 0.05 0.08 \
--sampling_probabilities 0.15 0.65 \
--X_C 10 --upsilon 0.2 --X_p 1 --max_notified_contacts 1 \
--nwk tree.nwk --log params.csv
```
The following command simulates a tree with 200-500 tips under a generic MTBD-CT(1)-Skyline model, 
with two states A and B, 
with μ<sub>ab</sub>=0.6, μ<sub>ba</sub>=0.7, 
λ<sub>aa</sub>=0.1, λ<sub>ab</sub>=0.2, λ<sub>ba</sub>=0.3, λ<sub>bb</sub>=0.4, 
ψ<sub>a</sub>=0.05, ψ<sub>b</sub>=0.08,
p=<sub>a</sub>0.15, p=<sub>b</sub>0.65,
υ=0.2, X<sub>C</sub>=10, X<sub>p</sub>=1  between t=0 and t=8,
and μ<sub>ab</sub>=1.6, μ<sub>ba</sub>=1.7, 
λ<sub>aa</sub>=1.1, λ<sub>ab</sub>=1.2, λ<sub>ba</sub>=1.3, λ<sub>bb</sub>=1.4, 
ψ<sub>a</sub>=1.05, ψ<sub>b</sub>=1.08,
p=<sub>a</sub>0.1, p=<sub>b</sub>0.6,
υ=0.4, X<sub>C</sub>=20, X<sub>p</sub>=1 starting at t=8,
and allowing to notify only the most recent contact of each sampled index case. 
The simulated tree is saved to the file tree.nwk, while the model parameters are saved to the comma-separated file params.csv:
```bash
generate_mtbd --min_tips 200 --max_tips 500 \
--states A B \
--transition_rates 0 0.6 0.7 0 0 1.6 1.7 0 \
--transmission_rates 0.1 0.2 0.3 0.4 1.1 1.2 1.3 1.4 \
--removal_rates 0.05 0.08 1.05 1.08 \
--sampling_probabilities 0.15 0.65 0.1 0.6 \
--X_C 10 20 --upsilon 0.2 0.4 --X_p 1 1 --max_notified_contacts 1 \
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
from treesimulator import save_forest
from treesimulator.generator import generate
from treesimulator.mtbd_models import Model, BirthDeathModel, BirthDeathExposedInfectiousModel, \
    BirthDeathWithSuperSpreadingModel, BirthDeathExposedInfectiousWithSuperSpreadingModel, CTModel

# 1. BD, BD-CT(1) and BD-CT(1)-Skyline
## BD model
bd_model = BirthDeathModel(p=0.5, la=0.5, psi=0.25)
bd_tree = generate([bd_model], min_tips=200, max_tips=500).sampled_forest[0]
save_forest([bd_tree], 'BD_tree.nwk')
## Adding -CT to the model above
bdct_model = CTModel(model=bd_model, upsilon=0.2, X_C=10, X_p=1)
bdct_tree = generate([bdct_model], min_tips=200, max_tips=500, max_notified_contacts=1).sampled_forest[0]
save_forest([bdct_tree], 'BDCT_tree.nwk')
## BD-CT(1)-Skyline models
bdct_model_1 = CTModel(BirthDeathModel(p=0.5, la=0.5, psi=0.25),
                       upsilon=0, X_C=10, X_p=1)
bdct_model_2 = CTModel(BirthDeathModel(p=0.75, la=1, psi=0.25),
                       upsilon=0.2, X_C=10, X_p=1)
bdct_skyline_tree = generate([bdct_model_1, bdct_model_2], skyline_times=[3],
                             min_tips=200, max_tips=500, max_notified_contacts=1).sampled_forest[0]
save_forest([bdct_skyline_tree], 'BDCTSkyline_tree.nwk')

# BDEI, BDEI-CT(2) and BDEI-CT(2)-Skyline
## BDEI model
bdei_model = BirthDeathExposedInfectiousModel(p=0.5, mu=1, la=0.5, psi=0.25)
print(bdei_model.get_epidemiological_parameters())
bdei_tree = generate([bdei_model], min_tips=200, max_tips=500).sampled_forest[0]
save_forest([bdei_tree], 'BDEI_tree.nwk')
## Adding -CT to the model above
bdeict_model = CTModel(model=bdei_model, upsilon=0.2, X_C=10, X_p=1)
bdeict_tree = generate([bdeict_model], min_tips=200, max_tips=500, max_notified_contacts=2).sampled_forest[0]
save_forest([bdeict_tree], 'BDEICT_tree.nwk')
## BDEI-CT(2)-Skyline with three time intervals
bdeict_model_1 = CTModel(model=BirthDeathExposedInfectiousModel(p=0.2, mu=1, la=0.5, psi=0.25),
                         upsilon=0.2, X_C=10, X_p=1)
bdeict_model_2 = CTModel(model=BirthDeathExposedInfectiousModel(p=0.3, mu=1, la=0.5, psi=0.3),
                         upsilon=0.3, X_C=10, X_p=1)
bdeict_model_3 = CTModel(model=BirthDeathExposedInfectiousModel(p=0.5, mu=1, la=0.5, psi=0.5),
                         upsilon=0.3, X_C=20, X_p=1)
bdeict_skyline_tree = generate([bdeict_model_1, bdeict_model_2, bdeict_model_3], skyline_times=[2, 3],
                               min_tips=200, max_tips=500, max_notified_contacts=2).sampled_forest[0]
save_forest([bdeict_skyline_tree], 'BDEICTSkyline_tree.nwk')

# BDSS, BDSS-CT(3) and BDSS-CT(3)-Skyline
## BDSS model
bdss_model = BirthDeathWithSuperSpreadingModel(p=0.5, la_nn=0.1, la_ns=0.3, la_sn=0.5, la_ss=1.5, psi=0.25)
print(bdss_model.get_epidemiological_parameters())
bdss_tree = generate([bdss_model], min_tips=200, max_tips=500).sampled_forest[0]
save_forest([bdss_tree], 'BDSS_tree.nwk')
## Adding -CT to the model above
bdssct_model = CTModel(model=bdss_model, upsilon=0.2, X_C=10, X_p=1)
bdssct_tree = generate([bdssct_model], min_tips=200, max_tips=500, max_notified_contacts=3).sampled_forest[0]
save_forest([bdssct_tree], 'BDSSCT_tree.nwk')
## BDSS-CT(3)-Skyline with two time intervals, using the model above for the first interval
bdssct_model_2 = CTModel(
    model=BirthDeathWithSuperSpreadingModel(p=0.5, la_nn=0.1, la_ns=0.3, la_sn=1, la_ss=3, psi=0.25),
    upsilon=0.5, X_C=20, X_p=1)
bdssct_skyline_tree = generate([bdssct_model, bdssct_model_2], skyline_times=[2], min_tips=200, max_tips=500,
                               max_notified_contacts=3).sampled_forest[0]
save_forest([bdssct_skyline_tree], 'BDSSCTSkyline_tree.nwk')

# BDEISS, BDEISS-CT(1) and BDEISS-CT(1)-Skyline
## BDEISS model
bdeiss_model = BirthDeathExposedInfectiousWithSuperSpreadingModel(p=0.5, mu_n=0.1, mu_s=0.3, la_n=0.5, la_s=1.5,
                                                                  psi=0.25)
print(bdeiss_model.get_epidemiological_parameters())
bdeiss_tree = generate([bdeiss_model], min_tips=200, max_tips=500).sampled_forest[0]
save_forest([bdeiss_tree], 'BDEISS_tree.nwk')
## Adding -CT to the model above
bdeissct_model = CTModel(model=bdeiss_model, upsilon=0.2, X_C=10, X_p=1)
bdeissct_tree = generate([bdeissct_model], min_tips=200, max_tips=500, max_notified_contacts=1).sampled_forest[0]
save_forest([bdeissct_tree], 'BDEISSCT_tree.nwk')
## BDEISS-CT(1)-Skyline with two time intervals, using the model above for the first interval
bdeissct_model_2 = CTModel(
    model=BirthDeathExposedInfectiousWithSuperSpreadingModel(p=0.2, mu_n=0.1, mu_s=0.3, la_n=0.5, la_s=1.5, psi=0.25),
    upsilon=0.5, X_C=20, X_p=1)
bdeissct_skyline_tree = generate([bdeissct_model, bdeissct_model_2], skyline_times=[2], min_tips=200,
                                 max_tips=500, max_notified_contacts=1).sampled_forest[0]
save_forest([bdeissct_skyline_tree], 'BDEISSCTSkyline_tree.nwk')

# MTBD, MTBD-CT(1) and MTBD-CT(1)-Skyline
## MTBD model with two states: A and B
mtbd_model = Model(states=['A', 'B'], transition_rates=[[0, 0.6], [0.7, 0]],
                   transmission_rates=[[0.1, 0.2], [0.3, 0.4]],
                   removal_rates=[0.05, 0.08], ps=[0.15, 0.65])
mtbd_tree = generate([mtbd_model], min_tips=200, max_tips=500).sampled_forest[0]
save_forest([mtbd_tree], 'MTBD_tree.nwk')
## Adding -CT to the model above
mtbdct_model = CTModel(model=mtbd_model, upsilon=0.2, X_C=10, X_p=1)
mtbdct_tree = generate([mtbdct_model], min_tips=200, max_tips=500, max_notified_contacts=1).sampled_forest[0]
save_forest([mtbdct_tree], 'MTBDCT_tree.nwk')
## MTBD-CT(1)-Skyline with two time intervals, using the model above for the first interval
mtbdct_model_2 = CTModel(model=Model(states=['A', 'B'], transition_rates=[[0, 1.6], [1.7, 0]],
                                     transmission_rates=[[1.1, 1.2], [1.3, 1.4]],
                                     removal_rates=[1.05, 1.08], ps=[0.1, 0.6]),
                         upsilon=0.4, X_C=20, X_p=1)
mtbdct_skyline_tree = generate([mtbdct_model, mtbdct_model_2], skyline_times=[8],
                               min_tips=200, max_tips=500, max_notified_contacts=1).sampled_forest[0]
save_forest([mtbdct_skyline_tree], 'MTBDCTSkyline_tree.nwk')

```

### Run with apptainer

Once [apptainer](https://apptainer.org/docs/user/latest/quick_start.html#installation) is installed, 
run the following command:

```bash
apptainer run docker://evolbioinfo/treesimulator
```

This will launch a terminal session within the container, 
in which you can run treesimulator following the instructions for the command line ("Basic usage in a command line") above.





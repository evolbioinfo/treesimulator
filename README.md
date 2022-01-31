# treesimulator

Simulation of rooted phylogenetic trees under a given Multitype Birth–Death (MTBD) model.
The MTBD models were introduced by Stadler & Bonhoeffer [[_Philos. Trans. R. Soc. B_ 2013]](https://royalsocietypublishing.org/doi/10.1098/rstb.2012.0198).

We pay particular interest to the classical BD model, the BD Exposed-Infectious (BDEI) model, 
and BD with superspreading (BDSS), 
as they are described in [[Voznica _et al._ 2021]]((https://www.biorxiv.org/content/10.1101/2021.03.11.435006v1)). 

## BD
3 parameters:
* λ -- transmission rate
* ψ -- removal rate
* p -- sampling probability upon removal

Epidemiological parameters:
* R<sub>0</sub>=λ/ψ -- reproduction number
* 1/ψ -- infectious time

## BDEI
2 compartments: 
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

## Installation
To install treesimulator:
```bash
pip3 install treesimulator
```

## Usage
### Command line 

## BD
The following command simulates a tree with 200-500 tips under BD model, with λ=0.5, ψ=0.25, p=0.5, 
and saves it to a file tree.nwk, while saving the parameters to a comma-separated file params.csv:
```bash
generate_bd --min_tips 200 --max_tips 500 --la 0.5 --psi 0.25 --p 0.5 --nwk tree.nwk --log params.csv
```
To seee detailed options, run:
```bash
generate_bd --help
```

## BDEI
The following command simulates a tree with 200-500 tips under BDEI model, with μ=1, λ=0.5, ψ=0.25, p=0.5, 
and saves it to a file tree.nwk, while saving the parameters to a comma-separated file params.csv:
```bash
generate_bdei --min_tips 200 --max_tips 500 --mu 1 --la 0.5 --psi 0.25 --p 0.5 --nwk tree.nwk --log params.csv
```
To seee detailed options, run:
```bash
generate_bdei --help
```


## BDSS
The following command simulates a tree with 200-500 tips under BDSS model, 
with λ<sub>nn</sub>=0.1, λ<sub>ns</sub>=0.3, λ<sub>sn</sub>=0.5, λ<sub>ss</sub>=1.5, ψ=0.25, p=0.5, 
and saves it to a file tree.nwk, while saving the parameters to a comma-separated file params.csv:
```bash
generate_bdss --min_tips 200 --max_tips 500 --la_nn 0.1 --la_ns 0.3 --la_sn 0.5 --la_ss 1.5 --psi 0.25 --p 0.5 --nwk tree.nwk --log params.csv
```
To seee detailed options, run:
```bash
generate_bdss --help
```

## User-defined MTBD model
The following command simulates a tree with 200-500 tips under a generic MTBD model, with two states A and B, 
with μ<sub>aa</sub>=0.5, μ<sub>ab</sub>=0.6, μ<sub>ba</sub>=0.7, μ<sub>bb</sub>=0.8, 
λ<sub>aa</sub>=0.1, λ<sub>ab</sub>=0.2, λ<sub>ba</sub>=0.3, λ<sub>bb</sub>=0.4, 
ψ<sub>a</sub>=0.05, ψ<sub>b</sub>=0.08,
p=<sub>a</sub>0.15, p=<sub>b</sub>0.65,
and saves it to a file tree.nwk, while saving the parameters to a comma-separated file params.csv:
```bash
generate_mtbd --min_tips 200 --max_tips 500 --nwk tree.nwk --log params.csv \
--states A B \
--transition_rates 0.5 0.6 0.7 0.8 \
--transmission_rates 0.1 0.2 0.3 0.4 \
--removal_rates 0.05 0.08 \
--sampling_probabilities 0.15 0.65
```
To seee detailed options, run:
```bash
generate_mtbd --help
```


### Python3
To simulate trees with 200-500 tips under the above models and settings:
```python
from treesimulator.generator import generate
from treesimulator import save_forest
from treesimulator.mtbd_models import Model, BirthDeathModel, BirthDeathExposedInfectiousModel, BirthDeathWithSuperSpreadingModel

bd_model = BirthDeathModel(p=0.5, la=0.5, psi=0.25)
print(bd_model.get_epidemiological_parameters())
[bd_tree], (bd_total_tips, _, bd_total_time) = generate(bd_model, min_tips=200, max_tips=500)
save_forest([bd_tree], 'BD_tree.nwk')

bdei_model = BirthDeathExposedInfectiousModel(p=0.5, mu=1, la=0.5, psi=0.25)
print(bdei_model.get_epidemiological_parameters())
[bdei_tree], (bdei_total_tips, _, bdei_total_time) = generate(bdei_model, min_tips=200, max_tips=500)
save_forest([bdei_tree], 'BDEI_tree.nwk')

bdss_model = BirthDeathWithSuperSpreadingModel(p=0.5, la_nn=0.1, la_ns=0.3, la_sn=0.5, la_ss=1.5, psi=0.25)
print(bdss_model.get_epidemiological_parameters())
[bdss_tree], (bdss_total_tips, _, bdss_total_time) = generate(bdss_model, min_tips=200, max_tips=500)
save_forest([bdss_tree], 'BDSS_tree.nwk')

mtbd_model = Model(states=['A', 'B'], transition_rates=[[0.5, 0.6], [0.7, 0.8]], 
                   transmission_rates=[[0.1, 0.2], [0.3, 0.4]],
                   removal_rates=[0.05, 0.08], ps=[0.15, 0.65])
[mtbd_tree], (mtbd_total_tips, _, mtbd_total_time) = generate(mtbd_model, min_tips=200, max_tips=500)
save_forest([mtbd_tree], 'MTBD_tree.nwk')
```
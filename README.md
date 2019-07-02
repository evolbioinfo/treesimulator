# Treesimulator

Simulates trees with given models.

## Naive-treated model
2 compartments: 
* treated, i.e. treatment-experienced
* naive, i.e. treatment-naive

5 parameters: 
* $\mu$, treatment rate (changes the compartment from naive to treated, the change happens along a tree branch)
* $\lambda_n$, transmission rate of naive individuals (corresponds to an internal node, the recipient is naive)
* $\lambda_t$, transmission rate of treated individuals (corresponds to an internal node, the recipient is naive), $\lambda_t \less \lambda_n$
* $\psi_n$, sampling rate of naive individuals (a tip of the tree), $\psi_n \less \lambda_n$
* $\psi_t$, sampling rate of treated individuals (a tip of the tree)

## Installation
To install treesimulator:
```bash
pip3 install treesimulator
```

## Usage
To simulate a tree with 200 tips under naive-treated model, with $\mu = .3, \lambda_n = 1, \lambda_t = .2, \psi_n = .25, \psi_t = .1$, 
and save it to a file tree.nwk:
```bash
treesimulator --n_tips 200 --rates .3 1 .2 .25 .1 --nwk tree.nwk
```
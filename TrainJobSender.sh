#!/bin/bash

sbatch train_adab.sh
sbatch train_bayesglm.sh
sbatch train_nn.sh
sbatch train_adatree.sh
sbatch train_cforest.sh
sbatch train_mlnn.sh
sbatch train_obl.sh
#!/bin/bash
#bsub -q cljymgi -n 1 -o out_0 -e err_0 -J rxnp_0 bash run.sh 0
bsub -q cljymgi -n 24 -o out_eva -e err_eva -J rxnp_eva bash run_evaluate.sh
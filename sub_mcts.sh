#!/bin/bash
bsub -q cljymgi -n 2 -o out_mcts_0 -e err_mcts_0 -J rxnp_0 bash run_mcts.sh 0
bsub -q cljymgi -n 2 -o out_mcts_1 -e err_mcts_1 -J rxnp_1 bash run_mcts.sh 1
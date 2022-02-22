#!/bin/bash
source activate /home/hfnl/zbc/.conda/envs/rxnp
echo "source activate /home/hfnl/zbc/.conda/envs/rxnp" >> log_mcts
/home/hfnl/zbc/.conda/envs/rxnp/bin/python mcts_process.py $1
echo "python mcts_process.py $1" >> log_mcts
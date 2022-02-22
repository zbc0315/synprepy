#!/bin/bash
source activate /home/hfnl/zbc/.conda/envs/rxnp
echo "source activate /home/hfnl/zbc/.conda/envs/rxnp" >> log_evaluate
/home/hfnl/zbc/.conda/envs/rxnp/bin/python mcts_evaluate_process.py
echo "python mcts_evaluate_process.py" >> log_evaluate
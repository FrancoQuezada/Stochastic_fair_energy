# Multistage Fair Energy Allocation

This repository contains a Julia implementation of a multistage stochastic energy allocation model with fairness criteria for residential energy communities. The code builds scenario-tree instances from photovoltaic production data, generates stochastic demand and price inputs, and solves linear optimization models under several fairness notions.

## Main features

- Multistage stochastic energy allocation model on a scenario tree
- Linear programming formulations solved with CPLEX and JuMP
- Fairness modes:
  - `NONE`
  - `PEA` / proportional photovoltaic energy allocation
  - `SA` / proportional savings allocation
  - `LEXMMFPEA` / lexicographic max-min fairness for PEA
  - `LEXMMFSA` / lexicographic max-min fairness for SA
- CSV reporting at system level and per-house level
- Shell scripts for single runs and parallel experiment batches
- Deterministic instance generation for the same input configuration

## Repository structure

- [multi.jl](/home/franco/Documentos/Multistage%20Energy/Code_stochastique/multi.jl): main multistage model, reporting utilities, and experiment entry points
- [parametersMS.jl](/home/franco/Documentos/Multistage%20Energy/Code_stochastique/parametersMS.jl): multistage instance generation, stochastic PV, prices, and demands
- [structuresMulti.jl](/home/franco/Documentos/Multistage%20Energy/Code_stochastique/structuresMulti.jl): scenario-tree and solution data structures
- [mmf_pea.jl](/home/franco/Documentos/Multistage%20Energy/Code_stochastique/mmf_pea.jl): lexicographic max-min routine for photovoltaic allocation
- [mmf_sa.jl](/home/franco/Documentos/Multistage%20Energy/Code_stochastique/mmf_sa.jl): lexicographic max-min routine for savings allocation
- [script.sh](/home/franco/Documentos/Multistage%20Energy/Code_stochastique/script.sh): generic execution and reporting driver
- [run_parallel_fairness.sh](/home/franco/Documentos/Multistage%20Energy/Code_stochastique/run_parallel_fairness.sh): launches multiple fairness and tree configurations in parallel
- `script_<FAIRNESS>_S<S>_C<C>_P<P>.sh`: dedicated scripts for a specific fairness mode and scenario-tree structure
- `inst/inst2020/`: input photovoltaic data instances

## Requirements

- Julia
- CPLEX
- Julia packages used by the project:
  - `JuMP`
  - `CPLEX`
  - `CSV`
  - `DataFrames`
  - `Distributions`
  - `Statistics`
  - `Random`

## Run one configuration

Example: run one instance with proportional savings fairness on tree `S=6, C=2, P=4`.

```bash
INSTANCE_FROM=7 INSTANCE_TO=7 \
J_SET=5 THETA_SET=0.6 AVG_D_SET=100.0 DEV_D_SET=10.0 \
./script_SA_S6_C2_P4.sh
```

Example: run one instance with lexicographic max-min fairness for PEA.

```bash
INSTANCE_FROM=1 INSTANCE_TO=1 \
J_SET=5 THETA_SET=0.2 AVG_D_SET=100.0 DEV_D_SET=10.0 \
./script_LEXMMFPEA_S8_C2_P3.sh
```

## Run a full experiment batch

Run the default parallel batch:

```bash
./run_parallel_fairness.sh
```

Print the commands without executing them:

```bash
DRY_RUN=1 ./run_parallel_fairness.sh
```

Run a custom batch:

```bash
INST_FOLDER=inst/inst2020 \
INSTANCE_FROM=1 \
INSTANCE_TO=100 \
FAIRNESS_LIST='NONE,PEA,SA,LEXMMFPEA,LEXMMFSA' \
TREE_LIST='6:2:4,6:6:4,8:2:3,8:6:3' \
J_SET='5,10' \
THETA_SET='0.2,0.6' \
AVG_D_SET='100.0' \
DEV_D_SET='10.0,20.0' \
./run_parallel_fairness.sh
```

## Outputs

Each fairness mode writes two CSV reports:

- summary report: one row per solved configuration
- house report: one row per house and solved configuration

Current summary columns include:

- `InstanceNo`
- `InstanceFile`
- `NBstage`, `Childs`, `Periods`
- `Theta`, `Avg_d`, `Dev_d`, `J`
- `Fairness`
- `InstFolder`
- `Status`
- `ExpectedCost`
- `ScenarioCostExpected`, `ScenarioCostMin`, `ScenarioCostMax`
- `RegretExpected`, `RegretMin`, `RegretMax`
- `RunTimeSec`, `ModelSolveTimeSec`

## Reproducibility

The multistage instance generator uses a deterministic seed derived from:

- input file name
- number of stages
- tree branching factor
- periods per stage
- number of houses
- stochastic parameters `theta`, `avg`, and `dev`

This means that the same configuration generates the same instance across runs and across fairness modes.

## Notes

- The base multistage model is solved as a linear program.
- `LEXMMFPEA` and `LEXMMFSA` are implemented as sequences of LPs, not as a single model.
- The main model currently uses a one-hour CPLEX time limit.
- The lexicographic SA routine includes additional console logging to help diagnose infeasibility between lexicographic steps.

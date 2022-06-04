# Restricted-TRL
Restricted Tree-based Reinforcement Learning Simulation

- `TRL Functions.R` stores functions for original TRL method.
- `Assignx3.R` stores true treatment assignment and true optimal regimes for the simulation.

/simulation_code
  - `simulation_true_pi_RTRL_NaivePatial.R` is the main simulation code comparing the RT-RL method with naive T-RL with partial non-viable records and uses true pi models.
  - `simulation_wrong_pi_RTRL_NaivePatial.R` is the main simulation code comparing the RT-RL method with naive T-RL with partial non-viable records and uses wrong pi models.

/additional_simulations
  - `simulation_true_pi_RTRL_NaiveFull.R` is the main simulation code comparing the RT-RL method with naive T-RL deleting all non-viable records and uses true pi models.
  - `simulation_wrong_pi_RTRL_NaiveFull.R` is the main simulation code comparing the RT-RL method with naive T-RL deleting all non-viable records and uses wrong pi models.

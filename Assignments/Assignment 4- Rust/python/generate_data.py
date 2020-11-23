import numpy as np
import pandas as pd
from rust import RustModel

# For simulations
n_buses = 50
n_periods= 150

# Setup the problem with a grid
rust=RustModel(x_min=0,x_max=200,n_grid=201,beta=0.975)
# Initial parameters
rust.set_beta(0.975)
rust.set_parameters([2.4569, 0.03  ,11.7257])
rust.update_transition_probs([ 0.0937, 0.4475, 0.4459, 0.0127, 0.0002 ])
rust.solve_value_function()
rust.compute_choice_probabilities()

df=pd.concat([pd.DataFrame(rust.simulate_bus(n_periods),columns=['period_id','y_it','x_it']) for i in range(0,n_buses)],axis=0)
df['bus_id']=(df.period_id < df.period_id.shift()).cumsum()
df.to_csv('../rust_data_2020.csv')

# Check that the transitions make sense
dx=(df.x_it-df.x_it.shift(1))
dx[dx>=0].value_counts(normalize=True).sort_index()
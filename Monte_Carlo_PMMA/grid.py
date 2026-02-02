import numpy as np


#%%
EE = np.logspace(0, 4.4, 1000)
EE_prec = np.logspace(-1, 4.4, 1000)

THETA_deg = np.linspace(0, 180, 1000)
THETA_rad = np.deg2rad(THETA_deg)

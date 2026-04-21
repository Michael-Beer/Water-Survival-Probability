import os
from pathlib import Path
import MDAnalysis
from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP
import matplotlib.pyplot as plt
import numpy as np

# Remove existing output
output_file = Path("SP_combined.txt")
if output_file.is_file():
    os.remove(output_file)

# Input: topology files and trajectories for ConfA and ConfB
topology_A = "ConfA/Run1/system.parm7"
trajectories_A = [
    "ConfA/Run1/prod.nc",
    "ConfA/Run2/prod.nc",
    "ConfA/Run3/prod.nc",
]

topology_B = "ConfB/Run1/system.parm7"
trajectories_B = [
    "ConfB/Run1/prod.nc",
    "ConfB/Run2/prod.nc",
    "ConfB/Run3/prod.nc",
]

# Merge into one list (with matching topology)
datasets = [(topology_A, t) for t in trajectories_A] + \
           [(topology_B, t) for t in trajectories_B]

# Selection string
selection = "resname WAT and sphzone 3.5 ((resid 43 and name C29) or (resid 139 and name CD))"

# Parameters
start = 0
stop = 15001
tau_max = 1000

# Storage
all_sp = []
tau_timeseries = None

# Loop over datasets
for topo, traj in datasets:
    print(f"Processing trajectory: {traj}")
    universe = MDAnalysis.Universe(topo, traj)
    sp = SP(universe, selection, verbose=True)
    sp.run(start=start, stop=stop, tau_max=tau_max)

    if tau_timeseries is None:
        tau_timeseries = sp.tau_timeseries
    all_sp.append(sp.sp_timeseries)

# Convert list to array
all_sp = np.array(all_sp)
# Mean survival probability
sp_avg = np.mean(all_sp, axis=0)

# Standard error of mean (SEM) per timepoint
N = all_sp.shape[0]
sp_sem = np.std(all_sp, axis=0, ddof=1) / np.sqrt(N)

# Save results (tau, mean, sem at each point)
with open("SP_combined.txt", "w") as file:
    for tau, mean, sem in zip(tau_timeseries, sp_avg, sp_sem):
        file.write(f"{tau} {mean} {sem}\n")

# -------- Plotting with shaded SEM band --------
plt.figure(figsize=(6,4))

# Plot mean line
plt.plot(tau_timeseries, sp_avg, '-', color='blue', linewidth=1.5, label='Average')

# Shaded area for SEM
plt.fill_between(
    tau_timeseries,
    sp_avg - sp_sem,
    sp_avg + sp_sem,
    color='blue',
    alpha=0.2,
    label='± SEM'
)

plt.xlabel('Time (ps)')
plt.ylabel('Survival Probability')
#plt.title('Survival Probability (ConfA + ConfB)')
plt.legend()
plt.tight_layout()
plt.savefig("SP.png", dpi=300)
plt.show()


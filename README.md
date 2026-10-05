# Water Survival Probability Analysis

A Python script that calculates the **survival probability (SP) of water molecules** in a defined region of an enzyme active site across multiple molecular dynamics (MD) replicates, using [MDAnalysis](https://www.mdanalysis.org/). Results from several independent runs (and multiple conformations) are pooled into a mean survival probability curve with the standard error of the mean (SEM) at each time lag.

## What it does

1. Loads AMBER topology (`.parm7`) and trajectory (`.nc`) files for each replicate of two systems (`ConfA` and `ConfB`).
2. Selects water molecules (`WAT`) within a sphere (default 3.5 Å) around a user-defined set of atoms.
3. Computes the survival probability of those waters for each trajectory using `MDAnalysis.analysis.waterdynamics.SurvivalProbability`.
4. Averages across all trajectories and computes the SEM per time lag.
5. Writes the combined results to a text file and produces a publication-ready plot with a shaded SEM band.

## Requirements

- Python 3.8+
- [MDAnalysis](https://docs.mdanalysis.org/stable/documentation_pages/analysis/waterdynamics.html)
- NumPy
- Matplotlib

```bash
pip install MDAnalysis numpy matplotlib
```

## Directory layout

The script expects the following structure, run from the parent directory:

```
.
├── survival_probability.py
├── ConfA/
│   ├── Run1/
│   │   ├── system.parm7
│   │   └── prod.nc
│   ├── Run2/
│   │   └── prod.nc
│   └── Run3/
│       └── prod.nc
└── ConfB/
    ├── Run1/
    │   ├── system.parm7
    │   └── prod.nc
    ├── Run2/
    │   └── prod.nc
    └── Run3/
        └── prod.nc
```

The topology from `Run1` of each configuration is used for all replicates of that configuration.

## Usage

```bash
python survival_probability.py
```

Progress is printed for each trajectory. Any existing `SP_combined.txt` is deleted at the start of each run.

## Configuration

All settings are defined at the top of the script.

| Variable | Default | Description |
|---|---|---|
| `topology_A`, `topology_B` | `ConfX/Run1/system.parm7` | Topology file for each configuration |
| `trajectories_A`, `trajectories_B` | `ConfX/Run{1-3}/prod.nc` | Trajectory files for each configuration |
| `selection` | see below | MDAnalysis selection string for the water region |
| `start` | `0` | First frame analysed |
| `stop` | `15001` | Last frame (exclusive) analysed |
| `tau_max` | `1000` | Maximum time lag (in frames) |

### Selection string

```python
selection = "resname WAT and sphzone 3.5 ((resid 43 and name C29) or (resid 139 and name CD))"
```

This selects water molecules within 3.5 Å of the union of two reference atoms (residue 43 atom `C29` and residue 139 atom `CD`). Edit the residue numbers, atom names and radius to suit your system. Note that the selection is re-evaluated at every frame, so the water population is dynamic.

## Output

### `SP_combined.txt`

Space-separated text file with one row per time lag:

```
tau  mean_SP  SEM
```

### `SP.png`

A 300 dpi plot of the mean survival probability versus time lag, with the SEM shown as a shaded band.

## Statistics

- The **mean** is taken across all trajectories (by default 3 × ConfA + 3 × ConfB = 6).
- The **SEM** is calculated per time lag as `std(ddof=1) / sqrt(N)`, where `N` is the number of trajectories.
- Trajectories are treated as independent replicates; ConfA and ConfB are pooled rather than analysed separately. To compare the two configurations, run the script once per configuration, or split the `datasets` list.

## Notes

- `tau_timeseries` returned by MDAnalysis is in **frames**. To plot in picoseconds, multiply by the time between saved frames (e.g. `tau * dt`, where `dt` is `universe.trajectory.dt`). Check that the x-axis label matches your units.
- All trajectories must have the same number of frames in the analysed window so that the SP arrays have equal length.
- `tau_max` should not exceed the number of frames analysed.

## References

- Michaud-Agrawal, N. et al. *MDAnalysis: A toolkit for the analysis of molecular dynamics simulations.* J. Comput. Chem. 32, 2319–2327 (2011).
- Gowers, R. J. et al. *MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations.* Proc. 15th Python in Science Conf. (2016).
- Liu, P., Harder, E. & Berne, B. J. *On the calculation of diffusion coefficients in confined fluids and interfaces with an application to the liquid–vapor interface of water.* J. Phys. Chem. B 108, 6595–6602 (2004).

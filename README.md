# Motion-Robust Gradient-Echo T2 Mapping in the Liver Using Readout Alternation (ROA)

This repository provides simulation, analysis, and visualization code accompanying the paper:

> **"Motion-robust gradient-echo T2 mapping in the liver using readout alternation"**  
> Daiki Tamada, Amirhossein Roshanshad, Lukas Müller, Jitka Starekova, Diego Hernando, Scott B. Reeder  
> *Department of Radiology, University of Wisconsin-Madison, Madison, WI, USA*  
> Correspondence: Daiki Tamada, PhD (dtamada@wisc.edu)

---

## Overview

Quantitative liver $T_2$ mapping provides critical biomarkers for hepatic iron overload, fibrosis, and inflammation. Rapid 3D gradient-echo (GRE)-based techniques, notably **double-echo steady-state (DESS)** and **phase-based $T_2$ mapping (PBT2)**, enable volumetric whole-liver coverage within a single breath-hold. However, unbalanced gradient moments make these steady-state sequences exquisitely sensitive to residual physiological motion (e.g., cardiac pulsation and respiratory drift), resulting in severe phase dispersion, signal loss, and biased $T_2$ estimation.

**Readout Alternation (ROA)** is an acquisition strategy that alternates the readout gradient polarity between consecutive TRs in combination with an orthogonal spoiler gradient. This mitigates motion-induced dephasing and stabilizes signal formation without acquisition time penalty.

This repository contains:
- **Monte Carlo Bloch simulations** evaluating $T_2$ estimation accuracy across motion velocities (Figure 3).
- **Transient and steady-state Bloch simulations** examining orthogonal spoiler moment requirements (Supporting Information Figure S1).
- **Extended Phase Graph (EPG) pathway-by-pathway analysis** quantifying dominant echo pathways and gradient first moments ($M_1$) (Supporting Information Table S1).
- **Shared simulation tools and dictionary-generation utilities** for DESS and PBT2.

---

## Repository Structure

```
ROA_paper_code/
├── Figure3/
│   ├── MC_simulation.jl          # Monte-Carlo Bloch simulations of motion in Julia
│   ├── RF_spoiling_2D_fast.jl     # High-performance 2D Bloch isochromat simulation in Julia
│   └── Figure3.m                 # MATLAB script to generate Figure 3 plots
├── Figure_S1/
│   └── FigureS1.m                # MATLAB script to simulate and plot Figure S1
├── TableS1/
│   ├── Table_S1.m                # MATLAB script for EPG pathway tracing
│   ├── epg_steady_state_last8.m  # EPG simulation engine (traces last 8 TRs of steady state)
│   ├── summarize_steady_v4.py    # Python script summarizing dominant pathways & moments
│   ├── steady_initialized_pathways_conventional.csv  # Precomputed EPG pathways (Conventional)
│   └── steady_initialized_pathways_roa.csv           # Precomputed EPG pathways (ROA)
├── tools/                        # Shared MATLAB utility library
│   ├── RF_spoiling_2D.m          # 2D Bloch simulation with RF spoiling
│   ├── RF_spoiling_2D_motion.m   # 2D Bloch simulation with motion profiles and gradient moments
│   ├── build_LUTs_DESS_Bloch.m   # Bloch-based dictionary (LUT) generation for DESS & PBT2
│   ├── T2map_phase.m             # Quantitative T2 mapping via LUT inversion (magnitude/phase)
│   ├── build_motion.m            # Motion profile generator (linear, sine, random)
│   ├── freeprecess.m             # Free precession matrix operator
│   ├── throt.m                   # Arbitrary-axis rotation matrix operator
│   ├── xrot.m                    # X-axis rotation matrix operator
│   └── zrot.m                    # Z-axis rotation matrix operator
└── README.md
```

---

## Prerequisites and Environment Setup

### 1. MATLAB
- Tested on MATLAB R2020b or later.
- Core functions require base MATLAB and the Parallel Computing Toolbox (for `parfor` speedup in dictionary generation, though it will run sequentially if unavailable).

### 2. Julia
- Julia 1.6 or higher.
- Required Julia packages:
  ```julia
  using Pkg
  Pkg.add(["Statistics", "MAT"])
  ```

### 3. Python
- Python 3.8 or higher.
- Required packages:
  ```bash
  pip install numpy pandas
  ```

---

## Detailed Directory and Code Descriptions

### 1. `Figure3/`: Monte Carlo Simulations of Motion Effects

This folder contains code to reproduce **Figure 3** of the manuscript: *"Bloch equation simulations of motion effects on DESS and PBT2 T2 estimation"*.

#### Methodological Summary
- Nominal baseline tissue parameters: $T_1 = 809\text{ ms}$, $T_2 = 34\text{ ms}$ (representative of liver tissue at 3.0T).
- Motion is modeled along the readout direction as a sinusoidal velocity profile with peak velocities $v_0$ ranging from $0$ to $1.0\text{ cm/s}$ in steps of $0.025\text{ cm/s}$.
- At each velocity step, $N = 100$ random trials are evaluated with random respiratory frequencies ($f \in [0.2, 1.2]\text{ Hz}$) and random initial phases ($t_0 \in [0, 2\pi]$).
- Compares conventional acquisition vs. ROA for:
  - **DESS**: Signal ratio $S_R = |F_1| / |F_0|$ mapped to $T_2$ via magnitude dictionary lookup.
  - **PBT2**: Phase of $S_0$ mapped to $T_2$ via phase dictionary lookup.

#### Files
- `Figure3/MC_simulation.jl`: Performs the Monte Carlo Bloch equation simulations in Julia. It evaluates both DESS and PBT2 with and without ROA across all motion realizations and exports the resulting signal arrays (`PH`, `PH_ROA`, `SR`, `SR_ROA`) to a MATLAB `.mat` file (`Simulation_Results.mat`).
- `Figure3/RF_spoiling_2D_fast.jl`: Core Julia function implementing vectorized 2D isochromat Bloch simulation under RF spoiling, gradient moments, and sinusoidal motion.
- `Figure3/Figure3.m`: MATLAB visualization script. It:
  1. Loads `simulation_results.mat` (or `Simulation_Results.mat`).
  2. Generates the Bloch lookup tables (LUTs) using `tools/build_LUTs_DESS_Bloch.m`.
  3. Reconstructs $T_2$ estimates using `tools/T2map_phase.m`.
  4. Plots mean estimated $T_2$ vs. peak velocity with shaded bands representing $\pm 1$ standard deviation for both DESS and PBT2.

#### How to Run
1. Run the Julia Monte Carlo simulation:
   ```bash
   cd Figure3
   julia MC_simulation.jl
   ```
   *(Note: Ensure the exported file is named `simulation_results.mat` for MATLAB, or rename `Simulation_Results.mat` to `simulation_results.mat` if on case-sensitive filesystems).*
2. In MATLAB, run:
   ```matlab
   cd('Figure3')
   Figure3
   ```

---

### 2. `Figure_S1/`: Orthogonal Spoiler Moment and Signal Evolution

This folder contains code to reproduce **Figure S1** of the Supporting Information: *"Effect of orthogonal spoiler moment on signal evolution"*.

#### Methodological Summary
- Evaluates the transient and steady-state behavior of DESS (FISP $F_0$ and PSIF $F_1$) and PBT2 as a function of the orthogonal spoiler gradient moment.
- Evaluates four spoiler moments per TR: $0$, $0.25 \times 2\pi$, $0.5 \times 2\pi$, and $1.0 \times 2\pi$.
- Demonstrates that for DESS, smaller moments alter the steady-state amplitudes, whereas $\ge 0.5 \times 2\pi$ produces stable signal evolution. For PBT2, spoiler moments $< 0.5 \times 2\pi$ produce persistent magnitude and phase oscillations, while $1.0 \times 2\pi$ guarantees stable, oscillation-free steady-state formation, validating the $1.0 \times 2\pi$ spoiler used in the ROA sequence design.

#### Files
- `Figure_S1/FigureS1.m`: Standalone MATLAB script that simulates 600 RF excitations across 961 isochromats using `tools/RF_spoiling_2D_motion.m` and plots signal intensity vs. RF pulse index for each spoiler moment condition.

#### How to Run
In MATLAB:
```matlab
cd('Figure_S1')
FigureS1
```

---

### 3. `TableS1/`: EPG Pathway-by-Pathway Analysis

This folder contains code to reproduce **Table S1** of the Supporting Information: *"Dominant steady-state signal pathways and corresponding first gradient moments for conventional and ROA acquisitions"*.

#### Methodological Summary
- Uses the Extended Phase Graph (EPG) formalism to decompose steady-state magnetization after 1500 TRs into individual coherence pathways, tracing the detailed evolution across the final 8 TRs.
- Computes:
  - Echo type (SE-like vs. STE-like) and effective echo time.
  - Relative pathway power (%) and cumulative pathway power (accounting for $\ge 97.7\%$ and $\ge 98.0\%$ of total power for conventional and ROA, respectively).
  - Effective first gradient moments ($M_1$) along the readout ($M_{1,\text{RO}}$) and slice ($M_{1,\text{SL}}$) directions in units of $10^{-3}\text{ T}\cdot\text{s}^2/\text{m}$.
  - Root-sum-of-squares total first moment ($M_{1,\text{RSS}}$).
- Highlights the mechanism of ROA: the dominant ROA echo component exhibits near-complete cancellation of the readout first moment ($M_{1,\text{RO}} \approx -0.595 \times 10^{-3}\text{ T}\cdot\text{s}^2/\text{m}$ compared to $892.742 \times 10^{-3}\text{ T}\cdot\text{s}^2/\text{m}$ in conventional DESS).

#### Files
- `TableS1/Table_S1.m`: MATLAB script that runs `epg_steady_state_last8.m` for conventional (no slice spoiler) and ROA (with slice spoiler) configurations, saving detailed pathway lists to:
  - `steady_initialized_pathways_conventional.csv`
  - `steady_initialized_pathways_roa.csv`
- `TableS1/epg_steady_state_last8.m`: Core EPG simulation engine implementing pathway tracing, state transitions, relaxation, and gradient moment tracking.
- `TableS1/summarize_steady_v4.py`: Python script that reads the pathway CSV files, sorts and filters pathways by power, computes scaled moments, and exports formatted summary tables matching Table S1:
  - `table_s1_dominant_pathways_conventional.csv`
  - `table_s1_dominant_pathways_roa.csv`
- Precomputed CSV files:
  - `steady_initialized_pathways_conventional.csv`
  - `steady_initialized_pathways_roa.csv`

#### How to Run
1. Generate the pathway CSV files using MATLAB (or use the provided precomputed files):
   ```matlab
   cd('TableS1')
   Table_S1
   ```
2. Summarize and format the dominant pathway tables using Python:
   ```bash
   cd TableS1
   python3 summarize_steady_v4.py --n-dominant 12 --cumulative-power 0.90
   ```
   Command-line options:
   - `--conventional`: Path to conventional pathway CSV (default: `./steady_initialized_pathways_conventional.csv`).
   - `--roa`: Path to ROA pathway CSV (default: `./steady_initialized_pathways_roa.csv`).
   - `--output-dir`: Output directory for generated summary tables (default: `./`).
   - `--n-dominant`: Number of top dominant pathways to include (default: `12`).
   - `--cumulative-power`: Cumulative power threshold cutoff (default: `0.90`).
   - `--allow-conventional-sl`: Retain slice moment for conventional sequence (by default forced to 0).

---

### 4. `tools/`: Shared MATLAB Utilities

The `tools/` directory provides core simulation functions and dictionary models shared across scripts:

| File | Description |
|---|---|
| [`RF_spoiling_2D.m`](tools/RF_spoiling_2D.m) | 2D isochromat Bloch simulation with quadratic/linear RF phase increments, yielding transient and steady-state signals ($F_0$ and $F_1$). Supports alternating spoiler polarity. |
| [`RF_spoiling_2D_motion.m`](tools/RF_spoiling_2D_motion.m) | 2D isochromat Bloch simulation incorporating 1D/2D motion profiles, gradient moments, and phase modulation. |
| [`build_LUTs_DESS_Bloch.m`](tools/build_LUTs_DESS_Bloch.m) | Generates Bloch simulation lookup tables (LUTs) across a range of $T_2$ values for DESS magnitude ratios ($|F_1|/|F_0|$) and PBT2 phases ($\angle S_0$). |
| [`T2map_phase.m`](tools/T2map_phase.m) | Inverts measured DESS signal ratios (magnitude mode) or PBT2 phase angles (phase mode) into quantitative $T_2$ maps using the precomputed LUTs. |
| [`build_motion.m`](tools/build_motion.m) | Generates displacement and velocity trajectories for linear, sinusoidal, and random motion models as a function of time. |
| [`freeprecess.m`](tools/freeprecess.m) | Computes the state transition matrix $A$ and recovery vector $B$ for free precession under relaxation ($T_1, T_2$) and off-resonance ($\Delta f$). |
| [`throt.m`](tools/throt.m) | Calculates the 3D rotation matrix for a flip angle $\phi$ about an arbitrary phase angle $\theta$. |
| [`xrot.m`](tools/xrot.m) | Calculates the 3D rotation matrix for an RF pulse applied along the $x$-axis. |
| [`zrot.m`](tools/zrot.m) | Calculates the 3D rotation matrix for precession about the $z$-axis. |

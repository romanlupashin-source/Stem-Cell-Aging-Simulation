# Stem Cell Aging Simulation 

This project implements a stochastic simulation of stem cell aging based on
telomere dynamics across multiple cell divisions.

The model compares normal somatic cells with telomerase-active
(cancer-like) cells using repeated Monte Carlo simulations.

---

## Scientific Background

Telomeres shorten during cell division due to the end-replication problem.
When telomeres reach a critical length (Hayflick limit), cells enter
senescence or apoptosis.

Some cells (e.g. cancer cells) activate telomerase, partially or fully
compensating telomere loss and enabling replicative immortality.

---

## Model Description

- Random telomere loss per division
- Multiple independent simulation runs
- Mean trajectory and interquartile range (IQR)
- Comparison between:
  - Normal stem cells
  - Telomerase-active cells

---

## Output

The simulation produces:
- Mean telomere length trajectories
- Variability bands (IQR)
- Example single-cell trajectories
- Estimated divisions to senescence

---

## Technologies

- Python 3.9+
- NumPy
- Matplotlib

---

## Usage

```bash
pip install -r requirements.txt
python stem_cell_aging_simulation.py

⚠️ Limitations

This is a simplified stochastic model and does not capture
DNA damage response, cell cycle checkpoints, or tissue-level effects.

Author

Roman Lupashin

import numpy as np
import matplotlib.pyplot as plt


# ---------------- PARAMETERS ----------------
SEED = 42
START_LENGTH = 10_000          # Initial telomere length (bp)
CRITICAL_LIMIT = 500           # Hayflick limit
MAX_DIVISIONS = 100
LOSS_MIN, LOSS_MAX = 50, 150
N_RUNS = 60

# Telomerase parameters (cancer-like model)
RESTORE_FACTOR = 1.0
RESTORE_NOISE = 30

rng = np.random.default_rng(SEED)


def simulate_normal(start, critical, max_divs, rng):
    """Simulate one normal somatic cell lifecycle."""
    tel = start
    hist = [tel]

    for i in range(max_divs):
        loss = rng.integers(LOSS_MIN, LOSS_MAX + 1)
        tel -= loss

        if tel <= critical:
            tel = critical
            hist.append(tel)
            hist += [tel] * (max_divs - i - 1)
            break

        hist.append(tel)

    if len(hist) < max_divs + 1:
        hist += [tel] * (max_divs + 1 - len(hist))

    return np.array(hist, dtype=float)


def simulate_cancer(start, max_divs, rng,
                    restore_factor=1.0, restore_noise=30):
    """Simulate a telomerase-active (cancer-like) cell."""
    tel = start
    hist = [tel]

    for _ in range(max_divs):
        loss = rng.integers(LOSS_MIN, LOSS_MAX + 1)
        noise = rng.integers(-restore_noise, restore_noise + 1)
        restore = int(loss * restore_factor) + noise
        restore = max(restore, 0)

        tel = max(tel - loss + restore, 0)
        hist.append(tel)

    return np.array(hist, dtype=float)


def main():
    normal_runs = np.vstack([
        simulate_normal(START_LENGTH, CRITICAL_LIMIT, MAX_DIVISIO

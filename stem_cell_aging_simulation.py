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
        simulate_normal(START_LENGTH, CRITICAL_LIMIT, MAX_DIVISIONS, rng)
        for _ in range(N_RUNS)
    ])

    cancer_runs = np.vstack([
        simulate_cancer(START_LENGTH, MAX_DIVISIONS, rng,
                        RESTORE_FACTOR, RESTORE_NOISE)
        for _ in range(N_RUNS)
    ])

    x = np.arange(0, MAX_DIVISIONS + 1)

    normal_mean = normal_runs.mean(axis=0)
    normal_p25 = np.percentile(normal_runs, 25, axis=0)
    normal_p75 = np.percentile(normal_runs, 75, axis=0)

    cancer_mean = cancer_runs.mean(axis=0)
    cancer_p25 = np.percentile(cancer_runs, 25, axis=0)
    cancer_p75 = np.percentile(cancer_runs, 75, axis=0)

    plt.figure(figsize=(10, 6))

    plt.plot(x, normal_mean, label="Normal cell (mean)", color="blue", lw=2)
    plt.fill_between(x, normal_p25, normal_p75,
                     color="blue", alpha=0.2, label="Normal IQR")

    plt.plot(x, cancer_mean, label="Cancer cell (mean)", color="red", lw=2)
    plt.fill_between(x, cancer_p25, cancer_p75,
                     color="red", alpha=0.2, label="Cancer IQR")

    for i in range(3):
        plt.plot(x, normal_runs[i], color="blue", alpha=0.4, lw=0.8)
        plt.plot(x, cancer_runs[i], color="red", alpha=0.4, lw=0.8, ls=":")

    plt.axhline(y=CRITICAL_LIMIT, color="black", ls="--",
                label="Hayflick limit")

    plt.title("Stem Cell Aging: Senescence vs Telomerase Activity")
    plt.xlabel("Cell divisions")
    plt.ylabel("Telomere length (bp)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    death_divs = [
        np.where(run <= CRITICAL_LIMIT)[0][0]
        if np.any(run <= CRITICAL_LIMIT) else MAX_DIVISIONS
        for run in normal_runs
    ]

    print(
        f"Mean divisions to senescence: "
        f"{np.mean(death_divs):.1f} ± {np.std(death_divs):.1f} "
        f"(n={N_RUNS})"
    )


if __name__ == "__main__":
    main()

import numpy as np
import matplotlib.pyplot as plt

# --- ПАРАМЕТРЫ ---
seed = 42
start_length = 10000    # Начальная длина теломер (bp)
critical_limit = 500    # Предел Хейфлика (смерть)
max_divisions = 100     # Число делений для каждого прогона
loss_min, loss_max = 50, 150   # Потеря теломер за деление (равномерно)
n_runs = 60             # Сколько независимых прогонов для усреднения

# Параметры теломеразы для "раковой" модели:
# restore_factor = 1.0 -> в среднем восстанавливает столько же, сколько теряется (поддержание)
# restore_noise — случайный шум восстановления (в bp)
restore_factor = 1.0
restore_noise = 30

rng = np.random.default_rng(seed)

def simulate_normal(start, critical, max_divs, rng):
    """Одна симуляция соматической клетки. Возвращает список длины max_divs+1."""
    tel = start
    hist = [tel]
    for i in range(max_divs):
        loss = rng.integers(loss_min, loss_max + 1)
        tel -= loss
        if tel <= critical:
            tel = critical
            hist.append(tel)
            # заполняем остаток, чтобы все истории были одной длины
            hist += [tel] * (max_divs - i - 1)
            break
        hist.append(tel)
    # если не умерла за max_divs, гарантируем длину списка max_divs+1
    if len(hist) < max_divs + 1:
        hist += [tel] * (max_divs + 1 - len(hist))
    return np.array(hist, dtype=float)

def simulate_cancer(start, critical, max_divs, rng, restore_factor=1.0, restore_noise=30):
    """Одна симуляция 'раковой' клетки с действующей теломеразой."""
    tel = start
    hist = [tel]
    for i in range(max_divs):
        loss = rng.integers(loss_min, loss_max + 1)
        # восстановление: пропорция от потери + небольшой шум
        noise = rng.integers(-restore_noise, restore_noise + 1)
        restore = int(loss * restore_factor) + noise
        if restore < 0:
            restore = 0
        tel = tel - loss + restore
        # на практике теломеры не уходят в отрицательное без клетки-механики; просто ограничим
        if tel < 0:
            tel = 0.0
        hist.append(tel)
    return np.array(hist, dtype=float)

# --- СБОР ДАННЫХ (несколько прогонов) ---
normal_runs = np.vstack([simulate_normal(start_length, critical_limit, max_divisions, rng) for _ in range(n_runs)])
# для раковой модели: несколько прогонов с теми же параметрами
cancer_runs = np.vstack([simulate_cancer(start_length, critical_limit, max_divisions, rng, restore_factor, restore_noise) for _ in range(n_runs)])

# --- СТАТИСТИКА ---
x = np.arange(0, max_divisions + 1)  # поколения (включая 0)
normal_mean = normal_runs.mean(axis=0)
normal_p25 = np.percentile(normal_runs, 25, axis=0)
normal_p75 = np.percentile(normal_runs, 75, axis=0)

cancer_mean = cancer_runs.mean(axis=0)
cancer_p25 = np.percentile(cancer_runs, 25, axis=0)
cancer_p75 = np.percentile(cancer_runs, 75, axis=0)

# --- РИСУНОК ---
plt.figure(figsize=(10, 6))

# среднее и интерквартильная область для обычной клетки
plt.plot(x, normal_mean, label='Обычная клетка (mean)', color='blue', linewidth=2)
plt.fill_between(x, normal_p25, normal_p75, color='blue', alpha=0.2, label='Обычные: IQR')

# среднее и интерквартильная область для раковой клетки
plt.plot(x, cancer_mean, label='Раковая клетка (mean)', color='red', linewidth=2)
plt.fill_between(x, cancer_p25, cancer_p75, color='red', alpha=0.2, label='Раковые: IQR')

# Покажем ещё несколько отдельных траекторий (для интуиции)
for i in range(3):
    plt.plot(x, normal_runs[i], color='blue', alpha=0.5, linewidth=0.8)
    plt.plot(x, cancer_runs[i], color='red', alpha=0.5, linewidth=0.8, linestyle=':')

# Линия смерти
plt.axhline(y=critical_limit, color='black', linestyle='--', label='Предел Хейфлика')

plt.title("Сравнение: Старение vs Бессмертие (модель с теломеразой)")
plt.xlabel("Количество делений (поколения)")
plt.ylabel("Длина теломеры (bp)")
plt.ylim(bottom=0)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# --- ПРОСТЫЕ ВЫВОДЫ ---
# Сколько в среднем делений до смерти для обычной клетки в этих параметрах:
# (найдём индекс первого вхождения критического предела в к��ждой симуляции)
death_divisions = [(np.where(run <= critical_limit)[0][0] if np.any(run <= critical_limit) else max_divisions) for run in normal_runs]
print(f"Среднее число делений до сенесценции (обычные): {np.mean(death_divisions):.1f} ± {np.std(death_divisions):.1f} (std), из {n_runs} прогонов")

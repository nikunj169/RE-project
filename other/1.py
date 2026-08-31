# Sweep temperature effect
import numpy as np
import matplotlib.pyplot as plt

# Use mean scaled values
base = X_test_scaled.mean(axis=0)

temps = np.linspace(-3, 3, 200)
preds = []

for t in temps:
    x = base.copy()
    x[1] = t  # temperature index
    sal, temp, aou = x[0], x[1], x[2]

    y = (sal/(3.054312 - aou) + temp - (aou + 6.475268))**2 + 10.15233
    preds.append(y)

plt.plot(temps, preds)
plt.xlabel("Scaled Temperature")
plt.ylabel("Predicted CO₂")
plt.title("Nonlinear Response Curve")
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error


# ================================
# 1️⃣ Load Data
# ================================
df = pd.read_csv("data/processed/southern_ocean_training.csv")

features = ['salinity','temperature','aou','stratification_index','c_atm']
X = df[features].values
y = df['co2_aq'].values

# ================================
# 2️⃣ RECREATE TRAIN SPLIT
# (CRITICAL — scaler trained on this)
# ================================
X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    random_state=42
)

# ================================
# 3️⃣ RECREATE SCALER
# ================================
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Extract scaled variables
sal = X_test_scaled[:,0]
temp = X_test_scaled[:,1]
aou = X_test_scaled[:,2]

# ================================
# 4️⃣ Apply Symbolic Equation
# ================================
pred_full = (
    (sal / (3.054312 - aou) + temp - (aou + 6.475268))**2
    + 10.15233
)

pred_simple = (
    (11.368984 - temp) * (aou - temp)
    + 52.853294
)

# ================================
# 5️⃣ Metrics
# ================================
def metrics(y, yhat):
    return r2_score(y, yhat), np.sqrt(mean_squared_error(y, yhat))

r2_full, rmse_full = metrics(y_test, pred_full)
r2_simple, rmse_simple = metrics(y_test, pred_simple)

print(f"\nFull Model   -> R²: {r2_full:.4f}  RMSE: {rmse_full:.4f}")
print(f"Simple Model -> R²: {r2_simple:.4f}  RMSE: {rmse_simple:.4f}")

# ================================
# 6️⃣ Plot
# ================================
plt.figure(figsize=(8,8))

plt.scatter(y_test, pred_full, alpha=0.4, s=10,
            label=f"Full (R²={r2_full:.3f})")

plt.scatter(y_test, pred_simple, alpha=0.4, s=10,
            label=f"Simple (R²={r2_simple:.3f})")

lims = [
    min(y_test.min(), pred_full.min()),
    max(y_test.max(), pred_full.max())
]

plt.plot(lims, lims, 'r--')

plt.xlabel("Actual CO₂(aq)")
plt.ylabel("Predicted CO₂(aq)")
plt.title("Symbolic Regression Performance (Correctly Scaled)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

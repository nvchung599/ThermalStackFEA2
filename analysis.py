import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv('data.csv')
df['T_ambient'] = 0  # Ambient temperature is 0 in training data

# Compute time step
time = df['time'].values
dt = time[1] - time[0]

# Compute derivatives
dT_A_dt = (df['T_avg_thermistor'].diff() / dt).shift(-1).values[:-1]
dT_B_dt = (df['T_avg_aluminum'].diff() / dt).shift(-1).values[:-1]
dT_C_dt = (df['T_avg_resistor'].diff() / dt).shift(-1).values[:-1]

n = len(df) - 1  # Last index with valid derivatives

# Regression for Mass A
X_A = (df['T_avg_aluminum'] - df['T_avg_thermistor']).values[:n].reshape(-1, 1)
y_A = dT_A_dt
reg_A = LinearRegression(fit_intercept=False).fit(X_A, y_A)
theta_1 = reg_A.coef_[0]

# Regression for Mass B
X_B = np.column_stack((
    df['T_avg_thermistor'] - df['T_avg_aluminum'],
    df['T_avg_resistor'] - df['T_avg_aluminum'],
    -df['T_avg_aluminum']
))[:n]
y_B = dT_B_dt
reg_B = LinearRegression(fit_intercept=False).fit(X_B, y_B)
theta_2, theta_3, theta_4 = reg_B.coef_

# Regression for Mass C
X_C = np.column_stack((df['T_avg_aluminum'] - df['T_avg_resistor'], df['P_resistor']))[:n]
y_C = dT_C_dt
reg_C = LinearRegression(fit_intercept=False).fit(X_C, y_C)
theta_5, theta_6 = reg_C.coef_

# Parameter vector
theta = [theta_1, theta_2, theta_3, theta_4, theta_5, theta_6]
print("Estimated theta:", [f"{t:.4f}" for t in theta])

# Recover original parameters
C_C = 1 / theta_6
R_BC = theta_6 / theta_5
C_B = theta_5 / theta_3 * C_C
R_AB = 1 / (theta_2 * C_B)
R_B_amb = 1 / (theta_4 * C_B)
C_A = 1 / (theta_1 * R_AB)

print("\nOriginal Parameters:")
print(f"C_A = {C_A:.4f}, C_B = {C_B:.4f}, C_C = {C_C:.4f}")
print(f"R_AB = {R_AB:.4f}, R_BC = {R_BC:.4f}, R_B_amb = {R_B_amb:.4f}")

# State-space matrices
A = np.array([
    [-theta_1, theta_1, 0],
    [theta_2, -(theta_2 + theta_3 + theta_4), theta_3],
    [0, theta_5, -theta_5]
])
B = np.array([
    [0, 0],
    [0, theta_4],
    [theta_6, 0]
])
C = np.array([[1, 0, 0]])

# Discretize
A_d = np.eye(3) + A * dt
B_d = B * dt

# Observer gain (tune as needed)
L = np.array([0.5, 0.5, 0.5])

# Observer simulation
hat_x = np.zeros((len(df), 3))
hat_x[0] = [df['T_avg_thermistor'][0], df['T_avg_thermistor'][0], df['T_avg_thermistor'][0]]

for k in range(len(df) - 1):
    u_k = [df['P_resistor'][k], df['T_ambient'][k]]
    y_k = df['T_avg_thermistor'][k]
    hat_x[k + 1] = A_d @ hat_x[k] + B_d @ u_k + L * (y_k - C @ hat_x[k])

# Add estimated T_C
df['T_C_observed'] = hat_x[:, 2]

# Plot
plt.figure(figsize=(12, 6))
plt.plot(df['time'], df['T_avg_resistor'], label='True T_C')
plt.plot(df['time'], df['T_C_observed'], label='Observed T_C')
plt.plot(df['time'], df['T_avg_thermistor'], label='Measured T_A')
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('Observer Performance')
plt.legend()
plt.grid(True)
plt.show()
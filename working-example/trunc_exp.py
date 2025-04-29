import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

rate = 1
maxT = 2
N = 1000000

FMin = np.exp(-rate * maxT)
FMax = 1 - np.exp(-rate * maxT)
U_lower_bound = np.random.uniform(FMin, 1.0, N)
U_upper_bound = np.random.uniform(0.0, FMax, N)
trunc_exp_lower = -np.log(U_lower_bound) / rate
trunc_exp_upper = -np.log(1 - U_upper_bound) / rate

print(f"Lower bounded error: {maxT - np.max(trunc_exp_lower)}")
print(f"Upper bounded error: {maxT - np.max(trunc_exp_upper)}")

sns.histplot(x=trunc_exp_lower, stat="density")
sns.histplot(x=trunc_exp_upper, stat="density")
plt.show()

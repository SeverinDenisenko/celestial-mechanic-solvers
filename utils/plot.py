from matplotlib import pyplot as plt
import pandas as pd

df1 = pd.read_csv("result_rk.dat", sep="\s+")
df2 = pd.read_csv("result.dat", sep="\s+")

print(df1["x"][:-1])
print(df2["x"][:-1])

plt.plot(df1["x"], df1["y"])
plt.plot(df2["x"], df2["y"])
plt.show()

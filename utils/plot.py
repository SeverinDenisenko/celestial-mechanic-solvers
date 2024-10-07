from matplotlib import pyplot as plt
import pandas as pd

df = pd.read_csv("result_n_body.dat", sep="\s+")
df1 = pd.read_csv("result_rk4_n_body.dat", sep="\s+")

plt.plot(df1["x"], df1["y"])
plt.plot(df["x"], df["y"])
plt.show()

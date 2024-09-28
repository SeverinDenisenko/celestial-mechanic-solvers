from matplotlib import pyplot as plt
import pandas as pd

df = pd.read_csv("result_n_body.dat", sep="\s+")

plt.plot(df["x1"], df["y1"])
plt.plot(df["x2"], df["y2"])
plt.plot(df["x3"], df["y3"])
plt.show()

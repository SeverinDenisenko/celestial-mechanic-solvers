from matplotlib import pyplot as plt
import pandas as pd

df = pd.read_csv("result_n_body.dat", sep="\s+")

plt.plot(df["x"], df["y"])
plt.show()

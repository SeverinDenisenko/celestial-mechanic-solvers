from matplotlib import pyplot as plt
import pandas as pd

df = pd.read_csv("result.dat", sep="\s+")

print(df["x"])

plt.plot(df["x"], df["y"])
plt.show()

from matplotlib import pyplot as plt
import pandas as pd

df = pd.read_csv("newton_fractal.dat", sep="\s+")

fig, ax = plt.subplots()

ax.imshow(df)

plt.show()

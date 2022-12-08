import numpy as np
import plotly.express as px
import pandas as pd

import plotly.io as pio
pio.templates.default = "plotly_dark"

path = "build/output.csv"

df = pd.read_csv(path)
print(df.head())

new_df = pd.DataFrame(columns=["particle_id"])

# Transform the dataframe from 3 columns per particle
# into three column with the positions of all particles + one column with the particle_id
for i in range(len(df.columns)//3):
    print("processing particle ", i)
    coords = ["x", "y", "z"]
    colnames = ["px_"+str(i), "py_"+str(i), "pz_"+str(i)]
    fcolnames = colnames+["particle_id", "time"]
    df["particle_id"] = i
    particle_df = df[fcolnames].rename(columns=dict(zip(colnames, coords)))

    new_df = pd.concat([new_df, particle_df])

stepsize = int(len(new_df)//5e3)
stepsize = max(stepsize, 1)
fig1 = px.line_3d(new_df.iloc[::stepsize, :], x="x",
                y="y", z="z", color="particle_id")
# fig2 = px.scatter_3d(new_df.iloc[::10, :], x="x", y="y", z="z",
#                   animation_frame = "time", color = "particle_id",
#                   range_x = [new_df["x"].min(), new_df["x"].max()],
#                   range_y = [new_df["y"].min(), new_df["y"].max()],
#                   range_z = [new_df["z"].min(), new_df["z"].max()])

fig1.show()

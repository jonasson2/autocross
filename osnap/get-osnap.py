import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset("osnap.nc")

moc_east = ds["T_EAST"].max(dim="LEVEL")
moc_west = ds["T_WEST"].max(dim="LEVEL")

east_mean = float(moc_east.mean())
east_sd = float(moc_east.std())
west_mean = float(moc_west.mean())
west_sd = float(moc_west.std())

fig, ax = plt.subplots(figsize=(12, 5))

ax.plot(moc_east["TIME"], moc_east, color="tab:orange", lw=2,
        label=f"OSNAP East   {east_mean:.2f}±{east_sd:.2f} Sv")
ax.plot(moc_west["TIME"], moc_west, color="teal", lw=2,
        label=f"OSNAP West   {west_mean:.2f}±{west_sd:.2f} Sv")

ax.set_ylabel("MOC [Sv]")
ax.set_xlabel("")
ax.set_title("OSNAP overturning strength")
ax.grid(True, alpha=0.3)
ax.legend(frameon=False, loc="upper left")
plt.tight_layout()
plt.show()

with open("osnap-east.txt", "w") as f:
    f.write("time moc-sv\n")
    for t, x in zip(moc_east["TIME"].values, moc_east.values):
        ym = str(t)[:7]
        f.write(f"{ym} {x:.6f}\n")

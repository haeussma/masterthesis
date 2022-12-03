import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

def plot(df: pd.DataFrame, xdata: str, ydata: str, xerror: str, yerror: str, colors: str, ax: plt.Axes = None):

    if ax is None:
        ax = plt.gca()
    else:
        ax =ax

    if xerror == None:
        xerror = np.zeros(len(xdata))
    if yerror == None:
        yerror = np.zeros(len(ydata))

    ph_map = sns.color_palette("husl")
    ph_colors = [ph_map[(x)] for x in np.arange(6)]
    pH_dict = dict(zip([3.0, 3.5, 4.0, 4.5, 5.0, 5.5], ph_colors))
    pH_tuples = [pH_dict[x] for x in df["pH"].values]
    temperature_map = matplotlib.cm.get_cmap('coolwarm')
    hue = [temperature_map(x) for x in np.linspace(0,1,5)]
    color_dict = dict(zip([25, 30, 35, 40 ,45], hue))
    temperature_tuples = [color_dict[x] for x in df["temperature [C]"]]

    dict_colors = {"pH": pH_tuples, "temperature [C]":temperature_tuples}

    for i, (x, y, xerr, yerr, color, label) in enumerate(zip(df[xdata], df[ydata], df[xerror], df[yerror], dict_colors[colors], df[colors].values)):
        if colors == "pH":
            ax.scatter(x, y, color = color, label = label if not i % 5 else "")
        if colors == "temperature [C]":
            ax.scatter(x, y, color = color, label = label if i in [0,1,2,3,4] else "")
        ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt=".", ecolor=color, markersize=0)

    ax.legend(title = colors)
    ax.set_ylabel(ydata)
    ax.set_xlabel(xdata)
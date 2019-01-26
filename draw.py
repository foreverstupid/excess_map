import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import griddata

if len(sys.argv) > 1:
    dimension = sys.argv[1]
else:
    dimension = "1"

input_file_name = "surface" + dimension + "d.plt"
output_file_name = "plots" + dimension + "d.png"

# size of plot matrix
m = n = 5

# count of data values for the one plot
shift = 225

# boundary of data values
if len(sys.argv) > 3:
    zmin = float(sys.argv[2])
    zmax = float(sys.argv[3])
else:
    zmin = 0.0
    zmax = 8.0

# data columns
x_column = 2;
y_column = 3;
z_column = 4;

# smooth coefficient for color gradient
smooth_coeff = 0.05

# titles making
if len(sys.argv) > 5:
    k_min = float(sys.argv[4])
    k_max = float(sys.argv[5])
else:
    k_min = -1.0
    k_max = 1.0

km_subtitles = []
k_step = (k_max - k_min) / (m - 1)
tmp = k_min
for i in range(m):
    km_subtitles.append(r"$k_m = %5.2f$" % (tmp))
    tmp += k_step

kw_subtitles = []
k_step = (k_max - k_min) / (n - 1)
tmp = k_min
for i in range(n):
    kw_subtitles.append(r"$k_w = %5.2f$" % (tmp))
    tmp += k_step

# ticks count
xt_cnt = 5
yt_cnt = 5

# size of matrix
matrix_size = (35, 40)

# font configuration
font = { 'size': 25 }
plt.rc('font', **font)
tck_size = 20

# data getting
dat = np.genfromtxt(input_file_name, delimiter = ' ', skip_header = 0)
X_dat = dat[:,x_column]
Y_dat = dat[:,y_column]
Z_dat = dat[:,z_column]

# creating matrix
fig, axes = plt.subplots(figsize = matrix_size, nrows = m, ncols = n)
fig.suptitle(r"$N(k_m, k_w, \sigma_m, \sigma_w)$ (" + dimension + "D case)", fontsize = 35)

plot_number = 0
for ax in axes[::-1].flat:
    # convert from pandas dataframes to numpy arrays
    X, Y, Z, = np.array([]), np.array([]), np.array([])
    for i in range(plot_number * shift, plot_number * shift + shift):
        X = np.append(X, X_dat[i])
        Y = np.append(Y, Y_dat[i])
        Z = np.append(Z, Z_dat[i])

    # create x-y points to be used in heatmap
    xi = np.linspace(X.min(), X.max(), 1000)
    yi = np.linspace(Y.min(), Y.max(), 1000)

    # Z is a matrix of x-y values
    zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method = 'cubic')

    # current heatmap building
    clev = np.arange(zmin, zmax, smooth_coeff)
    hm = ax.contourf(xi, yi, zi, clev, cmap = plt.cm.inferno, vmax = zmax,
        vmin = zmin)
    ax.set_title(km_subtitles[plot_number % 5] + ", " +
        kw_subtitles[plot_number // 5], y = 1.1, fontsize = 22)
    ax.set_xlabel(r"$\sigma_m$")
    ax.set_ylabel(r"$\sigma_w$")

    # ticks settings
    ax.xaxis.set_major_locator(ticker.MaxNLocator(xt_cnt))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(xt_cnt))
    ax.tick_params(axis = "both", which = "both", labelsize = tck_size)

    plot_number = plot_number + 1

fig.subplots_adjust(right = 0.8, top = 0.9, wspace = 0.4, hspace = 0.6)
cax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

# color bar settings
cbar = fig.colorbar(hm, cax = cax)
cbar.locator = ticker.MaxNLocator(10)
cbar.formatter = ticker.FormatStrFormatter("%3.1f")
cbar.update_ticks()

fig.savefig(output_file_name)

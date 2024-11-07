

import numpy as np
import corner
import matplotlib.pyplot as plt

# Load the datasets for three different data sets (data1, data2, data3)
#data1, data2, and data3 are the posterior distribution of the EOS parameters 
flat_samples_data1 = np.loadtxt('equal_weighted_data1.txt', delimiter=' ')
flat_samples_data2 = np.loadtxt('equal_weighted_data2.txt', delimiter=' ')
flat_samples_data3 = np.loadtxt('equal_weighted_data3.txt', delimiter=' ')


# Extract only the first two columns of each dataset (assumed to be the relevant parameters)
flat_samples_subset_data1 = flat_samples_data1[:, :2]
flat_samples_subset_data2 = flat_samples_data2[:, :2]
flat_samples_subset_data3 = flat_samples_data3[:, :2]

# Define the labels for the parameters in the corner plot
labels_subset = [r"$\epsilon$", r"$n_{\rm sur}$"]

# Create the initial corner plot for data1 without automatic titles
figure2 = corner.corner(
    flat_samples_subset_data1,
    labels=labels_subset,
    smooth=1.0,
    label_kwargs=dict(fontsize=22),
    quantiles=[0.16, 0.5, 0.84],  # 16th, 50th (median), and 84th percentiles
    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),  # Contour levels
    plot_density=False,
    plot_datapoints=False,
    fill_contours=True,
    show_titles=False,  # Disable automatic titles for flexibility in custom titles
    color='steelblue'
)

# Overlay the corner plot for data2
figure2 = corner.corner(
    flat_samples_subset_data2,
    labels=labels_subset,
    smooth=1.0,
    label_kwargs=dict(fontsize=22),
    quantiles=[0.16, 0.5, 0.84],
    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
    plot_density=False,
    plot_datapoints=False,
    fill_contours=True,
    show_titles=False,
    color='darkorange',
    fig=figure2  # Add to the existing plot
)

# Overlay the corner plot for data3
figure2 = corner.corner(
    flat_samples_subset_data3,
    labels=labels_subset,
    smooth=1.0,
    label_kwargs=dict(fontsize=22),
    quantiles=[0.16, 0.5, 0.84],
    levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
    plot_density=False,
    plot_datapoints=False,
    fill_contours=True,
    show_titles=False,
    color='seagreen',
    fig=figure2  # Add to the existing plot
)

# Remove any automatic titles from diagonal axes to avoid redundancy
axes = np.array(figure2.axes).reshape((2, 2))
for ax in [axes[0, 0], axes[1, 1]]:
    for line in ax.lines:
        line.set_visible(False)

# Add custom titles (quantiles and medians) for each dataset manually on the plot
for i, (dataset, color) in enumerate(
    [(flat_samples_subset_data1, 'steelblue'),   
     (flat_samples_subset_data2, 'darkorange'), 
     (flat_samples_subset_data3, 'seagreen')]):  

    # Calculate quantiles for epsilon (first parameter)
    q16, q50, q84 = np.percentile(dataset[:, 0], [16, 50, 84])
    title_text = f"{q50:.2f}$^{{+{q84 - q50:.2f}}}_{{-{q50 - q16:.2f}}}$"
    # Place the calculated median and quantiles on the plot for the first parameter
    axes[0, 0].text(
        0.2, 1.46 - i * 0.15,  # Adjust position for each dataset's text
        f"{title_text}",
        transform=axes[0, 0].transAxes,
        color=color,
        fontsize=16,
        verticalalignment="top"
    )

# Repeat for the second parameter (n_sur)
for i, (dataset, color) in enumerate(
    [(flat_samples_subset_data1, 'steelblue'), 
     (flat_samples_subset_data2, 'darkorange'), 
     (flat_samples_subset_data3, 'seagreen')]):   

    # Calculate quantiles for n_sur (second parameter)
    q16, q50, q84 = np.percentile(dataset[:, 1], [16, 50, 84])
    title_text = f"{q50:.2f}$^{{+{q84 - q50:.2f}}}_{{-{q50 - q16:.2f}}}$"
    # Place the calculated median and quantiles on the plot for the second parameter
    axes[1, 1].text(
        0.2, 1.46 - i * 0.15,  # Adjust position for each dataset's text
        f"{title_text}",
        transform=axes[1, 1].transAxes,
        color=color,
        fontsize=16,
        verticalalignment="top"
    )

# Create a custom legend for the different datasets
lines = [plt.Line2D([0], [0], color='steelblue', lw=2.5),
         plt.Line2D([0], [0], color='darkorange', lw=2.5),
         plt.Line2D([0], [0], color='seagreen', lw=2.5)
         ]  
labels = [r"$data1$", r"$data2$", r"$data3$"]
plt.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.8, 1.8), frameon=False)

# Display and save the plot as a PDF
plt.show()
figure2.savefig("posterior_plot.pdf", bbox_inches='tight')



import numpy as np
import corner
import matplotlib.pyplot as plt

def corner_three(array1, array2, array3, params_labels=None, sample_labels=None, save_plot=False, filename="posterior_corner_plot.pdf"):
    arrays = [array1, array2, array3]
    n_parmas = array1.shape[1]

    # Create the initial corner plot for data1 without automatic titles
    figure = corner.corner(
        array1,
        labels=params_labels,
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
    figure = corner.corner(
        array2,
        labels=params_labels,
        smooth=1.0,
        label_kwargs=dict(fontsize=22),
        quantiles=[0.16, 0.5, 0.84],
        levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
        plot_density=False,
        plot_datapoints=False,
        fill_contours=True,
        show_titles=False,
        color='darkorange',
        fig=figure  # Add to the existing plot
    )

    # Overlay the corner plot for data3
    figure = corner.corner(
        array3,
        labels=params_labels,
        smooth=1.0,
        label_kwargs=dict(fontsize=22),
        quantiles=[0.16, 0.5, 0.84],
        levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
        plot_density=False,
        plot_datapoints=False,
        fill_contours=True,
        show_titles=False,
        color='seagreen',
        fig=figure  # Add to the existing plot
    )

    # Remove any automatic titles from diagonal axes to avoid redundancy
    axes = np.array(figure.axes).reshape((n_parmas, n_parmas))
    for i in range(n_parmas):
        ax = axes[i,i]
        for line in ax.lines:
            line.set_visible(False)

    # Add custom titles (quantiles and medians) for each dataset manually on the plot
        for j, (dataset, color) in enumerate(
            [(array1, 'steelblue'),   
            (array2, 'darkorange'), 
            (array3, 'seagreen')]):  

            # Calculate quantiles for epsilon (first parameter)
            q16, q50, q84 = np.percentile(dataset[:, i], [16, 50, 84])
            title_text = f"{q50:.2f}$^{{+{q84 - q50:.2f}}}_{{-{q50 - q16:.2f}}}$"
            # Place the calculated median and quantiles on the plot for the first parameter
            ax.text(
                0.2, 1.46 - j * 0.15,  # Adjust position for each dataset's text
                f"{title_text}",
                transform=ax.transAxes,
                color=color,
                fontsize=16,
                verticalalignment="top"
            )

    # Create a custom legend for the different datasets
    lines = [plt.Line2D([0], [0], color='steelblue', lw=2.5),
            plt.Line2D([0], [0], color='darkorange', lw=2.5),
            plt.Line2D([0], [0], color='seagreen', lw=2.5)
            ]  
    labels = sample_labels
    plt.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.8, 1.8), frameon=False)

    # Display and save the plot as a PDF
    if save_plot:
        figure.savefig(filename, bbox_inches='tight')
        print(f"Corner plot of posterior saved as {filename}.")
    else:
        plt.show()


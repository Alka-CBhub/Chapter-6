import matplotlib as mpl

# Plot format
def set_plot_style(dpi=500):
    """
    Call once at the top of any script to apply global rcParams.
    Pass a different dpi if you want figures at lower/higher resolution.
    """
    mpl.rcParams.update({
        # Use serif fonts throughout
        "text.usetex": False,
        "font.family": "serif",
        # Axes labels and titles
        "axes.labelsize": 16,
        "axes.titlesize": 16,
        # Tick‐label sizes
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        # Legend font size
        "legend.fontsize": 13,
        # Default line thickness and marker size
        "lines.linewidth": 2.2,
        "lines.markersize": 7,
        # Grid appearance
        "grid.linestyle": "--",
        "grid.alpha": 0.6,
        # Respect the dpi argument
        "figure.dpi": dpi
    })


# Spline format
def set_spines_black(ax):
    """
    Given a Matplotlib Axes instance, make every spine black and linewidth=2.
    """
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_edgecolor('k')

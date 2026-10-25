import matplotlib.pyplot as plt

def setup_thesis_plots():
    """
    Sets up the global matplotlib layout for a 24x17 cm thesis.
    Applies a modern style with Fira Sans / Roboto typography.
    """
    plt.rcParams.update({
        # Dimensions recalibrated for a 24x17 cm format (usable text space ~12 cm / 4.7 inches)
        "figure.figsize": (4.7, 3.2), 
        "figure.dpi": 300,
        "savefig.bbox": "tight",
        
        "font.family": "sans-serif",
        "font.sans-serif": ["Fira Sans", "Roboto", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "normal",    
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "legend.frameon": False,
        
        # Axes and grid (modern style)
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,    
        "grid.color": "#E0E0E0",
        "grid.linestyle": "--",
        "grid.linewidth": 0.8,
    })
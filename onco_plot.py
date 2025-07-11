import pcdl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load simulation frame
mcds = pcdl.TimeStep("output00000050.xml", output_path="output")
df = mcds.get_cell_df()

# Ensure required columns exist
required_cols = ["cell_type", "position_x", "position_y", "current_cycle_phase_exit_rate"]
if not all(col in df.columns for col in required_cols):
    raise ValueError(f"Missing one of the required columns: {required_cols}")

# Filter sensitive and resistant cells
sensitive = df[df["cell_type"] == "sensitive"]
resistant = df[df["cell_type"] == "resistant"]

# Function to plot by type
def plot_growth_by_type(data, title, filename):
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x="position_x",
        y="position_y",
        data=data,
        hue="current_cycle_phase_exit_rate",
        palette="viridis",
        s=10,
        linewidth=0,
        alpha=0.7
    )
    plt.title(title)
    plt.xlabel("Position X")
    plt.ylabel("Position Y")
    plt.legend(title="Growth rate", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

# Plot and save
plot_growth_by_type(sensitive, "Sensitive Cells: Spatial Growth Rate", "sensitive_growthrate.png")
plot_growth_by_type(resistant, "Resistant Cells: Spatial Growth Rate", "resistant_growthrate.png")

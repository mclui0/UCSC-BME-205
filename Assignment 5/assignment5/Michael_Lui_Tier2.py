import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.manifold import MDS

# Define CPK colors based on atomic element names
cpk_colors = {
    'H': 'white',       # Hydrogen
    'C': 'black',       # Carbon
    'N': 'blue',        # Nitrogen
    'O': 'red',         # Oxygen
    'F': 'green',       # Fluorine
    'Cl': 'green',      # Chlorine
    'Br': 'darkred',    # Bromine
    'I': 'purple',      # Iodine
    'P': 'orange',      # Phosphorus
    'S': 'yellow',      # Sulfur
    'B': 'salmon',      # Boron
    'Li': 'purple',     # Lithium
    # Add more elements as needed
}

# Define atomic weights for marker sizes
atomic_weights = {
    'H': 1, 'C': 12, 'N': 14, 'O': 16, 'F': 19, 'Cl': 35, 'Br': 80, 'I': 127,
    'P': 31, 'S': 32, 'B': 11, 'Li': 7
}

# Part 3 Extra Credit: 3D Visualization of the molecular structure
def visualize_molecule_3d():
    # Load the 3D coordinates and elements from MDS output
    coordinates_df = pd.read_csv("molecule_coordinates.csv")
    elements = coordinates_df['Element']
    coordinates = coordinates_df[['X', 'Y', 'Z']].values

    # Load molecular distance data, assuming the first row may be a header
    molecule_data = pd.read_csv("molecule_distances.tsv", sep='\t', header=0)

    # Check if the first column is 'Element' and extract distances correctly
    if molecule_data.columns[0] == 'Element':
        distances = molecule_data.iloc[:, 1:].values  # Get all distance values except the first column
    else:
        distances = molecule_data.values  # Use the entire matrix if no extra columns

    # Ensure we are working with a square matrix (25x25)
    distances = distances[:len(elements), :len(elements)]

    # Verify if the matrix is square; if not, raise an error
    if distances.shape[0] != distances.shape[1]:
        raise ValueError(f"Distance matrix is not square: shape = {distances.shape}")
    
    # Make the matrix symmetric if necessary
    if distances.shape[0] == distances.shape[1]:  # Check if it is square
        distances = (distances + distances.T) / 2  # Symmetrize
    else:
        raise ValueError(f"Cannot symmetrize distance matrix of shape: {distances.shape}")

    # Create 3D plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot each atom as a scatter point
    for i, (x, y, z) in enumerate(coordinates):
        element = elements[i]
        color = cpk_colors.get(element, 'gray')  # Default color if element not in CPK
        size = atomic_weights.get(element, 10) * 5  # Scale marker size
        ax.scatter(x, y, z, color=color, s=size, alpha=0.8)

    # Draw bonds as line segments for pairs with distance < 1.6
    threshold_distance = 1.6
    for i in range(len(coordinates)):
        for j in range(i + 1, len(coordinates)):
            if distances[i, j] < threshold_distance:
                ax.plot([coordinates[i, 0], coordinates[j, 0]],
                        [coordinates[i, 1], coordinates[j, 1]],
                        [coordinates[i, 2], coordinates[j, 2]], color='grey', alpha=0.5)

    # Set plot labels and save the figure
    ax.set_title("3D Visualization of Molecular Structure")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.savefig("molecule_3D_plot.png")

# Call the function to generate and save the plot
visualize_molecule_3d()

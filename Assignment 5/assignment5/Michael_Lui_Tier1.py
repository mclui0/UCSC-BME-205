import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from mpl_toolkits.mplot3d import Axes3D

# Part 1: PCA on MNIST dataset
def pca_mnist():
    # Load the MNIST subset
    MNIST_X = np.load("MNIST_X_subset.npy", allow_pickle=True)  # Shape (6000, 784)
    MNIST_y = np.load("MNIST_y_subset.npy", allow_pickle=True)  # Shape (6000,)
    
    # Step 1: Perform PCA to reduce dimensions to 2
    pca = PCA(n_components=2, random_state=42)
    MNIST_X_pca = pca.fit_transform(MNIST_X)  # Transformed shape (6000, 2)

    # Step 2: Visualize the 2D PCA-transformed data
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(MNIST_X_pca[:, 0], MNIST_X_pca[:, 1], c=MNIST_y, cmap='tab10', alpha=0.6)
    plt.colorbar(scatter, label="Digit Label")
    plt.title("2D PCA of MNIST Data")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.savefig("MNIST_PCA_2D.png")  # Save the plot
    plt.close()
    
    # Step 3: Reconstruct an image using the first 2 principal components
    MNIST_X_reduced = pca.transform([MNIST_X[0]])  # Take the first image, shape (1, 2)
    MNIST_X_reconstructed = pca.inverse_transform(MNIST_X_reduced)  # Reconstruct, shape (1, 784)
    
    # Save the original and reconstructed images
    plt.imsave("MNIST_original.png", MNIST_X[0].reshape(28, 28), cmap='gray')  # Original image
    plt.imsave("MNIST_reconstructed_2PC.png", MNIST_X_reconstructed.reshape(28, 28), cmap='gray')  # Reconstructed image
    
    # Step 4: Reconstruct an image from a selected 2D point
    # Example coordinates for digit "1"
    chosen_point = np.array([-1.5, 2])  # Example point in PCA space
    point_reconstructed = pca.inverse_transform([chosen_point])  # Reconstruct the image from the chosen point
    plt.imsave("MNIST_reconstructed_1_from_coord.png", point_reconstructed.reshape(28, 28), cmap='gray')  # Save the reconstructed image

# Part 2: PCA on Dogs SNP dataset
def pca_dogs_snp():
    # Load Dogs SNP data
    dogs_X = np.load("dogs_X.npy", allow_pickle=True)
    dogs_clades = np.load("dogs_clades.npy", allow_pickle=True)
    
    # Convert clade labels to numeric values
    unique_clades = np.unique(dogs_clades)
    clade_to_numeric = {clade: idx for idx, clade in enumerate(unique_clades)}
    numeric_clades = np.array([clade_to_numeric[clade] for clade in dogs_clades])
    
    # Perform PCA to reduce to 2D
    pca = PCA(n_components=2, random_state=42)
    dogs_X_pca = pca.fit_transform(dogs_X)
    
    # Save 2D PCA scatter plot with color-coded clades
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(dogs_X_pca[:, 0], dogs_X_pca[:, 1], c=numeric_clades, cmap='viridis', alpha=0.6)
    colorbar = plt.colorbar(scatter, ticks=range(len(unique_clades)))
    colorbar.set_ticklabels(unique_clades)
    colorbar.set_label("Clade Label")
    plt.title("2D PCA of Dogs SNP Data")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.savefig("Dogs_PCA_2D.png")
    plt.close()

# Part 3: MDS on molecular distance matrix
def mds_molecule():
    # Load the data
    molecule_data = pd.read_csv("molecule_distances.tsv", sep='\t', header=0)
    
    # # Ensure the DataFrame has the expected shape
    # print(f"Original DataFrame shape: {molecule_data.shape}")
    
    # Assuming the first column is atom index and the second column is element
    distances = molecule_data.iloc[:, 2:].values  # Keep only the numeric distance values
    elements = molecule_data.iloc[:, 1].values     # Keep the element names from the second column
    
    # # Check the shape of the distance matrix and elements
    # print(f"Shape of distance matrix: {distances.shape}")
    # print(f"Number of elements: {len(elements)}")

    # Ensure the distance matrix is square
    if distances.shape[0] != distances.shape[1]:
        raise ValueError(f"Distance matrix is not square: shape = {distances.shape}")
    
    # Symmetrize the matrix if necessary
    distances = (distances + distances.T) / 2

    # Perform MDS to reconstruct 3D coordinates
    mds = MDS(n_components=3, dissimilarity='precomputed', random_state=42)
    molecule_coords = mds.fit_transform(distances)
    
    # Save coordinates to CSV with 'Element' as the first column
    coordinates_df = pd.DataFrame(molecule_coords, columns=['X', 'Y', 'Z'])
    coordinates_df.insert(0, 'Element', elements)  # Insert 'Element' as the first column
    coordinates_df.to_csv("molecule_coordinates.csv", index=False)
    
    # Plot 3D scatter of molecule with element types
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(molecule_coords[:, 0], molecule_coords[:, 1], molecule_coords[:, 2],
                         c=pd.factorize(elements)[0], cmap='cool', alpha=0.8)
    ax.set_title("3D MDS of Molecular Structure")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.colorbar(scatter, label="Element Type")
    plt.savefig("Molecule_MDS_3D.png")
    plt.close()

# Run all parts
if __name__ == "__main__":
    pca_mnist()
    pca_dogs_snp()
    mds_molecule()

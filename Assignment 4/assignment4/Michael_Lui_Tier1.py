import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import confusion_matrix
import pandas as pd

# Load MNIST data
MNIST_X = np.load('MNIST_X_subset.npy', allow_pickle=True)
MNIST_y = np.load('MNIST_y_subset.npy', allow_pickle=True)

# Reshape and visualize the first example
first_example = MNIST_X[0].reshape(28, 28)
plt.imshow(first_example, cmap='gray')
plt.title("First MNIST Example")
plt.axis('off')
plt.show()

# Perform K-Means clustering on MNIST dataset
def perform_kmeans(X, y, K):
    kmeans = KMeans(n_clusters=K, random_state=42)
    y_pred = kmeans.fit_predict(X)

    # Reshape and visualize centroids
    centroids = kmeans.cluster_centers_.reshape(K, 28, 28)
    rows = K // 2 + K % 2  # Calculate the number of rows needed
    plt.figure(figsize=(10, 5))
    for i in range(K):
        plt.subplot(rows, 2, i + 1)
        plt.imshow(centroids[i], cmap='gray')
        plt.axis('off')
    plt.suptitle(f"Centroids for K={K}")
    plt.savefig(f'centroids_k{K}.png')
    plt.show()

    # Compute clustering error
    clustering_error = 0
    for cluster in range(K):
        cluster_labels = y[y_pred == cluster]
        if len(cluster_labels) > 0:
            majority_label = np.bincount(cluster_labels).argmax()
            error_count = len(cluster_labels) - np.sum(cluster_labels == majority_label)
            clustering_error += error_count
    
    return clustering_error

# Perform K-Means for K=10 and K=11 and report errors
for K in [10, 11]:
    error = perform_kmeans(MNIST_X, MNIST_y, K)
    print(f'K={K} Error={error}')

# Load Dogs SNP Dataset (Tier 3)
dogs_X = np.load('dogs_X.npy', allow_pickle=True)
dogs_clades = np.load('dogs_clades.npy', allow_pickle=True)

# Optional: You can apply a similar K-Means clustering process for the Dogs SNP dataset.

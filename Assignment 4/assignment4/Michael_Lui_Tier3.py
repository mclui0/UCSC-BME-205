import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from collections import Counter

# Load Dogs SNP data and clade information
dogs_X = np.load('dogs_X.npy', allow_pickle=True)
dogs_clades = np.load('dogs_clades.npy', allow_pickle=True)

# Perform hierarchical clustering on the Dogs dataset
Z = linkage(dogs_X, method='average')

# Assign cluster labels for K=30
cluster_labels = fcluster(Z, t=30, criterion='maxclust')

# Visualize the dendrogram
plt.figure(figsize=(10, 7))
dendro = dendrogram(Z, truncate_mode='lastp', p=30, show_leaf_counts=True)
plt.title('Hierarchical Clustering Dendrogram (K=30)')
plt.xlabel('Sample Index')
plt.ylabel('Distance')
plt.savefig('Dogs_dendrogram.png')
plt.show()

# Find the majority clade for each cluster and compute clustering error
def compute_clustering_error(cluster_labels, true_labels, num_clusters):
    clustering_error = 0
    for cluster in range(1, num_clusters + 1):
        cluster_indices = np.where(cluster_labels == cluster)[0]
        if len(cluster_indices) > 0:
            true_cluster_labels = true_labels[cluster_indices]
            majority_clade = Counter(true_cluster_labels).most_common(1)[0][0]
            error_count = len(true_cluster_labels) - np.sum(true_cluster_labels == majority_clade)
            clustering_error += error_count
    return clustering_error

# Calculate clustering error
error = compute_clustering_error(cluster_labels, dogs_clades, 30)
print(f'Clustering Error: {error}')

# Label each terminal node in the dendrogram with the most common clade
# (This part is implicitly handled in the visualization process, assuming manual inspection)

### 1. **K-Means Clustering on MNIST (Tier 1)**:
- **Function**: This script performs K-Means clustering on the MNIST dataset (images of handwritten digits). K-Means tries to group the data into 10 or 11 clusters (representing digits), then calculates how often it groups the images incorrectly.
- **Steps**:
  1. Load a small subset of MNIST data (6000 images).
  2. Perform K-Means clustering for `K=10` and `K=11`.
  3. Visualize the cluster "centroids" (average of each group) as images, and save them.
  4. Calculate how many images in each cluster don't match the most common digit in that cluster (this is called the clustering error).
- **Main takeaway**: The script groups the images based on pixel similarity, then calculates how accurate those groups are.

### 2. **Hierarchical Clustering on MNIST (Tier 2 Extra Credit)**:
- **Function**: This script uses hierarchical clustering to group the MNIST images. Instead of a fixed number of clusters, hierarchical clustering builds a "tree" of clusters where smaller clusters merge into larger ones, and you can cut the tree at different levels (here, we cut it at 10 clusters).
- **Steps**:
  1. Load the same MNIST data as before.
  2. Perform hierarchical clustering using "average linkage," which groups images based on average similarity.
  3. Display a dendrogram (a tree diagram showing how clusters merge) and label the 10 terminal clusters with the most common digit in each cluster.
  4. Calculate the clustering error (same method as in the K-Means script).
- **Main takeaway**: Instead of grouping the images directly into 10 or 11 clusters, hierarchical clustering shows a gradual merging process. The dendrogram helps visualize how the images are grouped at each stage.

### 3. **Hierarchical Clustering on Dogs SNP Dataset (Tier 3 Extra Credit)**:
- **Function**: This script performs hierarchical clustering on genetic data (SNPs) from 1355 dog samples. The goal is to group the dogs based on genetic similarity and compare the groups (clusters) with known clade information (genetic groups of breeds).
- **Steps**:
  1. Load the genetic data for the dogs.
  2. Perform hierarchical clustering using average linkage.
  3. Visualize the dendrogram with 30 terminal clusters, and label each cluster with the most common clade.
  4. Calculate the clustering error (same method as in the previous scripts).
- **Main takeaway**: The script groups dog samples based on their genetic features, then checks how well these groups align with known breed groupings (clades).

### **Differences between the scripts**:
- **Type of data**:
  - The first two scripts work with image data (MNIST dataset of handwritten digits).
  - The third script works with genetic data (SNPs from dog samples).
  
- **Clustering method**:
  - The first script uses **K-Means clustering**, where the number of clusters (K) is fixed in advance.
  - The second and third scripts use **hierarchical clustering**, which doesn't need a predefined number of clusters and shows a tree-like structure of how the data can be grouped at different levels.

- **Visualization**:
  - The first script visualizes cluster centroids as images.
  - The second and third scripts visualize dendrograms, which show the hierarchical structure of the clusters.

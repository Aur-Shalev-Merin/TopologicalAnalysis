import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# Load the CSV file
data = pd.read_csv("examples/centers_of_mass.csv")

# Make sure that you have the columns 'x', 'y', and 'z' written in the data
X = data[['x', 'y', 'z']].values

# num of clusters
k = 8

# Create a KMeans instance
kmeans = KMeans(n_clusters=k, random_state=0)

# Fit the model to the data
kmeans.fit(X)

# Get the cluster labels for each point
labels = kmeans.labels_

# Get the coordinates of the cluster centers
centers = kmeans.cluster_centers_

# Visualize the clusters in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the data points, colored by cluster
ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=labels)

# Plot the cluster centers
ax.scatter(centers[:, 0], centers[:, 1], centers[:, 2], marker='x', s=200, c='red')

plt.show()
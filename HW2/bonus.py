import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.cluster import KMeans


# Load the iris dataset
iris = load_iris()
data_points = iris.data

#the Ks are 1 to 10
K_range = range(1, 11)

# Run KMeans for k = 1 to 10, calculate the inertia with Kmeans and save it
inertias = []
for k in K_range:
    kmeans = KMeans(n_clusters=k, init='k-means++', random_state=0)
    kmeans.fit(data_points)
    inertias.append(kmeans.inertia_)

# Plot the elbow curve
plt.figure(figsize=(8, 5))
plt.plot(K_range, inertias, marker='o', linestyle='-')
plt.title("Elbow Method for selection of optimal 'K' clusters")
plt.xlabel("K")
plt.ylabel("Inertia")

# decide which k is the elbow point manually (based on visual inspection)
elbow_k = 3

# "print" arrow on the plot that point to the elbow point
elbow_x = elbow_k
elbow_y = inertias[elbow_k - 1]
plt.annotate("Elbow Point",
             xy=(elbow_x, elbow_y),
             xytext=(elbow_x + 1, elbow_y + 100),
             arrowprops=dict(facecolor='black', shrink=0.05),
             bbox=dict(boxstyle="round", fc="w"))


# Save plot
plt.tight_layout()
plt.savefig("elbow.png")
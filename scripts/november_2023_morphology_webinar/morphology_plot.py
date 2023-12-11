import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

"""
Script to create a 3D plot of a crystal morphology.
"""


def generate_morphology_plot(morphology, labels=None):
    fig = plt.figure(figsize=(10., 10.))  # 3D graph instance
    ax = fig.add_subplot(111, projection='3d')  # 3D Axes
    _ = ax.set_xlim3d(-15, 15)  # Set the axes limits
    _ = ax.set_ylim3d(-15, 15)
    _ = ax.set_zlim3d(-15, 15)
    _ = ax.grid(False)  # To hide the gridlines
    _ = plt.axis('off')  # To hide the axes
    for i, facet in enumerate(morphology.facets):
        for edge in facet.edges:
            Axes3D.plot(ax,
                        [coord[0] for coord in edge],  # The x coordinates
                        [coord[1] for coord in edge],  # The y coordinates
                        [coord[2] for coord in edge],  # The z coordinates
                        c='black',
                        linewidth=1.5)
        for edge in facet.edges:
            vertices = [(edge[0], edge[1], facet.centre_of_geometry)]
            Axes3D.add_collection3d(ax, Poly3DCollection(vertices, color='blue', linewidth=0, alpha=0.3))
        if labels:
            ax.text(facet.centre_of_geometry[0], facet.centre_of_geometry[1], facet.centre_of_geometry[2],
                    '''   {}
                {}'''.format(facet.miller_indices.hkl, labels[i]),
                    color='black')
        else:
            ax.text(facet.centre_of_geometry[0], facet.centre_of_geometry[1], facet.centre_of_geometry[2],
                    '''   {}
            {}'''.format(facet.miller_indices.hkl, round(morphology.relative_area(facet.miller_indices), 3)),
                    color='black')

        ax.scatter(facet.centre_of_geometry[0],
                   facet.centre_of_geometry[1],
                   facet.centre_of_geometry[2],
                   s=10,
                   color='black')
    plt.show()

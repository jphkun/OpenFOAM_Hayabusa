#                         #p6
#
#
#
#
#
#                         #p5
#        #p3            #    #                p8
#                     #        #p7
#           b1      #p4   y    #
#                  #      |    #
#     #p2          #p1    #->x #p10            p9

import matplotlib.pyplot as plt # type: ignore
import numpy as np              # type: ignore
import sys                      # type: ignore
import os                       # type: ignore
from scipy.spatial.transform import Rotation as R                   # type: ignore[code, ...]
from classy_blocks.classes.mesh import Mesh                         # type: ignore[code, ...]
from classy_blocks.classes.operations import Extrude, Face, Revolve # type: ignore[code, ...]
from classy_blocks.classes import operations                        # type: ignore[code, ...]
from classy_blocks.classes.block import Block                       # type: ignore[code, ...]
from classy_blocks.classes.primitives import Edge                   # type: ignore[code, ...]
from classy_blocks.util import functions as f                       # type: ignore[code, ...]

# OKTODO: Implement the cell inflation on the 3D mesh
# TODO: Implement the 2D mesh
# TODO: Implement a simple case with rhoSonicFoam
# TODO: Add mesh rotation for the 2D case
# TODO: Add mesh rotation for the 3D case
# TODO: README.md
# TODO: Drawing for the docuementation
# TODO: Correct the arc center problem
# TODO: Equilibrate the number of cell for each block
# TODO: BUG, solve the problem with -t =/= -p

# For the sake of doing it:
# TODO: Finish tests for point_on_ellipsoid
# TODO: point_on_straigt_front:
#       - add test
#       - add special case when x < 0.11
# TODO: point_on_straigt_back:
#       - add test
#       - add special case when x > 0.8

class CapsuleMesh:
    def __init__(self,
                 dimensions: str,
                 capsuleDia: float,
                 rNpoints: int,
                 thetaNpoints: int,
                 phiNpoints: int,
                 inflation: float) -> None:
        """
        Init function

        Parameters
        ----------
        dimensions: str
            Number of dimension in the mesh. could be 2D or 3D
        capsuleDia: float
            Capsules diameter, determines the capsule diameter and the mesh size
        rNpoints: int
            Number of mesh points in the "radial" direction
        thetaNpoints: int
            Number of mesh points in the "theta" direction
        phiNpoints: int
            Number of mesh points in the "phi" direction
        inflation: float
            Inflation layer parameter
        """
        # General mesh parameters
        self.dimensions = dimensions
        self.capsuleDia = 0.4
        self.capsuleRad = self.capsuleDia/2
        self.capsuleSphereCenter = [self.capsuleDia/5, 0, 0]
        self.rNpoints = rNpoints
        self.thetaNpoints = thetaNpoints
        self.phiNpoints = phiNpoints
        self.cellInflationRatio = inflation

        # a1, b1, c1, are the "front" ellipsoid semi-axis
        self.a1 = 10 * self.capsuleDia
        self.b1 = 10 * self.capsuleDia
        self.c1 = 10 * self.capsuleDia

        # a2, b2, c2, are the "back" ellipsoid semi-axis
        self.a2 = 4  * self.capsuleDia
        self.b2 = 10 * self.capsuleDia
        self.c2 = 10 * self.capsuleDia

    def point_on_sphere(self,
                        rho: float,
                        theta: float,
                        phi: float,
                        origin=[0, 0, 0]) -> list[float]:
        """
        Returns a point [x,y,z] array of the asked point on the sphere

        Parameters
        ----------
        rho: float
            radius of the sphere
        theta: float
            angle betweent eh axis X and the point on the XY plane
        phi: float
            angle between the axis Z and the point

        Return
        ------
        point: vector
            [x,y,z] corrdinates from the computed point

        """
        # Tests if phi is between 0 and 2pi
        if theta < 0 or theta > 2*np.pi:
            print("Theta is not between 0 and 2pi")
            sys.exit()

        # Tests if theta is between 0 and pi
        if phi < 0 or phi > np.pi:
            print("Phi is not between 0 and pi")
            sys.exit()

        x = origin[0] + rho * np.cos(theta) * np.sin(phi)
        y = origin[1] + rho * np.sin(theta) * np.sin(phi)
        z = origin[2] + rho * np.cos(phi)
        point = [x, y, z]
        return point

    def point_on_ellipsoid(self, theta: float, phi: float) -> list[float]:
        """
        Returns a surface point of the ellipsoid that constitutes the O mesh.

        Parameters
        ----------
        theta: float
            angle between the reference axis X and the point in the "XY" plane
        phi: float
            angle between the reference axis Z and the point

        Returns
        -------
        point: Vector
            X, Y, Z position on the surface of the ellipsoid
        """
        # Tests if phi is between 0 and 2pi
        if theta < 0 or theta > 2*np.pi:
            print("Theta is not between 0 and 2pi")
            sys.exit()

        # Tests if theta is between 0 and pi
        if phi < 0 or phi > np.pi:
            print("Phi is not between 0 and pi")
            sys.exit()

        x = 0
        y = 0
        z = 0
        # Front ellipsoid
        if theta <= np.pi/2 or theta > (3/2)*np.pi:
            cadran = 'Front'  # Used in testing
            x = self.a1 * np.cos(theta) * np.sin(phi)
            y = self.b1 * np.sin(theta) * np.sin(phi)
            z = self.c1 * np.cos(phi)
        # Back ellipsoid
        elif theta > np.pi/2 or theta >= (3/2)*np.pi:
            cadran = 'Back'  # Used in testing
            x = self.a2 * np.cos(theta) * np.sin(phi)
            y = self.b2 * np.sin(theta) * np.sin(phi)
            z = self.c2 * np.cos(phi)
        else:
            msg1 = 'Point_on_ellipsoid: Error with phi and theta, they are not'
            msg2 = ' in the correct range'
            print(msg1 + msg2)
            sys.exit()

        point: list[float] = [x, y, z]
        return point  # cadran

    def angle_conversion(self, alpha:float) -> float:
        """
        converts the angle form one reference frame to the other

        Parameter
        ---------
        alpha: float
            angle in the first/initial reference frame

        Return
        ------
        theta: float
            angle in the second/final reference frame
        """
        b = self.capsuleDia / 2
        c = self.capsuleSphereCenter[0]
        a = np.sqrt(b**2 + c**2 - 2*b*c*np.cos(alpha))
        theta = np.arcsin(b * np.sin(alpha) / a)
        return theta

    def init_3d_plot(self):
        """
        WARNING: Do no delete useful for debugging without using paraview!

        Initialises a plot instance. This permits to plot multtiple points
        accross multiple functions.
        """
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(projection='3d')

    def plot_points(self, pnts:list, c:str):
        """
        WARNING: Do no delete useful for debugging without using paraview!

        Plots the points with a color. Points are in list format

        Parameters
        ----------
        pnts: list
            Point in 3D space so something like [[x1, y1, z1], ...]
        c: str
            Color of the point. Some examples: 'r', 'g', 'b', 'k', 'c', 'm'
        """
        for pnt in pnts:
            # print(f'plot pnt: {pnt}')
            self.ax.scatter(pnt[0], pnt[1], pnt[2], color=c)

        # To make the graph kind of equal
        dst = 3.5
        N = 25
        x = np.linspace(0, dst, N)
        self.ax.scatter( dst, 0.0, 0.0, color='w')
        self.ax.scatter(-dst, 0.0, 0.0, color='w')
        self.ax.scatter( 0.0, dst, 0.0, color='w')
        self.ax.scatter( 0.0,-dst, 0.0, color='w')
        self.ax.scatter( 0.0, 0.0, dst, color='w')
        self.ax.scatter( 0.0, 0.0,-dst, color='w')
        self.ax.scatter( x, np.zeros(N), np.zeros(N), color='k')

    def show_plot(self):
        """
        WARNING: Do no delete useful for debugging without using paraview!

        Shows the plot when python code is launched through terminal
        """
        plt.show()

    def point_on_straight_front(self, x:float):
        """
        Function that permits to compute the points at the front of the capsule
        between the arc edge and the middle YZ plane centered at the origin.
        The input is the x position of the desired point. Will be generally used
        with a linspace array.

        Parameter
        ---------
        x: float
            X position of the computed point.
        """
        # TODO: Add test for x < -0.11 and x > 0
        if x > 0:
            print('Error in point_on_straight_front, x is too big')
        x = x
        y = 0
        z = self.capsuleRad + x
        point = [x, y, z]
        return point

    def point_on_straight_back(self, x:float):
        """
        Function that permits to compute the points at the back of the capsule.
        The function is valit from theYZ plane centered at the origin, up until
        the YZ plane 1/5 of the capsule diamter on the back

        Parameter
        ---------
        x: float
            X position of the computed point.
        """
        # TODO: Generalise the 0.08
        if x > 0.08:
            print('Error in the point_on_straight_back function, x is too big')
            sys.exit()
        elif x < 0:
            print('Error in the point_on_straight_back function, x is too small')
        x = x
        y = 0
        z = self.capsuleRad - x
        point = [x, y, z]
        return point

    def point_on_vertical_back(self, z:float) -> list[float]:
        """
        This function converts the vertical point z which is inserted as an
        input, to a 3D point.

        Parameters
        ----------
        z: float
            Vertical position of the point on the back part of the capsule.

        Return
        ------
        point: list[float]
            3D point of the back of the capusle
        """
        x = self.capsuleSphereCenter[0]
        y = 0
        z = z
        point = [x, y, z]
        return point

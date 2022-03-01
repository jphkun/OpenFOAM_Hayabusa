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

class Hayabusa:

    def __init__(self,
                 dimentions: str,
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
        self.dimensions = dimentions
        self.capsuleDia = 0.4
        self.capsuleRad = self.capsuleDia/2
        self.capsuleSphereCenter = [self.capsuleDia/5, 0, 0]
        self.rNpoints = rNpoints
        self.thetaNpoints = thetaNpoints
        self.phiNpoints = phiNpoints
        self.inflatio = inflation
        # a1, b1, c1, are the "front" ellipsoid semi-axis
        self.a1 = 10 * self.capsuleDia
        self.b1 = 10 * self.capsuleDia
        self.c1 = 10 * self.capsuleDia

        # a2, b2, c2, are the "back" ellipsoid semi-axis
        self.a2 = 4  * self.capsuleDia
        self.b2 = 10 * self.capsuleDia
        self.c2 = 10 * self.capsuleDia

    # TODO: Wrong type hinting output
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

    # TODO: Wrong type hinting output
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
        # TODO: add a test
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

    def circle_on_sphere(self,
                         rho: float,
                         alpha: float,
                         beta: float,
                         gamma: float,
                         t: float
                         ) -> list[float]:
        """
        https://math.stackexchange.com/questions/643130/circle-on-sphere
        https://mathworld.wolfram.com/SphericalCoordinates.html

        Parameters
        ----------
        rho: float
            radius of the sphere
        alpha:
            radius angle (cone aperture)
        beta: float
            angle between the north pole and the central point of the circle
        gamma: float
            angle from the x-axis on the XY plane
        t: float
            angle on the cricle

        Returns
        -------
        point: vector
            (x,y,z) coordinates
        """
        x = rho * ( (np.sin(alpha) * np.cos(beta) * np.cos(gamma) * np.cos(t)) \
                  - (np.sin(alpha) * np.sin(gamma)* np.sin(t)    ) \
                  + (np.cos(alpha) * np.sin(beta) * np.cos(gamma)) )

        y = rho * ( (np.sin(alpha) * np.cos(beta) * np.sin(gamma) * np.cos(t)) \
                  + (np.sin(alpha) * np.cos(gamma)* np.sin(t)) \
                  + (np.cos(alpha) * np.sin(beta) * np.sin(gamma)))

        z = rho * (-(np.sin(alpha)*np.sin(beta)*np.cos(t)) + (np.cos(alpha)*np.cos(beta)))

        point = [x, y, z]
        return point

    # TODO: Documentation
    # TODO: Type defintion
    def point_circle_xaxis(self, og, r, theta):
        """
        Computes a point on a circle whoes central point is on the x axis

        Parameters
        ----------
        og: float
            origin point of the cirlce

        """
        x = og[0]
        y = og[1] + r * np.sin(theta)
        z = og[2] + r * np.cos(theta)
        point = [x, y, z]
        return point

    # TODO: Generaliste the a, b, c variables
    # TODO: set them as inputs?
    @staticmethod
    def angle_conversion(alpha):
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
        # TODO: Generalize
        b = 200
        c = 80
        a = np.sqrt(b**2 + c**2 - 2*b*c*np.cos(alpha))
        theta = np.arcsin(b * np.sin(alpha) / a)
        return theta

    # TODO: make the return block documentation correct
    # TODO: Documentation, look at how to do a class
    # TODO: max 100 lines for this function
    def block1_3d(self):
        """
        Constructs block 1 and return a block with edges.
        """
        # Center point of the sphere
        sphereOg = [0.08, 0, 0]

        #######################################################################
        # WARNING all the following points are in the negative x direction
        # hence add 180 deg for the correct point description in spherical
        # coordinates.
        #######################################################################
        # Angle between the XY plane the sphere origin and vertex 1 of the block
        block1Alpha = np.arcsin(1/4)
        print(f'block1Alpha [deg]: {np.rad2deg(block1Alpha)}')
        # Angle between the XZ plane the sphere origin and vertex 1 of the block
        block1Beta  = np.arcsin(1/4)
        print(f'block1Beta [deg]: {np.rad2deg(block1Beta)}')
        # Angle between the XY plane the MESH origin and vertex 4 of the block
        block1Theta = self.angle_conversion(block1Alpha)
        print(f'block1Theta [deg]: {np.rad2deg(block1Theta)}')
        # Angle between the XZ plane the MESH origin and vertex 4 of the block
        block1Phi   = self.angle_conversion(block1Beta)
        print(f'block1Phi [deg]: {np.rad2deg(block1Phi)}')

        # Vertex 0
        alpha0 = np.pi   + block1Alpha
        beta0  = np.pi/2 - block1Beta
        blk1v0 = self.point_on_sphere(self.capsuleRad, alpha0, beta0, sphereOg)
        print(f'alpha0 [deg]: {np.rad2deg(alpha0)}')
        print(f'beta0 [deg]: {np.rad2deg(beta0)}')
        # Vertex 1
        alpha1 = np.pi   - block1Alpha
        beta1  = np.pi/2 - block1Beta
        blk1v1 = self.point_on_sphere(self.capsuleRad, alpha1, beta1, sphereOg)
        # Vertex 2
        alpha2 = np.pi   - block1Alpha
        beta2  = np.pi/2 + block1Beta
        blk1v2 = self.point_on_sphere(self.capsuleRad, alpha2, beta2, sphereOg)
        # Vertex 3
        alpha3 = np.pi   + block1Alpha
        beta3  = np.pi/2 + block1Beta
        blk1v3 = self.point_on_sphere(self.capsuleRad, alpha3, beta3, sphereOg)
        # Vertex 4
        theta4 = np.pi   + block1Theta
        phi4   = np.pi/2 - block1Phi
        blk1v4 = self.point_on_ellipsoid(theta4, phi4)
        # Vertex 5
        theta5 = np.pi   - block1Theta
        phi5   = np.pi/2 - block1Phi
        blk1v5 = self.point_on_ellipsoid(theta5, phi5)
        # Vertex 6
        theta6 = np.pi   - block1Theta
        phi6   = np.pi/2 + block1Phi
        blk1v6 = self.point_on_ellipsoid(theta6, phi6)
        # Vertex 7
        theta7 = np.pi   + block1Theta
        phi7   = np.pi/2 + block1Phi
        blk1v7 = self.point_on_ellipsoid(theta7, phi7)

        # Block point definition
        self.block1Vertices = [
            blk1v0, blk1v1, blk1v2, blk1v3,  # Face 1
            blk1v4, blk1v5, blk1v6, blk1v7,  # Face 2?
        ]

        # Arc defintions
        #######################################################################
        # WARNING: alpha is set to pi since the front of the capsule is set in
        # the negative x direction
        #######################################################################
        edge01pnt = self.point_on_sphere(self.capsuleRad, np.pi,  beta1,   sphereOg)
        edge12pnt = self.point_on_sphere(self.capsuleRad, alpha1, np.pi/2, sphereOg)
        edge23pnt = self.point_on_sphere(self.capsuleRad, np.pi,  beta3,   sphereOg)
        edge30pnt = self.point_on_sphere(self.capsuleRad, alpha3, np.pi/2, sphereOg)

        # Spline defintions
        # Number of points for each spline defintion
        N = 10
        # Computation of the angles at which to compute the interpolation point
        # for the spline.
        angles45 = np.linspace(theta4, theta5, N),
        edge45pnts = [self.point_on_ellipsoid(i, phi4) for i in angles45[0]]

        angles56 = np.linspace(phi5, phi6, N),
        edge56pnts = [self.point_on_ellipsoid(theta5, i) for i in angles56[0]]

        angles67 = np.linspace(theta6, theta7, N),
        edge67pnts = [self.point_on_ellipsoid(i, phi6) for i in angles67[0]]

        angles74 = np.linspace(phi7, phi4, N),
        edge74pnts = [self.point_on_ellipsoid(theta7, i) for i in angles74[0]]

        self.block1Edges = [
            Edge(0, 1, edge01pnt),  # Arc edge
            Edge(1, 2, edge12pnt),  # Arc edge
            Edge(2, 3, edge23pnt),  # Arc edge
            Edge(3, 0, edge30pnt),  # Arc edge
            Edge(4, 5, edge45pnts),  # Spline edge
            Edge(5, 6, edge56pnts),  # Spline edge
            Edge(6, 7, edge67pnts),  # Spline edge
            Edge(7, 4, edge74pnts),  # Spline edge
        ]

        self.block1 = Block.create_from_points(
            self.block1Vertices,
            self.block1Edges
        )
        self.block1.set_patch('top', 'inlet')
        self.block1.set_patch('bottom', 'wall')

        self.block1.chop(0, count=10, c2c_expansion=1)
        self.block1.chop(1, count=10, c2c_expansion=1)
        self.block1.chop(2, count=10, c2c_expansion=1)

    def init_3d_plot(self):
        """
        Initialises a plot instance. This permits to plot multtiple points
        accross multiple functions
        """
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(projection='3d')

    def plot_points(self, pnts, c):
        """
        Plots the points with a color

        Parameters:
        pnts: array
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
        Shows the plot when python code is launched through terminal
        """
        plt.show()

    # TODO: Documentation
    def point_on_straight_front(self, x):
        """
        """
        # TODO: Add test for x < -0.11 and x > 0
        if x > 0:
            print('Error in point_on_straight_front, x is too big')
        x = x
        y = 0
        z = self.capsuleRad + x
        point = [x, y, z]
        return point

    # TODO: Documentation
    def point_on_straight_back(self, x):
        """
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

    # TODO: Documentation + type defining
    def point_on_vertical_back(self, z: float):
        """
        This function converts the vertical point z which is inserted as an
        input, to a 3D point.

        Parameters
        ----------
        z: float
            Vertical position of the point on the back part of the capsule.

        Return
        ------
        point: array
            3D point of the back of the capusle
        """
        x = self.capsuleSphereCenter[0]
        y = 0
        z = z
        point = [x, y, z]
        return point

    # TODO: Documentation
    # TODO: Logic semplification
    def section_points(self):
        """

        """
        # TODO: Generalize the sphere origin
        # tODO: Clarify this mess
        sphereOg = [0.08, 0, 0]

        #######################################################################
        # WARNING all the following points are in the negative x direction
        # hence add 180 deg for the correct point description in spherical
        # coordinates.
        #######################################################################

        # Section 1 angles (alpha) bottom (bot)
        # Angle between the XY plane the sphere origin and vertex 3 of the block 1
        section1AlphaBot1 = np.arcsin(1/4)
        section1AlphaBot2 = np.pi + section1AlphaBot1

        # Angle between the XZ plane the sphere origin and vertex 3 of the block 1
        section1BetaBot1 = np.arcsin(1/4)
        section1BetaBot2 = np.pi/2 - section1BetaBot1

        # Inlet points
        # Angle between the XY plane the MESH origin and vertex 7? 4? of the block 1
        section1ThetaBot1 = self.angle_conversion(section1AlphaBot1)
        section1ThetaBot2 =  np.pi + section1ThetaBot1

        # Angle between the XZ plane the MESH origin and vertex 7? 4? of the block 1
        section1PhiBot1 = self.angle_conversion(section1BetaBot1)
        section1PhiBot2 = np.pi/2 - section1PhiBot1

        # Capsule points top
        # Angle between teh XY plane, the sphere origin, and the center point
        # of the edge between vertex 1 and vertex 2.
        section1BetaTop1 = np.pi - 2 * np.arctan(2)
        section1BetaTop2 = np.pi - 2 * np.arctan(3)
        section1PhiTop1 = np.pi/2 - section1BetaTop1
        section1PhiTop2 = np.pi/2 - section1BetaTop2

        #print(f'block2beta1       [deg]: {np.rad2deg(block2BetaTop1):.4f}')
        #print(f'section1BetaTop1  [deg]: {np.rad2deg(section1BetaTop1):.4f}')
        #print(f'block2Beta2       [deg]: {np.rad2deg(block2BetaTop2):.4f}')
        #print(f'section1BetaTop2  [deg]: {np.rad2deg(section1BetaTop2):.4f}')
        #print(f'block2Phi1        [deg]: {np.rad2deg(block2PhiTop1):.4f}')
        #print(f'section1PhiTop1   [deg]: {np.rad2deg(section1PhiTop1):.4f}')
        #print(f'block2Phi2        [deg]: {np.rad2deg(block2PhiTop2):.4f}')
        #print(f'section1PhiTop2   [deg]: {np.rad2deg(section1PhiTop2):.4f}')

        # Outlet points top
        # TODO: Block section 3 and not block 3
        l = self.capsuleSphereCenter[0] / np.sqrt(2)
        m = self.capsuleRad - l
        section3PhiBot = np.arctan(self.capsuleSphereCenter[0] / m)
        # TODO: Optimise angle
        section4PhiBot1 = np.deg2rad(60)
        section4PhiBot2 = np.pi/2 - section4PhiBot1
        section5PhiBot1 = np.deg2rad(90)
        section5PhiBot2 = np.pi/2 - section5PhiBot1

        # Computation of 10 angles "central point" between the edges01 and edges 45
        N = 10
        anglesCapsule0 = np.linspace(np.pi/2, section1BetaBot2, N)
        # From "Top" to "Bottom"
        anglesCapsule1 = np.linspace(section1BetaTop1, section1BetaBot2, N)
        #print(f'anglesCapsule0  [deg]:\n {np.rad2deg(anglesCapsule0)}')
        #print(f'anglesCapsule1  [deg]:\n {np.rad2deg(anglesCapsule1)}')

        anglesInlet0  = np.linspace(section1PhiBot2, np.pi/2,         N)
        anglesInlet1  = np.linspace(section1PhiTop2, section1PhiBot2, N)
        anglesInlet2  = np.linspace(0,               section1PhiTop2, N)
        anglesOutlet3 = np.linspace(section3PhiBot, 0,                N)
        anglesOutlet4 = np.linspace(section4PhiBot1, section3PhiBot,  N)
        anglesOutlet5 = np.linspace(section5PhiBot1, section4PhiBot1, N)
        print(f'anglesInlet0  [deg]:\n {np.rad2deg(anglesInlet0)}')
        print(f'anglesInlet1  [deg]:\n {np.rad2deg(anglesInlet1)}')
        print(f'anglesInlet2  [deg]:\n {np.rad2deg(anglesInlet2)}')
        print(f'anglesOutlet3 [deg]:\n {np.rad2deg(anglesOutlet3)}')
        print(f'anglesOutlet4 [deg]:\n {np.rad2deg(anglesOutlet4)}')
        print(f'anglesOutlet5 [deg]:\n {np.rad2deg(anglesOutlet5)}')
        # Central points on wall/capsule

        # Front of the capsule
        # Section 0
        self.ccPnts0 = [self.point_on_sphere( self.capsuleRad, np.pi, ang, sphereOg) for ang in anglesCapsule0]
        # Section 1
        self.ccPnts1 = [self.point_on_sphere( self.capsuleRad, np.pi, ang, sphereOg) for ang in anglesCapsule1]
        # Section 2
        xFront = np.linspace(0, self.ccPnts1[0][0], N)
        self.ccPnts2 = [self.point_on_straight_front(x) for x in xFront]

        # Back of the capsule
        # Section 3
        xBack = np.linspace(sphereOg[0], 0, N)
        self.ccPnts3 = [self.point_on_straight_back(x) for x in xBack]

        # Section 4
        zTop = self.ccPnts3[0][-1]
        zBot = self.capsuleSphereCenter[0] * np.tan(section4PhiBot2)
        zBackS4 = np.linspace(zBot, zTop, N)
        self.ccPnts4 = [self.point_on_vertical_back(z) for z in zBackS4]
        print(f'zTop: {zTop}')
        print(f'zBot: {zBot}')
        print(f'zBackS4: {zBackS4}')
        print(f'{self.ccPnts4}')

        # Section 5
        zTop = self.capsuleSphereCenter[0] * np.tan(section4PhiBot2)
        zBot = 0
        zBackS5 = np.linspace(zBot, zTop, N)
        self.ccPnts5 = [self.point_on_vertical_back(z) for z in zBackS5]

        #print(3*'\n')
        #print([x[0] for x in self.ccPnts0])
        #print(self.ccPnts0[0][:])

        plotting = False
        plot1 = False
        if plotting or plot1:
            plt.figure()
            plt.plot([x[0] for x in self.ccPnts0], [x[2] for x in self.ccPnts0], 'o')
            plt.plot([x[0] for x in self.ccPnts1], [x[2] for x in self.ccPnts1], 'o')
            plt.plot([x[0] for x in self.ccPnts2], [x[2] for x in self.ccPnts2], 'o')
            plt.plot([x[0] for x in self.ccPnts3], [x[2] for x in self.ccPnts3], 'o')
            plt.plot([x[0] for x in self.ccPnts4], [x[2] for x in self.ccPnts4], 'o')
            plt.plot([x[0] for x in self.ccPnts5], [x[2] for x in self.ccPnts5], 'o')

        # Outside central points, on external mesh
        # Front points
        # Section 0
        self.iocPnts0 = [self.point_on_ellipsoid(np.pi, ang) for ang in anglesInlet0]
        self.DEBUGGING_ANGLES = anglesInlet0
        # Section 1
        self.iocPnts1 = [self.point_on_ellipsoid(np.pi, ang) for ang in anglesInlet1]
        # Section 2
        self.iocPnts2 = [self.point_on_ellipsoid(np.pi, ang) for ang in anglesInlet2]
        # Back of the capsule
        # TODO: compute the final phi angle for section 3
        # Section 3
        print(f'anglesOutlet3: {np.rad2deg(anglesOutlet3)}')
        self.iocPnts3 = [self.point_on_ellipsoid(0, ang) for ang in anglesOutlet3]
        self.iocPnts4 = [self.point_on_ellipsoid(0, ang) for ang in anglesOutlet4]
        self.iocPnts5 = [self.point_on_ellipsoid(0, ang) for ang in anglesOutlet5]

        if plotting or plot1:
            #plt.figure()
            plt.plot([x[0] for x in self.iocPnts0], [x[2] for x in self.iocPnts0], 'o')
            plt.plot([x[0] for x in self.iocPnts1], [x[2] for x in self.iocPnts1], 'o')
            plt.plot([x[0] for x in self.iocPnts2], [x[2] for x in self.iocPnts2], 'o')
            plt.plot([x[0] for x in self.iocPnts3], [x[2] for x in self.iocPnts3], 'o')
            plt.plot([x[0] for x in self.iocPnts4], [x[2] for x in self.iocPnts4], 'o')
            plt.plot([x[0] for x in self.iocPnts5], [x[2] for x in self.iocPnts5], 'o')

            plt.axis('equal')
            plt.show()

        plot2 = False
        if plotting or plot2:
            self.init_3d_plot()
            self.plot_points(self.ccPnts0, 'r')
            self.plot_points(self.ccPnts1, 'g')
            self.plot_points(self.ccPnts2, 'b')
            self.plot_points(self.ccPnts3, 'c')
            self.plot_points(self.ccPnts4, 'm')
            self.show_plot()

        """
        self.plot_points(ccPnts0,'r')
        self.plot_points(ccPnts1,'g')

        self.plot_points(icPnts0,'r')
        self.plot_points(icPnts1,'g')
        self.plot_points(icPnts2,'b')
        """

    def rotation_vector(self, e, rot):
        """
        Returns a rotation vector for quaternions.

        Parameters
        ----------
        eX: array
            Vector around which the rotation happends
        rot: float
            Rotation angle

        Returns
        -------
        vector: array
            quaternion definition vector
        """
        vector = [
            np.cos(rot/2),
            e[0]*np.sin(rot/2),
            e[1]*np.sin(rot/2),
            e[2]*np.sin(rot/2),
        ]
        return vector

    def central_points_rotation(self):
        """
        Rotates the points into each quadrant
        """
        # Rotation construction
        quadrant1Angle = np.deg2rad( 45)
        quadrant1cAngle= np.deg2rad(  0)
        quadrant2Angle = np.deg2rad(135)
        quadrant2cAngle= np.deg2rad( 90)
        quadrant3Angle = np.deg2rad(225)
        quadrant3cAngle= np.deg2rad(180)
        quadrant4Angle = np.deg2rad(315)
        quadrant4cAngle= np.deg2rad(270)

        # Unit x vector
        eX = [0, 0, 1]

        # Rotation vector for quaternions
        quaternionRotation1  = self.rotation_vector(eX, quadrant1Angle)
        quaternionRotationC1 = self.rotation_vector(eX, quadrant1cAngle)
        quaternionRotation2  = self.rotation_vector(eX, quadrant2Angle)
        quaternionRotationC2 = self.rotation_vector(eX, quadrant2cAngle)
        quaternionRotation3  = self.rotation_vector(eX, quadrant3Angle)
        quaternionRotationC3 = self.rotation_vector(eX, quadrant3cAngle)
        quaternionRotation4  = self.rotation_vector(eX, quadrant4Angle)
        quaternionRotationC4 = self.rotation_vector(eX, quadrant4cAngle)

        # Quaternion self.rotation functions defined
        self.rotQuat1  = R.from_quat(quaternionRotation1)
        self.rotQuatC1 = R.from_quat(quaternionRotationC1)
        self.rotQuat2  = R.from_quat(quaternionRotation2)
        self.rotQuatC2 = R.from_quat(quaternionRotationC2)
        self.rotQuat3  = R.from_quat(quaternionRotation3)
        self.rotQuatC3 = R.from_quat(quaternionRotationC3)
        self.rotQuat4  = R.from_quat(quaternionRotation4)
        self.rotQuatC4 = R.from_quat(quaternionRotationC4)

        # Saves the self.rotated capsule construction points
        # Section 0
        self.capsuleSection0quad1 = self.rotQuat1.apply(self.ccPnts0)
        self.capsuleSection0quad2 = self.rotQuat2.apply(self.ccPnts0)
        self.capsuleSection0quad3 = self.rotQuat3.apply(self.ccPnts0)
        self.capsuleSection0quad4 = self.rotQuat4.apply(self.ccPnts0)
        # Section 1
        self.capsuleSection1quad1  = self.rotQuat1.apply(self.ccPnts1)
        self.capsuleSection1quadC1 = self.rotQuatC1.apply(self.ccPnts1)
        self.capsuleSection1quad2  = self.rotQuat2.apply(self.ccPnts1)
        self.capsuleSection1quadC2 = self.rotQuatC2.apply(self.ccPnts1)
        self.capsuleSection1quad3  = self.rotQuat3.apply(self.ccPnts1)
        self.capsuleSection1quadC3 = self.rotQuatC3.apply(self.ccPnts1)
        self.capsuleSection1quad4  = self.rotQuat4.apply(self.ccPnts1)
        self.capsuleSection1quadC4 = self.rotQuatC4.apply(self.ccPnts1)
        # Section 2
        self.capsuleSection2quad1  = self.rotQuat1.apply(self.ccPnts2)
        self.capsuleSection2quadC1 = self.rotQuatC1.apply(self.ccPnts2)
        self.capsuleSection2quad2  = self.rotQuat2.apply(self.ccPnts2)
        self.capsuleSection2quadC2 = self.rotQuatC2.apply(self.ccPnts2)
        self.capsuleSection2quad3  = self.rotQuat3.apply(self.ccPnts2)
        self.capsuleSection2quadC3 = self.rotQuatC3.apply(self.ccPnts2)
        self.capsuleSection2quad4  = self.rotQuat4.apply(self.ccPnts2)
        self.capsuleSection2quadC4 = self.rotQuatC4.apply(self.ccPnts2)
        # Section 3
        self.capsuleSection3quad1  = self.rotQuat1.apply(self.ccPnts3)
        self.capsuleSection3quadC1 = self.rotQuatC1.apply(self.ccPnts3)
        self.capsuleSection3quad2  = self.rotQuat2.apply(self.ccPnts3)
        self.capsuleSection3quadC2 = self.rotQuatC2.apply(self.ccPnts3)
        self.capsuleSection3quad3  = self.rotQuat3.apply(self.ccPnts3)
        self.capsuleSection3quadC3 = self.rotQuatC3.apply(self.ccPnts3)
        self.capsuleSection3quad4  = self.rotQuat4.apply(self.ccPnts3)
        self.capsuleSection3quadC4 = self.rotQuatC4.apply(self.ccPnts3)
        # Section 4
        self.capsuleSection4quad1  = self.rotQuat1.apply(self.ccPnts4)
        self.capsuleSection4quadC1 = self.rotQuatC1.apply(self.ccPnts4)
        self.capsuleSection4quad2  = self.rotQuat2.apply(self.ccPnts4)
        self.capsuleSection4quadC2 = self.rotQuatC2.apply(self.ccPnts4)
        self.capsuleSection4quad3  = self.rotQuat3.apply(self.ccPnts4)
        self.capsuleSection4quadC3 = self.rotQuatC3.apply(self.ccPnts4)
        self.capsuleSection4quad4  = self.rotQuat4.apply(self.ccPnts4)
        self.capsuleSection4quadC4 = self.rotQuatC4.apply(self.ccPnts4)
        # Section 5
        self.capsuleSection5quad1  = self.rotQuat1.apply(self.ccPnts5)
        self.capsuleSection5quadC1 = self.rotQuatC1.apply(self.ccPnts5)
        self.capsuleSection5quad2  = self.rotQuat2.apply(self.ccPnts5)
        self.capsuleSection5quadC2 = self.rotQuatC2.apply(self.ccPnts5)
        self.capsuleSection5quad3  = self.rotQuat3.apply(self.ccPnts5)
        self.capsuleSection5quadC3 = self.rotQuatC3.apply(self.ccPnts5)
        self.capsuleSection5quad4  = self.rotQuat4.apply(self.ccPnts5)
        self.capsuleSection5quadC4 = self.rotQuatC4.apply(self.ccPnts5)

        # Saves the rotated inlet construction points
        # Section 0
        self.ioSection0quad1 = self.rotQuat1.apply(self.iocPnts0)
        self.ioSection0quad2 = self.rotQuat2.apply(self.iocPnts0)
        self.ioSection0quad3 = self.rotQuat3.apply(self.iocPnts0)
        self.ioSection0quad4 = self.rotQuat4.apply(self.iocPnts0)
        # Section 1
        self.ioSection1quad1 = self.rotQuat1.apply(self.iocPnts1)
        self.ioSection1quadC1= self.rotQuatC1.apply(self.iocPnts1)
        self.ioSection1quad2 = self.rotQuat2.apply(self.iocPnts1)
        self.ioSection1quadC2= self.rotQuatC2.apply(self.iocPnts1)
        self.ioSection1quad3 = self.rotQuat3.apply(self.iocPnts1)
        self.ioSection1quadC3= self.rotQuatC3.apply(self.iocPnts1)
        self.ioSection1quad4 = self.rotQuat4.apply(self.iocPnts1)
        self.ioSection1quadC4= self.rotQuatC4.apply(self.iocPnts1)
        # Section 2
        self.ioSection2quad1 = self.rotQuat1.apply(self.iocPnts2)
        self.ioSection2quadC1= self.rotQuatC1.apply(self.iocPnts2)
        self.ioSection2quad2 = self.rotQuat2.apply(self.iocPnts2)
        self.ioSection2quadC2= self.rotQuatC2.apply(self.iocPnts2)
        self.ioSection2quad3 = self.rotQuat3.apply(self.iocPnts2)
        self.ioSection2quadC3= self.rotQuatC3.apply(self.iocPnts2)
        self.ioSection2quad4 = self.rotQuat4.apply(self.iocPnts2)
        self.ioSection2quadC4= self.rotQuatC4.apply(self.iocPnts2)
        # Section 3
        self.ioSection3quad1 = self.rotQuat1.apply(self.iocPnts3)
        self.ioSection3quadC1= self.rotQuatC1.apply(self.iocPnts3)
        self.ioSection3quad2 = self.rotQuat2.apply(self.iocPnts3)
        self.ioSection3quadC2= self.rotQuatC2.apply(self.iocPnts3)
        self.ioSection3quad3 = self.rotQuat3.apply(self.iocPnts3)
        self.ioSection3quadC3= self.rotQuatC3.apply(self.iocPnts3)
        self.ioSection3quad4 = self.rotQuat4.apply(self.iocPnts3)
        self.ioSection3quadC4= self.rotQuatC4.apply(self.iocPnts3)
        # Section 4
        self.ioSection4quad1 = self.rotQuat1.apply(self.iocPnts4)
        self.ioSection4quadC1= self.rotQuatC1.apply(self.iocPnts4)
        self.ioSection4quad2 = self.rotQuat2.apply(self.iocPnts4)
        self.ioSection4quadC2= self.rotQuatC2.apply(self.iocPnts4)
        self.ioSection4quad3 = self.rotQuat3.apply(self.iocPnts4)
        self.ioSection4quadC3= self.rotQuatC3.apply(self.iocPnts4)
        self.ioSection4quad4 = self.rotQuat4.apply(self.iocPnts4)
        self.ioSection4quadC4= self.rotQuatC4.apply(self.iocPnts4)
        # Section 5
        self.ioSection5quad1 = self.rotQuat1.apply(self.iocPnts5)
        self.ioSection5quadC1= self.rotQuatC1.apply(self.iocPnts5)
        self.ioSection5quad2 = self.rotQuat2.apply(self.iocPnts5)
        self.ioSection5quadC2= self.rotQuatC2.apply(self.iocPnts5)
        self.ioSection5quad3 = self.rotQuat3.apply(self.iocPnts5)
        self.ioSection5quadC3= self.rotQuatC3.apply(self.iocPnts5)
        self.ioSection5quad4 = self.rotQuat4.apply(self.iocPnts5)
        self.ioSection5quadC4= self.rotQuatC4.apply(self.iocPnts5)

    # TODO: Correct all the spline edges that are an arc edge
    def build_block1(self):
        """
        Computes vertices and point in order to build block 1.
        """
        # Computes vertices for block 1
        blk1v0 = self.capsuleSection0quad1[-1]
        blk1v1 = self.capsuleSection0quad2[-1]
        blk1v2 = self.capsuleSection0quad3[-1]
        blk1v3 = self.capsuleSection0quad4[-1]
        blk1v4 = self.ioSection0quad1[0]
        blk1v5 = self.ioSection0quad2[0]
        blk1v6 = self.ioSection0quad3[0]
        blk1v7 = self.ioSection0quad4[0]

        self.block1Points = [
            blk1v0, blk1v1, blk1v2, blk1v3,
            blk1v4, blk1v5, blk1v6, blk1v7,
        ]

        # Computes Alpha and Beta for the capsule surface
        # First we need the vector
        vCapsuleFront   = blk1v1 - self.capsuleSphereCenter
        vCapsuleFrontXY = [vCapsuleFront[0], vCapsuleFront[1], 0]
        # Since eX = [1, 0, 0] the angle can be easily found in the following way
        alphaCapsule = np.arccos(vCapsuleFrontXY[0] / np.linalg.norm(vCapsuleFrontXY))
        # Since eZ = [0, 0, 1] the angle can be easily found in the following way
        betaCapsule   = np.arccos(vCapsuleFront[2] / np.linalg.norm(vCapsuleFront))

        # Computes Theta and Phi for the inlet surface
        vInletFront   = blk1v6
        vInletFrontXY = [vInletFront[0], vInletFront[1], 0]
        # WARNING: works only in the first quadrand!
        # Since eX = [1, 0, 0] the angle can be easily found in the following way
        theta = np.arccos(vInletFront[2]   / self.c2)
        thetaInlet = np.arccos(vInletFront[0]/(self.a2*np.sin(theta)))
        # Since eZ = [0, 0, 1] the angle can be easily found in the following way
        phiInlet   = np.arccos(vInletFront[2]   / self.c2)

        # TODO: Implement a test for this

        """
        print(f'debugging angles:{np.rad2deg(self.DEBUGGING_ANGLES)}')
        print(f'thetaInlet: {np.rad2deg(thetaInlet)}')
        print(f'phiInlet:   {np.rad2deg(phiInlet)}')
        """

        """
        # Test to see if the angles found in spherical coordiantes are correct
        p = self.point_on_sphere(
            rho   = self.capsuleRad,
            theta = alphaCapsule,
            phi   = betaCapsule,
            origin= self.capsuleSphereCenter)
        print(f'p:             {p}')
        print(f'vCapsuleFront: {blk1v2}')
        """

        """
        p = self.point_on_ellipsoid(thetaInlet, phiInlet)
        print(f'p:          {p}')
        print(f'vInletBack: {blk1v6}')
        #pp = self.point_on_ellipsoid(np.pi, phiInlet)
        #print(f'pp:         {pp}')
        for pnt in self.iocPnts0:
            print(pnt)
        """

        N = 10
        theta = np.pi - thetaInlet
        thetaMin = np.pi + theta
        thetaMax = np.pi - theta
        angles = np.linspace(thetaMin, thetaMax, N)

        # Arc edges on capsule surface
        block1CapsuleCenterPoint = self.point_on_sphere(
            rho   = self.capsuleRad,
            theta = np.pi,
            phi   = betaCapsule,
            origin= self.capsuleSphereCenter)
        # Rotates the point. This point will be used as center point of the
        # edges in order to build a arc edge on the capsule surface.
        c1 = self.rotQuatC1.apply(block1CapsuleCenterPoint)
        c2 = self.rotQuatC2.apply(block1CapsuleCenterPoint)
        c3 = self.rotQuatC3.apply(block1CapsuleCenterPoint)
        c4 = self.rotQuatC4.apply(block1CapsuleCenterPoint)

        # Spline edges for the inlet surface
        c5 = []  # Spline points for ege 56
        c6 = []  # Spline points for ege 67
        c7 = []  # Spline points for ege 74
        c8 = []  # Spline points for ege 45

        for theta in angles:
            point = self.point_on_ellipsoid(theta, phiInlet)

            tmp = self.rotQuatC1.apply(point)
            c5.append(tmp)
            tmp = self.rotQuatC2.apply(point)
            c6.append(tmp)
            tmp = self.rotQuatC3.apply(point)
            c7.append(tmp)
            tmp = self.rotQuatC4.apply(point)
            c8.append(tmp)

        self.block1Edges = [
            Edge(0, 1, c2),  # Arc Edge
            Edge(1, 2, c3),  # Arc Edge
            Edge(2, 3, c4),  # Arc Edge
            Edge(3, 0, c1),  # Arc Edge
            Edge(7, 4, c5),  # Spline Edge
            Edge(4, 5, c6),  # Spline Edge
            Edge(5, 6, c7),  # Spline Edge
            Edge(6, 7, c8),  # Spline Edge
        ]

        self.block1 = Block.create_from_points(
            self.block1Points,
            self.block1Edges
        )
        self.block1.set_patch('top', 'inlet')
        self.block1.set_patch('bottom', 'wall')

        self.block1.chop(0, count=10, c2c_expansion=1)
        self.block1.chop(1, count=10, c2c_expansion=1)
        self.block1.chop(2, count=10, c2c_expansion=1)

    def build_block2(self):
        """
        Computes block 2 vertices and points
        """
        # Computes vertices for block 2
        blk2v0 = self.capsuleSection1quad1[-1]
        blk2v1 = self.capsuleSection1quad4[-1]
        blk2v2 = self.capsuleSection1quad4[0]
        blk2v3 = self.capsuleSection1quad1[0]
        blk2v4 = self.ioSection1quad1[-1]
        blk2v5 = self.ioSection1quad4[-1]
        blk2v6 = self.ioSection1quad4[0]
        blk2v7 = self.ioSection1quad1[0]

        self.block2Points = [
            blk2v0, blk2v1, blk2v2, blk2v3,
            blk2v4, blk2v5, blk2v6, blk2v7,
        ]

        self.block2Edges = [
            Edge(1, 2, self.capsuleSection1quad4[4]),   # Arc Edge
            Edge(2, 3, self.capsuleSection1quadC1[0]),  # Arc Edge
            Edge(3, 0, self.capsuleSection1quad1[4]),   # Arc Edge
            Edge(5, 6, self.ioSection1quad4[::-1]),     # Spline Edge
            Edge(6, 7, self.ioSection1quadC1[0]),       # Arc Edge
            Edge(7, 4, self.ioSection1quad1),           # Spline Edge
        ]

        self.block2 = Block.create_from_points(
            self.block2Points,
            self.block2Edges
        )
        self.block2.set_patch('top', 'inlet')
        self.block2.set_patch('bottom', 'wall')

        self.block2.chop(0, count=10, c2c_expansion=1)
        self.block2.chop(1, count=10, c2c_expansion=1)
        self.block2.chop(2, count=10, c2c_expansion=1)

    def build_block3(self):
        """
        Computes block 3 vertices and points
        """
        # Computes vertices for block 2
        blk3v0 = self.capsuleSection1quad2[-1]
        blk3v1 = self.capsuleSection1quad1[-1]
        blk3v2 = self.capsuleSection1quad1[0]
        blk3v3 = self.capsuleSection1quad2[0]
        blk3v4 = self.ioSection1quad2[-1]
        blk3v5 = self.ioSection1quad1[-1]
        blk3v6 = self.ioSection1quad1[0]
        blk3v7 = self.ioSection1quad2[0]

        self.block3Points = [
            blk3v0, blk3v1, blk3v2, blk3v3,
            blk3v4, blk3v5, blk3v6, blk3v7,
        ]

        self.block3Edges = [
            # The missing edges are defined in block 1 and 2
            Edge(2, 3, self.capsuleSection1quadC2[0]),  # Arc Edge
            Edge(3, 0, self.capsuleSection1quad2[4]),   # Arc Edge
            Edge(6, 7, self.ioSection1quadC2[0]),       # Arc Edge
            Edge(7, 4, self.ioSection1quad2),           # Spline Edge
        ]

        self.block3 = Block.create_from_points(
            self.block3Points,
            self.block3Edges
        )
        self.block3.set_patch('top', 'inlet')
        self.block3.set_patch('bottom', 'wall')

        self.block3.chop(0, count=10, c2c_expansion=1)
        self.block3.chop(1, count=10, c2c_expansion=1)
        self.block3.chop(2, count=10, c2c_expansion=1)

    def build_block4(self):
        """
        Computes block 4 vertices and points
        """
        # Computes vertices for block 4
        blk4v0 = self.capsuleSection1quad3[-1]
        blk4v1 = self.capsuleSection1quad2[-1]
        blk4v2 = self.capsuleSection1quad2[0]
        blk4v3 = self.capsuleSection1quad3[0]
        blk4v4 = self.ioSection1quad3[-1]
        blk4v5 = self.ioSection1quad2[-1]
        blk4v6 = self.ioSection1quad2[0]
        blk4v7 = self.ioSection1quad3[0]

        self.block4Points = [
            blk4v0, blk4v1, blk4v2, blk4v3,
            blk4v4, blk4v5, blk4v6, blk4v7,
        ]

        self.block4Edges = [
            # The missing edges are defined in block 1, 2, 3
            Edge(2, 3, self.capsuleSection1quadC3[0]),  # Arc Edge
            Edge(3, 0, self.capsuleSection1quad3[4]),   # Arc Edge
            Edge(6, 7, self.ioSection1quadC3[0]),       # Arc Edge
            Edge(7, 4, self.ioSection1quad3),           # Spline Edge
        ]

        self.block4 = Block.create_from_points(
            self.block4Points,
            self.block4Edges
        )
        self.block4.set_patch('top', 'inlet')
        self.block4.set_patch('bottom', 'wall')

        self.block4.chop(0, count=10, c2c_expansion=1)
        self.block4.chop(1, count=10, c2c_expansion=1)
        self.block4.chop(2, count=10, c2c_expansion=1)

    def build_block5(self):
        """
        Computes block 5 vertices and points
        """
        # Computes vertices for block 2
        blk5v0 = self.capsuleSection1quad4[-1]
        blk5v1 = self.capsuleSection1quad3[-1]
        blk5v2 = self.capsuleSection1quad3[0]
        blk5v3 = self.capsuleSection1quad4[0]
        blk5v4 = self.ioSection1quad4[-1]
        blk5v5 = self.ioSection1quad3[-1]
        blk5v6 = self.ioSection1quad3[0]
        blk5v7 = self.ioSection1quad4[0]

        self.block5Points = [
            blk5v0, blk5v1, blk5v2, blk5v3,
            blk5v4, blk5v5, blk5v6, blk5v7,
        ]

        self.block5Edges = [
            # The missing edges are defined in previous blocks
            Edge(2, 3, self.capsuleSection1quadC4[0]),  # Arc Edge
            Edge(6, 7, self.ioSection1quadC4[0]),       # Arc Edge
        ]

        self.block5 = Block.create_from_points(
            self.block5Points,
            self.block5Edges
        )
        self.block5.set_patch('top', 'inlet')
        self.block5.set_patch('bottom', 'wall')

        self.block5.chop(0, count=10, c2c_expansion=1)
        self.block5.chop(1, count=10, c2c_expansion=1)
        self.block5.chop(2, count=10, c2c_expansion=1)

    def build_block6(self):
        """
        Computes block 6 vertices and points
        """
        # Computes vertices for block 6
        blk6v0 = self.capsuleSection2quad1[-1]
        blk6v1 = self.capsuleSection2quad4[-1]
        blk6v2 = self.capsuleSection2quad4[0]
        blk6v3 = self.capsuleSection2quad1[0]
        blk6v4 = self.ioSection2quad1[-1]
        blk6v5 = self.ioSection2quad4[-1]
        blk6v6 = self.ioSection2quad4[0]
        blk6v7 = self.ioSection2quad1[0]

        self.block6Points = [
            blk6v0, blk6v1, blk6v2, blk6v3,
            blk6v4, blk6v5, blk6v6, blk6v7,
        ]

        self.block6Edges = [
            # The missing edges are defined in previous blocks
            Edge(2, 3, self.capsuleSection2quadC1[0]),  # Arc Edge
            Edge(5, 6, self.ioSection2quad4[::-1]),     # Spline Edge
            Edge(6, 7, self.ioSection2quadC1[0]),       # Arc Edge
            Edge(7, 4, self.ioSection2quad1),           # Spline Edge
        ]

        self.block6 = Block.create_from_points(
            self.block6Points,
            self.block6Edges
        )
        self.block6.set_patch('top', 'inlet')
        self.block6.set_patch('bottom', 'wall')

        self.block6.chop(0, count=10, c2c_expansion=1)
        self.block6.chop(1, count=10, c2c_expansion=1)
        self.block6.chop(2, count=10, c2c_expansion=1)

    def build_block7(self):
        """
        Computes block 7 vertices and points
        """
        # Computes vertices for block 7
        blk7v0 = self.capsuleSection2quad2[-1]
        blk7v1 = self.capsuleSection2quad1[-1]
        blk7v2 = self.capsuleSection2quad1[0]
        blk7v3 = self.capsuleSection2quad2[0]
        blk7v4 = self.ioSection2quad2[-1]
        blk7v5 = self.ioSection2quad1[-1]
        blk7v6 = self.ioSection2quad1[0]
        blk7v7 = self.ioSection2quad2[0]

        self.block7Points = [
            blk7v0, blk7v1, blk7v2, blk7v3,
            blk7v4, blk7v5, blk7v6, blk7v7,
        ]

        self.block7Edges = [
            # The missing edges are defined in previous blocks
            Edge(2, 3, self.capsuleSection2quadC2[0]),  # Arc Edge
            Edge(6, 7, self.ioSection2quadC2[0]),       # Arc Edge
            Edge(7, 4, self.ioSection2quad2),           # Spline Edge
        ]

        self.block7 = Block.create_from_points(
            self.block7Points,
            self.block7Edges
        )
        self.block7.set_patch('top', 'inlet')
        self.block7.set_patch('bottom', 'wall')

        self.block7.chop(0, count=10, c2c_expansion=1)
        self.block7.chop(1, count=10, c2c_expansion=1)
        self.block7.chop(2, count=10, c2c_expansion=1)

    def build_block8(self):
        """
        Computes block 8 vertices and points
        """
        # Computes vertices for block 8
        blk8v0 = self.capsuleSection2quad3[-1]
        blk8v1 = self.capsuleSection2quad2[-1]
        blk8v2 = self.capsuleSection2quad2[0]
        blk8v3 = self.capsuleSection2quad3[0]
        blk8v4 = self.ioSection2quad3[-1]
        blk8v5 = self.ioSection2quad2[-1]
        blk8v6 = self.ioSection2quad2[0]
        blk8v7 = self.ioSection2quad3[0]

        self.block8Points = [
            blk8v0, blk8v1, blk8v2, blk8v3,
            blk8v4, blk8v5, blk8v6, blk8v7,
        ]

        self.block8Edges = [
            # The missing edges are defined in previous blocks
            Edge(2, 3, self.capsuleSection2quadC3[0]),  # Arc Edge
            Edge(6, 7, self.ioSection2quadC3[0]),       # Arc Edge
            Edge(7, 4, self.ioSection2quad3),           # Spline Edge
        ]

        self.block8 = Block.create_from_points(
            self.block8Points,
            self.block8Edges
        )
        self.block8.set_patch('top', 'inlet')
        self.block8.set_patch('bottom', 'wall')

        self.block8.chop(0, count=10, c2c_expansion=1)
        self.block8.chop(1, count=10, c2c_expansion=1)
        self.block8.chop(2, count=10, c2c_expansion=1)

    def build_block9(self):
        """
        Computes block 9 vertices and points
        """
        # Computes vertices for block 9
        blk9v0 = self.capsuleSection2quad4[-1]
        blk9v1 = self.capsuleSection2quad3[-1]
        blk9v2 = self.capsuleSection2quad3[0]
        blk9v3 = self.capsuleSection2quad4[0]
        blk9v4 = self.ioSection2quad4[-1]
        blk9v5 = self.ioSection2quad3[-1]
        blk9v6 = self.ioSection2quad3[0]
        blk9v7 = self.ioSection2quad4[0]

        self.block9Points = [
            blk9v0, blk9v1, blk9v2, blk9v3,
            blk9v4, blk9v5, blk9v6, blk9v7,
        ]

        self.block9Edges = [
            # The missing edges are defined in previous blocks
            Edge(2, 3, self.capsuleSection2quadC4[0]),  # Arc Edge
            Edge(6, 7, self.ioSection2quadC4[0]),       # Arc Edge
        ]

        self.block9 = Block.create_from_points(
            self.block9Points,
            self.block9Edges
        )
        self.block9.set_patch('top', 'inlet')
        self.block9.set_patch('bottom', 'wall')

        self.block9.chop(0, count=10, c2c_expansion=1)
        self.block9.chop(1, count=10, c2c_expansion=1)
        self.block9.chop(2, count=10, c2c_expansion=1)

    def build_block10(self):
        """
        Computes block 10 vertices and points
        """
        # Computes vertices for block 10
        blk10v0 = self.capsuleSection3quad1[-1]
        blk10v1 = self.capsuleSection3quad4[-1]
        blk10v2 = self.capsuleSection3quad4[0]
        blk10v3 = self.capsuleSection3quad1[0]
        blk10v4 = self.ioSection3quad1[-1]
        blk10v5 = self.ioSection3quad4[-1]
        blk10v6 = self.ioSection3quad4[0]
        blk10v7 = self.ioSection3quad1[0]

        self.block10Points = [
            blk10v0, blk10v1, blk10v2, blk10v3,
            blk10v4, blk10v5, blk10v6, blk10v7,
        ]
        print('\n block10')
        for pnt in self.block10Points:
            print(pnt)
        self.block10Edges = [
            # The missing edges are defined in previous blocks
            Edge(2, 3, self.capsuleSection3quadC1[0]),  # Arc Edge
            Edge(5, 6, self.ioSection3quad4[::-1]),     # Spline Edge
            Edge(6, 7, self.ioSection3quadC1[0]),       # Arc Edge
            Edge(7, 4, self.ioSection3quad1),           # Spline Edge
        ]

        self.block10 = Block.create_from_points(
            self.block10Points,
            self.block10Edges
        )
        self.block10.set_patch('top', 'outlet')
        self.block10.set_patch('bottom', 'wall')

        self.block10.chop(0, count=10, c2c_expansion=1)
        self.block10.chop(1, count=10, c2c_expansion=1)
        self.block10.chop(2, count=10, c2c_expansion=1)

    def build_block11(self):
        """
        Computes block 11 vertices and points
        """
        # Computes vertices for block 11
        blk11v0 = self.capsuleSection3quad2[-1]
        blk11v1 = self.capsuleSection3quad1[-1]
        blk11v2 = self.capsuleSection3quad1[0]
        blk11v3 = self.capsuleSection3quad2[0]
        blk11v4 = self.ioSection3quad2[-1]
        blk11v5 = self.ioSection3quad1[-1]
        blk11v6 = self.ioSection3quad1[0]
        blk11v7 = self.ioSection3quad2[0]

        self.block11Points = [
            blk11v0, blk11v1, blk11v2, blk11v3,
            blk11v4, blk11v5, blk11v6, blk11v7,
        ]

        self.block11Edges = [
            # The missing edges are defined in previous blocks
            Edge(2, 3, self.capsuleSection3quadC2[0]),  # Arc Edge
            Edge(6, 7, self.ioSection3quadC2[0]),       # Arc Edge
            Edge(7, 4, self.ioSection3quad2),           # Spline Edge
        ]

        self.block11 = Block.create_from_points(
            self.block11Points,
            self.block11Edges
        )
        self.block11.set_patch('top', 'outlet')
        self.block11.set_patch('bottom', 'wall')

        self.block11.chop(0, count=10, c2c_expansion=1)
        self.block11.chop(1, count=10, c2c_expansion=1)
        self.block11.chop(2, count=10, c2c_expansion=1)

    def build_block12(self):
        """
        Computes block 12 vertices and points
        """
        # Computes vertices for block 12
        blk12v0 = self.capsuleSection3quad3[-1]
        blk12v1 = self.capsuleSection3quad2[-1]
        blk12v2 = self.capsuleSection3quad2[0]
        blk12v3 = self.capsuleSection3quad3[0]
        blk12v4 = self.ioSection3quad3[-1]
        blk12v5 = self.ioSection3quad2[-1]
        blk12v6 = self.ioSection3quad2[0]
        blk12v7 = self.ioSection3quad3[0]

        self.block12Points = [
            blk12v0, blk12v1, blk12v2, blk12v3,
            blk12v4, blk12v5, blk12v6, blk12v7,
        ]

        self.block12Edges = [
            # The missing edges are defined in previous blocks
            Edge(2, 3, self.capsuleSection3quadC3[0]),  # Arc Edge
            Edge(6, 7, self.ioSection3quadC3[0]),       # Arc Edge
            Edge(7, 4, self.ioSection3quad3),           # Spline Edge
        ]

        self.block12 = Block.create_from_points(
            self.block12Points,
            self.block12Edges
        )
        self.block12.set_patch('top', 'outlet')
        self.block12.set_patch('bottom', 'wall')

        self.block12.chop(0, count=10, c2c_expansion=1)
        self.block12.chop(1, count=10, c2c_expansion=1)
        self.block12.chop(2, count=10, c2c_expansion=1)

    def build_block13(self):
        """
        Computes block 13 vertices and points
        """
        # Computes vertices for block 13
        blk13v0 = self.capsuleSection3quad4[-1]
        blk13v1 = self.capsuleSection3quad3[-1]
        blk13v2 = self.capsuleSection3quad3[0]
        blk13v3 = self.capsuleSection3quad4[0]
        blk13v4 = self.ioSection3quad4[-1]
        blk13v5 = self.ioSection3quad3[-1]
        blk13v6 = self.ioSection3quad3[0]
        blk13v7 = self.ioSection3quad4[0]

        self.block13Points = [
            blk13v0, blk13v1, blk13v2, blk13v3,
            blk13v4, blk13v5, blk13v6, blk13v7,
        ]

        self.block13Edges = [
            # The missing edges are defined in previous blocks
            Edge(2, 3, self.capsuleSection3quadC4[0]),  # Arc Edge
            Edge(6, 7, self.ioSection3quadC4[0]),       # Arc Edge
        ]

        self.block13 = Block.create_from_points(
            self.block13Points,
            self.block13Edges)
        self.block13.set_patch('top', 'outlet')
        self.block13.set_patch('bottom', 'wall')

        self.block13.chop(0, count=10, c2c_expansion=1)
        self.block13.chop(1, count=10, c2c_expansion=1)
        self.block13.chop(2, count=10, c2c_expansion=1)

    def build_block14(self):
        """
        Computes block 14 vertices and points
        """
        # Computes vertices for block 14
        blk14v0 = self.capsuleSection4quad1[-1]
        blk14v1 = self.capsuleSection4quad4[-1]
        blk14v2 = self.capsuleSection4quad4[0]
        blk14v3 = self.capsuleSection4quad1[0]
        blk14v4 = self.ioSection4quad1[-1]
        blk14v5 = self.ioSection4quad4[-1]
        blk14v6 = self.ioSection4quad4[0]
        blk14v7 = self.ioSection4quad1[0]

        self.block14Points = [
            blk14v0, blk14v1, blk14v2, blk14v3,
            blk14v4, blk14v5, blk14v6, blk14v7,
        ]

        self.block14Edges = [
            # The missing edges are defined in previous blocks
            Edge(5, 6, self.ioSection4quad4[::-1]),    # Spline Edge
            Edge(7, 4, self.ioSection4quad1),          # Spline Edge
        ]

        self.block14 = Block.create_from_points(
            self.block14Points,
            self.block14Edges
        )
        self.block14.set_patch('top', 'outlet')
        self.block14.set_patch('bottom', 'wall')

        self.block14.chop(0, count=10, c2c_expansion=1)
        self.block14.chop(1, count=10, c2c_expansion=1)
        self.block14.chop(2, count=10, c2c_expansion=1)

    def build_block15(self):
        """
        Computes block 15 vertices and points
        """
        # Computes vertices for block 15
        blk15v0 = self.capsuleSection4quad2[-1]
        blk15v1 = self.capsuleSection4quad1[-1]
        blk15v2 = self.capsuleSection4quad1[0]
        blk15v3 = self.capsuleSection4quad2[0]
        blk15v4 = self.ioSection4quad2[-1]
        blk15v5 = self.ioSection4quad1[-1]
        blk15v6 = self.ioSection4quad1[0]
        blk15v7 = self.ioSection4quad2[0]

        self.block15Points = [
            blk15v0, blk15v1, blk15v2, blk15v3,
            blk15v4, blk15v5, blk15v6, blk15v7,
        ]

        self.block15Edges = [
            # The missing edges are defined in previous blocks
            Edge(7, 4, self.ioSection4quad2),  # Spline Edge
        ]

        self.block15 = Block.create_from_points(
            self.block15Points,
            self.block15Edges
        )
        self.block15.set_patch('top', 'outlet')
        self.block15.set_patch('bottom', 'wall')

        self.block15.chop(0, count=10, c2c_expansion=1)
        self.block15.chop(1, count=10, c2c_expansion=1)
        self.block15.chop(2, count=10, c2c_expansion=1)

    def build_block16(self):
        """
        Computes block 16 vertices and points
        """
        # Computes vertices for block 16
        blk16v0 = self.capsuleSection4quad3[-1]
        blk16v1 = self.capsuleSection4quad2[-1]
        blk16v2 = self.capsuleSection4quad2[0]
        blk16v3 = self.capsuleSection4quad3[0]
        blk16v4 = self.ioSection4quad3[-1]
        blk16v5 = self.ioSection4quad2[-1]
        blk16v6 = self.ioSection4quad2[0]
        blk16v7 = self.ioSection4quad3[0]

        self.block16Points = [
            blk16v0, blk16v1, blk16v2, blk16v3,
            blk16v4, blk16v5, blk16v6, blk16v7,
        ]

        self.block16Edges = [
            # The missing edges are defined in previous blocks
            Edge(7, 4, self.ioSection4quad3),          # Spline Edge
        ]

        self.block16 = Block.create_from_points(
            self.block16Points,
            self.block16Edges
        )
        self.block16.set_patch('top', 'outlet')
        self.block16.set_patch('bottom', 'wall')

        self.block16.chop(0, count=10, c2c_expansion=1)
        self.block16.chop(1, count=10, c2c_expansion=1)
        self. block16.chop(2, count=10, c2c_expansion=1)

    def build_block17(self):
        """
        Computes block 13 vertices and points
        """
        # Computes vertices for block 13
        blk17v0 = self.capsuleSection4quad4[-1]
        blk17v1 = self.capsuleSection4quad3[-1]
        blk17v2 = self.capsuleSection4quad3[0]
        blk17v3 = self.capsuleSection4quad4[0]
        blk17v4 = self.ioSection4quad4[-1]
        blk17v5 = self.ioSection4quad3[-1]
        blk17v6 = self.ioSection4quad3[0]
        blk17v7 = self.ioSection4quad4[0]

        self.block17Points = [
            blk17v0, blk17v1, blk17v2, blk17v3,
            blk17v4, blk17v5, blk17v6, blk17v7,
        ]

        # Initialised for code calrity but not needed, all the edges are
        # already set into previous blocks
        self.block17Edges = []

        self.block17 = Block.create_from_points(
            self.block17Points,
            self.block17Edges)
        self.block17.set_patch('top', 'outlet')
        self.block17.set_patch('bottom', 'wall')

        self.block17.chop(0, count=10, c2c_expansion=1)
        self.block17.chop(1, count=10, c2c_expansion=1)
        self.block17.chop(2, count=10, c2c_expansion=1)

    def build_block18(self):
        """
        Computes vertices and point in order to build block 18.
        """
        # TODO: Understand why the order has tob shifted
        # TODO: Correct the flawed logic
        # Computes vertices for block 18
        blk18v0 = self.capsuleSection5quad1[-1]
        blk18v1 = self.capsuleSection5quad2[-1]
        blk18v2 = self.capsuleSection5quad3[-1]
        blk18v3 = self.capsuleSection5quad4[-1]
        blk18v4 = self.ioSection5quad1[-1]  # [0]
        blk18v5 = self.ioSection5quad2[-1]  # [0]
        blk18v6 = self.ioSection5quad3[-1]  # [0]
        blk18v7 = self.ioSection5quad4[-1]  # [0]

        self.block18Points = [
            blk18v4, blk18v5, blk18v6, blk18v7,
            blk18v0, blk18v1, blk18v2, blk18v3,
        ]

        """
        self.init_3d_plot()
        print('block 18')
        colors = ['r','g','b','c','m','y','k','b']
        for pnt, c in zip(self.block18Points, colors):
            print(pnt)
            self.plot_points([pnt], c)
        """

        # Computes Theta and Phi for the inlet surface
        vInletFront = blk18v6
        self.show_plot()
        vInletFrontXY = [vInletFront[0], vInletFront[1], 0]
        print(f'vInletFront:   {vInletFront}')
        print(f'vInletFrontXY: {vInletFrontXY}')

        # WARNING: works only in the first quadrand!
        # Since eX = [1, 0, 0] the angle can be easily found in the following way
        theta = np.arccos(vInletFront[2] / self.c2)
        thetaOutlet = np.arccos(vInletFront[0]/(self.a1*np.sin(theta)))
        # Since eZ = [0, 0, 1] the angle can be easily found in the following way
        phiInlet = np.arccos(vInletFront[2] / self.c2)

        print(f'theta    [deg]: {np.rad2deg(theta)})')
        print(f'phiInlet [deg]: {np.rad2deg(phiInlet)}')

        # Theta: apearture from the x_axis on the XY plane
        N = 10
        theta = thetaOutlet
        thetaMin = +theta
        thetaMax = -theta
        angles = np.linspace(thetaMin, thetaMax, N)
        print(f'thetaOutlet [deg]: {thetaOutlet}')
        print(f'angles [deg]: {angles}')

        # Spline edges for the inlet surface
        c03 = []  # Spline points for edge 03
        c10 = []  # Spline points for edge 10
        c21 = []  # Spline points for edge 21
        c32 = []  # Spline points for edge 32

        for theta in angles:
            print(f'Theta [deg]: {theta}')
            if theta < 0:
                theta = 2*np.pi + theta
            centralPoint = self.point_on_ellipsoid(theta, phiInlet)

            tmp = self.rotQuatC1.apply(centralPoint)
            c03.append(tmp)
            tmp = self.rotQuatC2.apply(centralPoint)
            c10.append(tmp)
            tmp = self.rotQuatC3.apply(centralPoint)
            c21.append(tmp)
            tmp = self.rotQuatC4.apply(centralPoint)
            c32.append(tmp)

        self.block18Edges = [
            Edge(0, 3, c03),  # Spline Edge
            Edge(1, 0, c10),  # Spline Edge
            Edge(2, 1, c21),  # Spline Edge
            Edge(3, 2, c32),  # Spline Edge
        ]

        self.block18 = Block.create_from_points(
            self.block18Points,
            self.block18Edges)
        self.block18.set_patch('bottom', 'outlet')
        self.block18.set_patch('top', 'wall')

        self.block18.chop(0, count=10, c2c_expansion=1)
        self.block18.chop(1, count=10, c2c_expansion=1)
        self.block18.chop(2, count=10, c2c_expansion=1)

    # TODO: Put this part into the main for code clarity
    def mesh_3D(self):
        """
        Calls all the functions
        """
        mesh = Mesh()

        # Computes the upper points of a mesh section
        self.section_points()

        # Rotates the computed points in order to get 4 quadrants
        self.central_points_rotation()

        # Compute vertices and points for each blocks
        # Section 0
        self.build_block1()
        # Section 1
        self.build_block2()
        self.build_block3()
        self.build_block4()
        self.build_block5()
        # Section 2
        self.build_block6()
        self.build_block7()
        self.build_block8()
        self.build_block9()
        # Section 3
        self.build_block10()
        self.build_block11()
        self.build_block12()
        self.build_block13()
        # Section 4
        self.build_block14()
        self.build_block15()
        self.build_block16()
        self.build_block17()
        # Section 5
        self.build_block18()

        # Assembles the blocks
        # Section 0
        mesh.add_block(self.block1)
        # Section 1
        mesh.add_block(self.block2)
        mesh.add_block(self.block3)
        mesh.add_block(self.block4)
        mesh.add_block(self.block5)
        # Section 2
        mesh.add_block(self.block6)
        mesh.add_block(self.block7)
        mesh.add_block(self.block8)
        mesh.add_block(self.block9)
        # Section 3
        mesh.add_block(self.block10)
        mesh.add_block(self.block11)
        mesh.add_block(self.block12)
        mesh.add_block(self.block13)
        # Section 4
        mesh.add_block(self.block14)
        mesh.add_block(self.block15)
        mesh.add_block(self.block16)
        mesh.add_block(self.block17)
        # Section 5
        mesh.add_block(self.block18)

        mesh.write(output_path=os.path.join('case', 'system', 'blockMeshDict'))
        os.system('case/Allrun.mesh')

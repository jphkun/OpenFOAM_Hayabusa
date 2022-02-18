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

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

from scipy.spatial.transform import Rotation as R

from classy_blocks.classes.mesh import Mesh
from classy_blocks.classes.operations import Extrude, Face, Revolve
from classy_blocks.classes import operations
from classy_blocks.classes.block import Block
from classy_blocks.classes.primitives import Edge

from classy_blocks.util import functions as f


# TODO create a logger class
# TODO: Finish the circle function
# TODO: Document the cirlce function
# TODO: Test 5 configurations of the circle function
# TODO: Finish block 2
# TODO: Finish block 3
# TODO: Finish block 4
# TODO: Finish block 5

class Hayabusa:

    # TODO: Connect the capsuleDia input into the class
    # TODO: Verify that the main.py file has the correct input

    def __init__(self,
                 capsuleDia: float,
                 rNpoints: int,
                 thetaNpoints: int,
                 phiNpoints: int) -> None:
        """
        Init function

        Parameters
        ----------
        capsuleDia: float
            Capsules diameter, determines the capsule diameter and the mesh size
        rNpoints: int
            Number of mesh points in the "radial" direction
        thetaNpoints: int
            Number of mesh points in the "theta" direction
        phiNpoints: int
            Number of mesh points in the "phi" direction
        """
        self.capsuleDia = 0.4
        self.capsuleRad = self.capsuleDia/2
        self.capsuleSphereCenter = [self.capsuleDia/5, 0, 0]
        self.rNpoints = rNpoints
        self.thetaNpoints = thetaNpoints
        self.phiNpoints = phiNpoints

        # a1, b1, c1, are the "back" ellipsoid semi-axis
        self.a1 = 10 * self.capsuleDia
        self.b1 = 10 * self.capsuleDia
        self.c1 = 10 * self.capsuleDia

        # a2, b2, c2, are the "front" ellipsoid semi-axis
        self.a2 =  4 * self.capsuleDia
        self.b2 = 10 * self.capsuleDia
        self.c2 = 10  * self.capsuleDia


    # TODO: Wring type hinting output

    def point_on_sphere(self, rho: float, theta: float, phi: float, origin=[0,0,0]) -> None:
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

    def point_on_ellipsoid(self, theta: float, phi: float) -> None:
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

        if theta > np.pi/2 or theta > (3/2)*np.pi:
            x = self.a1 * np.cos(theta) * np.sin(phi)
            y = self.b1 * np.sin(theta) * np.sin(phi)
            z = self.c1 * np.cos(phi)
        elif theta > np.pi/2:
            x = self.a2 * np.cos(theta) * np.sin(phi)
            y = self.b2 * np.sin(theta) * np.sin(phi)
            z = self.c2 * np.cos(phi)
        else:
            print('Point_on_ellipsoid: Error with phi and theta, they are not in the correct range')
            sys.exit()

        point = [x,y,z]
        return point


    # TODO: Change the None in a vector

    def circle_on_sphere(self, rho:float, alpha:float, beta:float, gamma:float, t:float) -> None:
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

        point = np.array([x,y,z])
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
        point = [x,y,z]
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
        self.block1Vertices= [
            blk1v0, blk1v1, blk1v2, blk1v3, # Face 1
            blk1v4, blk1v5, blk1v6, blk1v7, # Face 2?
        ]

        # Arc defintions
        #######################################################################
        # WARNING: alpha is set to pi since the front of the capsule is set in
        # the negative x direction
        #######################################################################
        edge01pnt = self.point_on_sphere( self.capsuleRad, np.pi,  beta1,   sphereOg)
        edge12pnt = self.point_on_sphere( self.capsuleRad, alpha1, np.pi/2, sphereOg)
        edge23pnt = self.point_on_sphere( self.capsuleRad, np.pi,  beta3,   sphereOg)
        edge30pnt = self.point_on_sphere( self.capsuleRad, alpha3, np.pi/2, sphereOg)

        # Spline defintions
        N = 10 # Number of points for each spline defintion
        # Computation of the angles at which to compute the interpolation point
        # for the spline.
        angles45 = np.linspace(theta4, theta5, N),
        edge45pnts = [self.point_on_ellipsoid( i, phi4) for i in angles45[0]]

        angles56 = np.linspace(phi5, phi6, N),
        edge56pnts = [self.point_on_ellipsoid( theta5, i) for i in angles56[0]]

        angles67 = np.linspace(theta6, theta7, N),
        edge67pnts = [self.point_on_ellipsoid( i, phi6) for i in angles67[0]]

        angles74 = np.linspace(phi7, phi4, N),
        edge74pnts = [self.point_on_ellipsoid( theta7, i) for i in angles74[0]]

        self.block1Edges = [
            Edge(0, 1, edge01pnt), # Arc edge
            Edge(1, 2, edge12pnt), # Arc edge
            Edge(2, 3, edge23pnt), # Arc edge
            Edge(3, 0, edge30pnt), # Arc edge
            Edge(4, 5, edge45pnts), # Spline edge
            Edge(5, 6, edge56pnts), # Spline edge
            Edge(6, 7, edge67pnts), # Spline edge
            Edge(7, 4, edge74pnts), # Spline edge
        ]

        self.block1 = Block.create_from_points(self.block1Vertices, self.block1Edges)
        self.block1.set_patch('top','inlet')
        self.block1.set_patch('bottom','wall')

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
            print(pnt)
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

    # TODO: On traversera le pont rendu à la rivière
    def point_vertical_back():
        pass

    # TODO: Documentation
    # TODO: Logic semplification
    def section_points(self):
        """

        """
        # TODO: Generalize the sphere origin
        sphereOg = [0.08, 0, 0]
        print(2*'\n')
        # New code
        #######################################################################
        # WARNING all the following points are in the negative x direction
        # hence add 180 deg for the correct point description in spherical
        # coordinates.
        #######################################################################
        # Capsules points bottom
        # Angle between the XY plane the sphere origin and vertex 3 of the block 1
        block2AlphaBot1 = np.arcsin(1/4)
        block2AlphaBot2 = np.pi + block2AlphaBot1
        print(f'block2AlphaBot  [deg]: {np.rad2deg(block2AlphaBot1):.4f}')
        print(f'block2ThetaBot  [deg]: {np.rad2deg(block2AlphaBot2):.4f}')
        # Angle between the XZ plane the sphere origin and vertex 3 of the block 1
        block2BetaBot1 = np.arcsin(1/4)
        block2BetaBot2 = np.pi/2 - block2BetaBot1 # To be used in point_on_sphere
        print(f'block2BetaBot   [deg]: {np.rad2deg(block2BetaBot1):.4f}')
        print(f'block2PhiBot    [deg]: {np.rad2deg(block2BetaBot2):.4f}')

        # Inlet points
        # Angle between the XY plane the MESH origin and vertex 7? 4? of the block 1
        block2ThetaBot1 = self.angle_conversion(block2AlphaBot1)
        block2ThetaBot2 = np.pi + block2ThetaBot1
        print(f'block2ThetaBot1 [deg]: {np.rad2deg(block2ThetaBot1):.4f}')
        print(f'block2ThetaBot2 [deg]: {np.rad2deg(block2ThetaBot2):.4f}')
        # Angle between the XZ plane the MESH origin and vertex 7? 4? of the block 1
        block2PhiBot1 = self.angle_conversion(block2BetaBot1)
        block2PhiBot2 = np.pi/2 - block2PhiBot1 # To use in point_on_ellipsoid
        print(f'block2PhiBot1   [deg]: {np.rad2deg(block2PhiBot1):.4f}')
        print(f'block2PhiBot2   [deg]: {np.rad2deg(block2PhiBot2):.4f}')

        # Capsule points top
        # Angle between teh XY plane, the sphere origin, and the center point
        # of the edge between vertex 1 and vertex 2.
        block2BetaTop1 = np.pi - 2 * np.arctan(2)
        block2BetaTop2 = np.pi - 2 * np.arctan(3)
        block2PhiTop1 = np.pi/2 - block2BetaTop1
        block2PhiTop2 = np.pi/2 - block2BetaTop2
        print(f'block2beta1     [deg]: {np.rad2deg(block2BetaTop1):.4f}')
        print(f'block2Beta2     [deg]: {np.rad2deg(block2BetaTop2):.4f}')
        print(f'block2Phi1      [deg]: {np.rad2deg(block2PhiTop1):.4f}')
        print(f'block2Phi2      [deg]: {np.rad2deg(block2PhiTop2):.4f}')

        # Computation of 10 angles "central point" between the edges01 and edges 45
        N = 10
        anglesCapsule0 = np.linspace(np.pi/2, block2BetaBot2, N)
        anglesCapsule1 = np.linspace(block2BetaTop1, block2BetaBot2, N) # From "Top" to "Bottom"
        print(f'anglesCapsule0  [deg]:\n {np.rad2deg(anglesCapsule0)}')
        print(f'anglesCapsule1  [deg]:\n {np.rad2deg(anglesCapsule1)}')

        anglesInlet0 = np.linspace(block2PhiBot2, np.pi/2, N)
        anglesInlet1 = np.linspace(block2PhiTop2, block2PhiBot2, N)
        anglesInlet2 = np.linspace(0, block2PhiTop2, N)
        print(f'anglesInlet    [deg]:\n {np.rad2deg(anglesInlet0)}')
        print(f'anglesInlet    [deg]:\n {np.rad2deg(anglesInlet1)}')
        print(f'anglesInlet    [deg]:\n {np.rad2deg(anglesInlet2)}')

        # Central points on wall/capsule

        # Front of the capsule
        # Section 0
        self.ccPnts0 = [self.point_on_sphere( self.capsuleRad, np.pi, ang, sphereOg) for ang in anglesCapsule0]
        # Section 1
        self.ccPnts1 = [self.point_on_sphere( self.capsuleRad, np.pi, ang, sphereOg) for ang in anglesCapsule1]
        # Section 2
        xFront = np.linspace(self.ccPnts1[0][0], 0, N)
        self.ccPnts2 = [self.point_on_straight_front(x) for x in xFront]

        # Back of the capsule
        # Section 3
        xBack  = np.linspace(0, sphereOg[0], N)
        self.ccPnts3 = [self.point_on_straight_back(x) for x in xBack]
        # TODO:
        self.ccPnts4 = []
        self.ccPnts5 = []

        """
        plt.figure()
        print(3*'\n')
        print([x[0] for x in ccPnts0])
        print(ccPnts0[0][:])
        plt.plot([x[0] for x in ccPnts0], [x[2] for x in ccPnts0], 'o')
        plt.plot([x[0] for x in ccPnts1], [x[2] for x in ccPnts1], 'o')
        plt.plot([x[0] for x in ccPnts2], [x[2] for x in ccPnts2], 'o')
        plt.plot([x[0] for x in ccPnts3], [x[2] for x in ccPnts3], 'o')
        plt.axis('equal')
        plt.show()
        """

        # Outside central points, on external mesh
        # Front points
        # Section 0
        self.iocPnts0 = [self.point_on_ellipsoid(np.pi, ang) for ang in anglesInlet0]
        # Section 1
        self.iocPnts1 = [self.point_on_ellipsoid(np.pi, ang) for ang in anglesInlet1]
        # Section 2
        self.iocPnts2 = [self.point_on_ellipsoid(np.pi, ang) for ang in anglesInlet2]
        # Back of the capsule
        # TODO: compute the final phi angle for section 3
        # Section 3
        self.iocPnts3 = []
        self.iocPnts4 = []
        self.iocPnts5 = []

        """
        self.init_3d_plot()
        self.plot_points(ccPnts0,'r')
        self.plot_points(ccPnts1,'g')
        self.plot_points(ccPnts2,'b')
        self.plot_points(ccPnts3,'c')
        self.plot_points(ccPnts4,'m')
        self.show_plot()
        """
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
        quaternionRotation1 = self.rotation_vector(eX, quadrant1Angle)
        quaternionRotationC1= self.rotation_vector(eX, quadrant1cAngle)
        quaternionRotation2 = self.rotation_vector(eX, quadrant2Angle)
        quaternionRotationC2= self.rotation_vector(eX, quadrant2cAngle)
        quaternionRotation3 = self.rotation_vector(eX, quadrant3Angle)
        quaternionRotationC3= self.rotation_vector(eX, quadrant3cAngle)
        quaternionRotation4 = self.rotation_vector(eX, quadrant4Angle)
        quaternionRotationC4= self.rotation_vector(eX, quadrant4cAngle)

        # Quaternion self.rotation functions defined
        self.rotQuat1 = R.from_quat(quaternionRotation1)
        self.rotQuatC1= R.from_quat(quaternionRotationC1)
        self.rotQuat2 = R.from_quat(quaternionRotation2)
        self.rotQuatC2= R.from_quat(quaternionRotationC2)
        self.rotQuat3 = R.from_quat(quaternionRotation3)
        self.rotQuatC3= R.from_quat(quaternionRotationC3)
        self.rotQuat4 = R.from_quat(quaternionRotation4)
        self.rotQuatC4= R.from_quat(quaternionRotationC4)

        # Saves the self.rotated capsule construction points
        # Section 0
        self.capsuleSection0quad1 = self.rotQuat1.apply(self.ccPnts0)
        self.capsuleSection0quad2 = self.rotQuat2.apply(self.ccPnts0)
        self.capsuleSection0quad3 = self.rotQuat3.apply(self.ccPnts0)
        self.capsuleSection0quad4 = self.rotQuat4.apply(self.ccPnts0)
        # Section 1
        self.capsuleSection1quad1 = self.rotQuat1.apply(self.ccPnts1)
        self.capsuleSection1quadC1= self.rotQuatC1.apply(self.ccPnts1)
        self.capsuleSection1quad2 = self.rotQuat2.apply(self.ccPnts1)
        self.capsuleSection1quadC2= self.rotQuatC2.apply(self.ccPnts1)
        self.capsuleSection1quad3 = self.rotQuat3.apply(self.ccPnts1)
        self.capsuleSection1quadC3= self.rotQuatC3.apply(self.ccPnts3)
        self.capsuleSection1quad4 = self.rotQuat4.apply(self.ccPnts1)
        self.capsuleSection1quadC4= self.rotQuatC4.apply(self.ccPnts1)
        # Section 2
        self.capsuleSection2quad1 = self.rotQuat1.apply(self.ccPnts2)
        self.capsuleSection2quad2 = self.rotQuat2.apply(self.ccPnts2)
        self.capsuleSection2quad3 = self.rotQuat3.apply(self.ccPnts2)
        self.capsuleSection2quad4 = self.rotQuat4.apply(self.ccPnts2)
        # Section 3
        self.capsuleSection3quad1 = self.rotQuat1.apply(self.ccPnts3)
        self.capsuleSection3quad2 = self.rotQuat2.apply(self.ccPnts3)
        self.capsuleSection3quad3 = self.rotQuat3.apply(self.ccPnts3)
        self.capsuleSection3quad4 = self.rotQuat4.apply(self.ccPnts3)
        # Section 4
        # Section 5

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
        self.ioSection2quad2 = self.rotQuat2.apply(self.iocPnts2)
        self.ioSection2quad3 = self.rotQuat3.apply(self.iocPnts2)
        self.ioSection2quad4 = self.rotQuat4.apply(self.iocPnts2)
        # Section 3
        # Section 4
        # Section 5

    def compute_angles(self, p1, e):
        """
        Returns
        -------
        angle: float
            angle between the point and the vector e
        """

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
        # First we need the angle
        vCapsuleFront  = blk1v1 - self.capsuleSphereCenter
        vCapsuleFrontXY = [vCapsuleFront[0], vCapsuleFront[1], 0]
        # Since eX = [1, 0, 0] the angle can be easily found in the following way
        alphaCapsule = np.arccos(vCapsuleFrontXY[0] / np.linalg.norm(vCapsuleFrontXY))
        # Since eZ = [0, 0, 1] the angle can be easily found in the following way
        betaCapsule   = np.arccos(vCapsuleFront[2] / np.linalg.norm(vCapsuleFront))

        # Computes Theta and Phi for the inlet surface
        vInletFront   = blk1v5
        vInletFrontXY = [vInletFront[0], vInletFront[1], 0]
        # Since eX = [1, 0, 0] the angle can be easily found in the following way
        thetaInlet = np.arccos(vInletFrontXY[0] / np.linalg.norm(vInletFrontXY))
        # Since eZ = [0, 0, 1] the angle can be easily found in the following way
        phiInlet   = np.arccos(vInletFront[2]   / np.linalg.norm(vInletFront))

        # TODO: Implement a test for this
        """
        # Test to see if the angles found in spherical coordiantes are correct
        p = self.point_on_sphere(
            rho   = self.capsuleRad,
            theta = alphaCapsule,
            phi   = betaCapsule,
            origin= self.capsuleSphereCenter)
        print(f'p:             {p}')
        print(f'vCapsuleFront: {blk1v2}')

        p = self.point_on_ellipsoid(thetaInlet, phiInlet)
        print(f'p:  {p}')
        print(f'vInletBack:    {vInletBack}')
        """

        N = 10
        theta  = np.pi - thetaInlet
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
        c5 = [] # Spline points for ege 56
        c6 = [] # Spline points for ege 67
        c7 = [] # Spline points for ege 74
        c8 = [] # Spline points for ege 45

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
            Edge(0, 1, c2), # Arc Edge
            Edge(1, 2, c3), # Arc Edge
            Edge(2, 3, c4), # Arc Edge
            Edge(3, 0, c1), # Arc Edge
            Edge(7, 4, c5), # Spline Edge
            Edge(4, 5, c6), # Spline Edge
            Edge(5, 6, c7), # Spline Edge
            Edge(6, 7, c8), # Spline Edge
        ]

        self.block1 = Block.create_from_points(self.block1Points, self.block1Edges)
        self.block1.set_patch('top','inlet')
        self.block1.set_patch('bottom','wall')

        self.block1.chop(0, count=10, c2c_expansion=1)
        self.block1.chop(1, count=10, c2c_expansion=1)
        self.block1.chop(2, count=10, c2c_expansion=1)
        self.show_plot()


    def build_block2(self):
        """

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
            Edge(1, 2, self.capsuleSection1quad4[4]),  # Arc Edge
            Edge(2, 3, self.capsuleSection1quadC1[0]), # Arc Edge
            Edge(3, 0, self.capsuleSection1quad1[4]),  # Arc Edge
            Edge(5, 6, self.ioSection1quad4[::-1]),    # Spline Edge
            Edge(6, 7, self.ioSection1quadC1[0]),      # Spline Edge
            Edge(7, 4, self.ioSection1quad1), # Spline Edge
        ]

        self.block2 = Block.create_from_points(self.block2Points, self.block2Edges)
        self.block2.set_patch('top','inlet')
        self.block2.set_patch('bottom','wall')

        self.block2.chop(0, count=10, c2c_expansion=1)
        self.block2.chop(1, count=10, c2c_expansion=1)
        self.block2.chop(2, count=10, c2c_expansion=1)
        self.show_plot()

    # TODO: Put this part into the main for code clarity
    def mesh_3D(self):
        """

        """
        mesh = Mesh()
        # Computes the upper points of a mesh section
        self.section_points()
        # Rotates the computed points in order to get 4 quadrants
        self.central_points_rotation()
        # Builds block1
        self.build_block1()
        # Builds block2
        self.build_block2()
        # Assembles the blocks
        mesh.add_block(self.block1)
        mesh.add_block(self.block2)

        #sys.exit()

        mesh.write(output_path=os.path.join('case','system','blockMeshDict'))
        os.system('case/Allrun.mesh')

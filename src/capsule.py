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

        point = np.array([x,y,z])
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

    def test_circle(self):
        p = []
        x = []
        y = []
        z = []
        X = []
        Y = []
        Z = []
        alpha = np.deg2rad(2*2*np.arctan(2))
        theta = np.deg2rad(180)
        phi = np.deg2rad(90)
        N = 21
        for i in range(0,360,N):
            angle = np.deg2rad(i)
            point = self.circle_on_sphere(1, alpha, theta, phi, angle)
            p.append(point)
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
            X.append((i/180)-1)
            Y.append((i/180)-1)
            Z.append((i/180)-1)

        p1 = self.circle_on_sphere(1, alpha, theta, phi, np.deg2rad( 45))
        p2 = self.circle_on_sphere(1, alpha, theta, phi, np.deg2rad(-45))
        print(' ')
        print(X)
        print(np.zeros(N-3)+X)
        plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(x, y, z)
        ax.scatter(p1[0],p1[1],p1[2],'or')
        ax.scatter(p2[0],p2[1],p2[2],'ob')
        ax.scatter(X, np.zeros(N-3), np.zeros(N-3))
        ax.scatter(np.zeros(N-3), X, np.zeros(N-3))
        ax.scatter(np.zeros(N-3), np.zeros(N-3), X)

        plt.show()


    # TODO: Documentation
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

    @staticmethod
    def angle_conversion(alpha):
        """
        converts the angle form one reference frame to the other
        """
        b = 200
        c = 80
        a = np.sqrt(b**2 + c**2 - 2*b*c*np.cos(alpha))
        theta = np.arcsin(b * np.sin(alpha) / a)
        return theta


    # TODO: make the return block documentation correct
    def block1_3d(self):
        """
        Constructs block 1 and return a block with edges.
        """
        rho = self.capsuleDia/2
        sphereOg = [0.08, 0, 0]
        block1Alpha = np.arcsin(1/4)
        block1Beta  = np.arcsin(1/4)
        block1Theta = self.angle_conversion(block1Alpha)
        block1Phi   = self.angle_conversion(block1Beta)

        # Block point definition
        self.block1Points = [
            self.point_on_sphere(rho, np.pi + block1Alpha, np.pi/2 - block1Beta, sphereOg),
            self.point_on_sphere(rho, np.pi - block1Alpha, np.pi/2 - block1Beta, sphereOg),
            self.point_on_sphere(rho, np.pi - block1Alpha, np.pi/2 + block1Beta, sphereOg),
            self.point_on_sphere(rho, np.pi + block1Alpha, np.pi/2 + block1Beta, sphereOg),
            self.point_on_ellipsoid( np.pi + block1Theta,  np.pi/2 - block1Phi),
            self.point_on_ellipsoid( np.pi - block1Theta,  np.pi/2 - block1Phi),
            self.point_on_ellipsoid( np.pi - block1Theta,  np.pi/2 + block1Phi),
            self.point_on_ellipsoid( np.pi + block1Theta,  np.pi/2 + block1Phi),
        ]

        for point in self.block1Points:
            print(point)

        # Edge defintions
        angle45 = np.linspace(np.pi   + block1Phi,   np.pi   - block1Phi,   10), # Spline Edge
        angle56 = np.linspace(np.pi/2 - block1Theta, np.pi/2 + block1Theta, 10), # Spline Edge
        angle67 = np.linspace(np.pi   - block1Phi,   np.pi   + block1Phi,   10), # Spline Edge
        angle74 = np.linspace(np.pi/2 + block1Theta, np.pi/2 - block1Theta, 10), # Spline Edge

        self.block1Edges = [
            Edge(0, 1, self.point_on_sphere(rho, np.pi,             np.pi/2-block1Beta, sphereOg)), # Arc Edge
            Edge(1, 2, self.point_on_sphere(rho, np.pi-block1Alpha , np.pi/2           , sphereOg)), # Arc Edge
            Edge(2, 3, self.point_on_sphere(rho, np.pi,             np.pi/2+block1Beta, sphereOg)), # Arc Edge
            Edge(3, 0, self.point_on_sphere(rho, np.pi+block1Alpha, np.pi/2           , sphereOg)), # Arc Edge
            Edge(4, 5, [self.point_on_ellipsoid( i, np.pi/2-block1Phi ) for i in angle45[0]]), # Spline Edge
            Edge(5, 6, [self.point_on_ellipsoid( np.pi-block1Theta, i ) for i in angle56[0]]), # Spline Edge
            Edge(6, 7, [self.point_on_ellipsoid( i, np.pi/2+block1Phi ) for i in angle67[0]]), # Spline Edge
            Edge(7, 4, [self.point_on_ellipsoid( np.pi+block1Theta, i ) for i in angle74[0]]), # Spline Edge
        ]

        self.block1 = Block.create_from_points(self.block1Points, self.block1Edges)
        #block.set_patch(['left','right','front','back'], '')
        self.block1.set_patch('top','inlet')
        self.block1.set_patch('bottom','wall')

        self.block1.chop(0, count=10, c2c_expansion=1)
        self.block1.chop(1, count=10, c2c_expansion=1)
        self.block1.chop(2, count=10,  c2c_expansion=1)


    def block2_3d(self):
        """

        """

        rho = self.capsuleDia/2
        sphereOg = [0.08, 0, 0]
        block1Alpha = np.arcsin(1/4)
        block1Beta  = np.arcsin(1/4)
        block1gamma = 2 * np.arctan(2)
        block1delta = 2 * np.arctan(2)
        block1Theta = self.angle_conversion(block1Alpha)
        block1Phi   = self.angle_conversion(block1Beta)
        # Aperture angle of the cone
        alpha = 2*np.arctan(2) - np.pi/2
        beta = np.deg2rad(90)
        gamma = np.deg2rad(180)
        t1 = np.deg2rad(-45+180)
        t2 = np.deg2rad(45+180)

        alpha2 = self.angle_conversion(alpha)
        # TODO: change the point_on_sphere function by circle_on_sphere
        self.block2Points = [
            self.circle_on_sphere(rho, alpha, beta, gamma, t1) + sphereOg,
            self.circle_on_sphere(rho, alpha, beta, gamma, t2) + sphereOg,
            self.point_on_sphere(rho, np.pi - block1Alpha, np.pi/2 - block1Beta, sphereOg),
            self.point_on_sphere(rho, np.pi + block1Alpha, np.pi/2 - block1Beta, sphereOg),

            self.point_on_ellipsoid((5/4)*np.pi, alpha2),
            self.point_on_ellipsoid((3/4)*np.pi, alpha2),
            self.point_on_ellipsoid( np.pi - block1Theta,  np.pi/2 - block1Phi),
            self.point_on_ellipsoid( np.pi + block1Theta,  np.pi/2 - block1Phi),
        ]
        for point in self.block2Points:
            print(point)

        # angle on the XZ plane between the leading edge of block2 and the
        # trailing edge of block2 (leading edge of block 6 7 8 or 9 TBD)
        angles = np.linspace(alpha2, alpha ,3)
        center_points = []
        radiuses = []
        edge_1 = []
        edge_2 = []

        for ang in angles:
            p = self.point_on_ellipsoid(np.deg2rad(180), ang)
            r = np.linalg.norm(p[1:])
            radiuses.append(r)
            c = [p[0], 0, 0]
            center_points.append(c)
            pe1 = self.point_circle_xaxis(c, r, np.deg2rad(45))
            edge_1.append(pe1)

            pe2 = self.point_circle_xaxis(c, r, np.deg2rad(-45))
            edge_2.append(pe2)

        self.block2Edges = [
            Edge(0,1, self.circle_on_sphere(rho,     alpha, beta, gamma, np.deg2rad(180))+sphereOg), # Arc Edge
            Edge(1,2, self.circle_on_sphere(rho,0.75*alpha, beta, gamma, np.deg2rad(180+45))+sphereOg), # Arc Edge
            Edge(3,0, self.circle_on_sphere(rho,0.75*alpha, beta, gamma, np.deg2rad(180-45))+sphereOg), # Arc Edge
            Edge(4,5, self.point_circle_xaxis(center_points[-1], r, 0)), # Arc Edge
            Edge(5,6, [edge_1[1]]), # Spline Edge
            Edge(7,4, [edge_2[1]]), # Spline Edge
        ]
        self.block2 = Block.create_from_points(self.block2Points, self.block2Edges)
        self.block2.set_patch('top','inlet')
        self.block2.set_patch('bottom','wall')

        self.block2.chop(0, count=10, c2c_expansion=1)
        self.block2.chop(1, count=10, c2c_expansion=1)
        self.block2.chop(2, count=10, c2c_expansion=1)


    def mesh_3D(self):
        """

        """
        mesh = Mesh()
        self.block1_3d()
        self.block2_3d()
        mesh.add_block(self.block1)
        mesh.add_block(self.block2)
        #self.test_circle()
        #sys.exit()
        mesh.write(output_path=os.path.join('case','system','blockMeshDict'))
        os.system('case/Allrun.mesh')

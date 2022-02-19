import unittest
import numpy as np
import src.capsule as cpsl

class Test_capsule(unittest.TestCase):

    capsule = cpsl.Hayabusa(0.2,1,1,1)

    def test_point_on_sphere_0_0_0(self):
        rho = 0
        theta = np.deg2rad(0)
        phi = np.deg2rad(0)

        point = self.capsule.point_on_sphere(rho, theta, phi)

        solution = [0,0,0]
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)


    def test_point_on_sphere_1_0_0(self):
        rho = 1
        theta = np.deg2rad(0)
        phi = np.deg2rad(0)

        point = self.capsule.point_on_sphere(rho, theta, phi)

        solution = [0,0,1]
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)


    def test_point_on_sphere_1_0_90(self):
        rho = 1
        theta = np.deg2rad(0)
        phi = np.deg2rad(90)

        point = self.capsule.point_on_sphere(rho, theta, phi)
        solution = [1,0,0]
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)


    def test_point_on_sphere_1_90_90(self):
        rho = 1
        theta = np.deg2rad(90)
        phi = np.deg2rad(90)

        point = self.capsule.point_on_sphere(rho, theta, phi)
        solution = [0,1,0]
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)


    def test_point_on_sphere_1_180_90(self):
        rho = 1
        theta = np.deg2rad(180)
        phi = np.deg2rad(90)

        point = self.capsule.point_on_sphere(rho, theta, phi)
        solution = [-1,0,0]
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)


    def test_circle_on_sphere_0_0_0_0_0(self):
        rho = 0
        alpha = 0
        beta = 0
        gamma = 0
        t = 0

        point = self.capsule.circle_on_sphere(rho, alpha, beta, gamma, t)
        solution = np.array([0.,0.,0.])
        #print(' ')
        #print(point)
        #print(solution)
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)

    def test_circle_on_sphere_1_0_0_0_0(self):
        rho = 1
        alpha = 0
        beta = 0
        gamma = 0
        t = 0

        point = self.capsule.circle_on_sphere(rho, alpha, beta, gamma, t)
        solution = np.array([0.,0.,1.])

        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)


    def test_circle_on_sphere_1_90_0_0_0(self):
        rho = 1
        alpha = np.deg2rad(90)
        beta = 0
        gamma = 0
        t = 0
        point = self.capsule.circle_on_sphere(rho, alpha, beta, gamma, t)
        solution = np.array([1.,0.,0.])
        #print(' ')
        #print(point)
        #print(solution)
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)

    def test_circle_on_sphere_1_90_0_0_90(self):
        rho = 1
        alpha = np.deg2rad(90)
        beta = 0
        gamma = 0
        t = np.deg2rad(90)
        point = self.capsule.circle_on_sphere(rho, alpha, beta, gamma, t)
        solution = np.array([0.,1.,0.])
        #print(' ')
        #print(point)
        #print(solution)
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)

    def test_circle_on_sphere_1_90_0_0_180(self):
        rho = 1
        alpha = np.deg2rad(90)
        beta = 0
        gamma = 0
        t = np.deg2rad(180)
        point = self.capsule.circle_on_sphere(rho, alpha, beta, gamma, t)
        solution = np.array([-1.,0.,0.])
        #print(' ')
        #print(point)
        #print(solution)
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)

    def test_circle_on_sphere_1_90_90_180_0(self):
        rho = 1
        alpha = np.deg2rad( 90)
        beta  = np.deg2rad( 90)
        gamma = np.deg2rad(180)
        t     = np.deg2rad(  0)
        point = self.capsule.circle_on_sphere(rho, alpha, beta, gamma, t)
        solution = np.array([0.,0.,-1.])
        #print(' ')
        #print(point)
        #print(solution)
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)

    def test_circle_on_sphere_1_90_90_180_180(self):
        rho = 1
        alpha = np.deg2rad( 90)
        beta  = np.deg2rad( 90)
        gamma = np.deg2rad(180)
        t     = np.deg2rad(180)
        point = self.capsule.circle_on_sphere(rho, alpha, beta, gamma, t)
        solution = np.array([0.,0.,1.])
        #print(' ')
        #print(point)
        #print(solution)
        for coord, sol in zip(point, solution):
            self.assertAlmostEqual(coord, sol, 5)

    def test_point_on_ellipsoid(self):
        """
        Tests if theta is selected correctly
        """
        # Test angles
        theta_000_090 = np.linspace(  np.deg2rad(0),  np.deg2rad(90), 91)
        theta_091_180 = np.linspace( np.deg2rad(91), np.deg2rad(180), 90)
        theta_181_270 = np.linspace(np.deg2rad(181), np.deg2rad(270), 90)
        theta_271_360 = np.linspace(np.deg2rad(271), np.deg2rad(360), 90)

        #print(' ')
        #print(theta_000_090)
        #print(theta_091_180)
        #print(theta_181_270)
        #print(theta_271_360)

        for ang in theta_000_090:
            point, cadran = self.capsule.point_on_ellipsoid(ang, 0)
            self.assertEqual(cadran, 'Front')

        for ang in theta_091_180:
            point, cadran = self.capsule.point_on_ellipsoid(ang, 0)
            self.assertEqual(cadran, 'Back')

        for ang in theta_181_270:
            point, cadran = self.capsule.point_on_ellipsoid(ang, 0)
            self.assertEqual(cadran, 'Back')

        for ang in theta_271_360:
            point, cadran = self.capsule.point_on_ellipsoid(ang, 0)
            self.assertEqual(cadran, 'Front')


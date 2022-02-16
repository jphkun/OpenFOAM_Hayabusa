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

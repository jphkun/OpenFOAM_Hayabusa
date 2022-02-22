from src.capsule import Hayabusa

def main():
    mesh = Hayabusa(
        capsuleDia=0.2,
        rNpoints=10,
        thetaNpoints=10,
        phiNpoints=10
    )
    mesh.mesh_3D()


if __name__ == "__main__":
    main()

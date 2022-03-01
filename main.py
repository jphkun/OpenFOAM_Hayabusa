import argparse
from typing import Union, Dict, Optional, Sequence
from src.capsule import Hayabusa

def cli(argv: Optional[Sequence[str]] = None ) -> dict[str,object]:
    """
    Command line interface, permets to the user to choose which type of mesh he
    desires.

    Parameter
    --------
    None

    Return
    ------
    meshParameters: dict
        Dictionnary wich conatinas the necessary information to compute the mesh
    """
    parser = argparse.ArgumentParser()
    # TODO: General helper and code description

    # Mesh dimensions
    parser.add_argument(
        '-d', '--mesh-dimensions',
        default = '3D',
        type = str,
        help = 'Number of mesh dimensions could be 2D or 3D'
    )
    # Capsule diameter
    parser.add_argument(
        '-dia', '--capsule-diameter',
        default = 0.4,
        type = float,
        help = 'Number of mesh dimensions could be 2D or 3D'
    )

    # Number of points
    parser.add_argument(
        '-r', '-radial-number-points',
        default = 10,
        type = int,
        help = 'Number of mesh points in the radial direction'
    )
    parser.add_argument(
        '-t', '-theta-number-points',
        default = 10,
        type = int,
        help = 'Number of mesh points in the theta direction'
    )
    parser.add_argument(
        '-p', '-phi-number-points',
        default = 10,
        type = int,
        help = 'If the mesh is 2D, will be ignored'
    )

    # Inflation layer parameter
    parser.add_argument(
        '-i', '-inflation-layer-parameter',
        default = 1.2,
        type = float,
        help = 'Infaltion layer parameter. Always in the radial diection'
    )

    # Retreives all command line inputs
    args = parser.parse_args()

    # Spherical coordinates use the mathematical convention in this case
    meshParameters = {
        'dimentions':           args.mesh_dimensions,
        'capsuleDiameter':      args.capsule_diameter,
        'radialNumberOfPoints': args.r,
        'thetaNumberOfPoints':  args.t,
        'phiNumberOfPoints':    args.p,
        'inflationLayerParam':  args.i,
    }
    return meshParameters

def main(inputParameters:dict) -> None:
    """
    Main function that give an overview of the code architecture

    Parameter
    ---------
    inputParameters: dict
        Input parameters from the command line, contains all the necessary
        information to compute the mesh.

    Return
    ------
    """
    mesh = Hayabusa(
        dimentions =   inputParameters['dimentions'],
        capsuleDia =   inputParameters['capsuleDiameter'],
        rNpoints =     inputParameters['radialNumberOfPoints'],
        thetaNpoints = inputParameters['thetaNumberOfPoints'],
        phiNpoints =   inputParameters['phiNumberOfPoints'],
        inflation =    inputParameters['inflationLayerParam'],
    )
    mesh.mesh_3D()

if __name__ == "__main__":
    meshParameters = cli()
    main(meshParameters)

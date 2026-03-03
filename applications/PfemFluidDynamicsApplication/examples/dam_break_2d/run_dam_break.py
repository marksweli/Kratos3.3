"""
Run script for the 2D Dam Break PFEM example.

This is the first PFEM case in GID: a 2D dam break simulation using a
Bingham fluid (density 1600 kg/m³, dynamic viscosity 300 Pa·s, yield
stress 50 Pa) in a rectangular container (0.64 m × 0.40 m).

The initial water column occupies the left half of the container
(width 0.32 m, height 0.267 m).

Usage
-----
Run from the directory containing this script (important, because the
solver resolves file paths relative to the current working directory):

    cd applications/PfemFluidDynamicsApplication/examples/dam_break_2d
    python run_dam_break.py

Prerequisites
-------------
KratosMultiphysics must be compiled and its Python bindings must be
importable. Follow the build instructions in the repository INSTALL.md.
"""

import os
import sys

# Ensure that imports work when the script is executed directly.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import KratosMultiphysics
from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis import (
    PfemFluidDynamicsAnalysis,
)


def run_simulation(parameter_file_name: str = "ProjectParameters.json") -> None:
    """
    Read the project parameters and run the PFEM fluid dynamics analysis.

    Parameters
    ----------
    parameter_file_name:
        Path to the JSON project parameters file (relative to the
        current working directory).
    """
    with open(parameter_file_name, "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = PfemFluidDynamicsAnalysis(model, parameters)
    simulation.Run()


if __name__ == "__main__":
    # Change to the script directory so that all relative paths in
    # ProjectParameters.json resolve correctly.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    run_simulation()

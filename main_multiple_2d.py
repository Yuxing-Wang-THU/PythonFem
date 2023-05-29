from fem_2d import run_fem_simulation_multiple
from matplotlib import pyplot as plt
from matplotlib.patches import Wedge
import seaborn as sns
import numpy as np
import argparse

sns.set_style("darkgrid")
# Rectangle obstacle
rectangle_pos = [(4,0),(2,1)]
rectangle = plt.Rectangle(rectangle_pos[0],*rectangle_pos[1],facecolor='black')

# Semicircle obstacle
ball_pos, ball_radius = np.array([10, 0.0]), 2
semicircle=Wedge((ball_pos[0],ball_pos[1]),ball_radius,0,180,facecolor='green')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')

    # Parameters of Afine transformation
    parser.add_argument('--voxelx', type=int, default=8,
                        help='Number of voxels (x-axis)')
    parser.add_argument('--voxely', type=int, default=3,
                        help='Number of voxels (y-axis)')
    parser.add_argument('--em', type=list, default=[2,3,6,7,10,11,14,15],
                        help='Empty voxels:[],[2,3,6,7,10,11,14,15],[2,3,6,7]')
    
    # Parameters of FEM
    parser.add_argument('--mode', type=str, default="fall",
                        help='Mode: bend,fall')
    parser.add_argument('--model', type=int, default=3,
                        help='Constitutive model, 0:Linear, 1:STVK, 2:Co-rotated, 3: Neohookean')
    parser.add_argument('--Y', type=int, default=10000,
                        help='Youngs modulus')
    parser.add_argument('--V', type=float, default=0.1,
                        help='Poisson ratio, [-1, 0.5]')
    parser.add_argument('--dt', type=float, default=0.001,
                        help='Time horizon')
    parser.add_argument('--max_iter', type=int, default=200,
                        help='Iteration')
    parser.add_argument('--mass', type=float, default=1.0,
                        help='Mass of the point') 
    parser.add_argument('--implict', action="store_true",
                        help='Implict Euler Method, for STVK and Neohookean')                               
    args = parser.parse_args()

    args.rectangle = rectangle
    args.semicircle = semicircle
    args.ball_pos = ball_pos
    args.ball_radius = ball_radius
    args.rectangle_pos = rectangle_pos

    run_fem_simulation_multiple(args)


import argparse
import numpy as np
from fem_2d import run_fem_simulation_single

# Init position of the triangular
ref_node_pos = np.array([[0,0],[1,0],[0,1]], dtype=float)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')

    # Parameters of Afine transformation
    parser.add_argument('--theta', type=float, default=-25.0,
                        help='Counterclockwise rotation angle')
    parser.add_argument('--tx', type=float, default=3.0,
                        help='Movement along the x-axis')
    parser.add_argument('--ty', type=float, default=3.0,
                        help='Movement along the y-axis')
    parser.add_argument('--sx', type=float, default=2,
                        help='Stretching along the x-axis')
    parser.add_argument('--sy', type=float, default=2,
                        help='Stretching along the y-axis')
    
    # Parameters of FEM
    parser.add_argument('--model', type=int, default=3,
                        help='Constitutive model, 0:Linear, 1:STVK, 2:Co-rotated, 3: Neohookean')
    parser.add_argument('--Y', type=int, default=100,
                        help='Youngs modulus')
    parser.add_argument('--V', type=float, default=0.1,
                        help='Poisson ratio, [-1, 0.5]')
    parser.add_argument('--dt', type=float, default=0.01,
                        help='Time horizon')
    parser.add_argument('--max_iter', type=int, default=200,
                        help='Iteration')
    parser.add_argument('--mass', type=float, default=1.0,
                        help='Mass of the point') 
    parser.add_argument('--gravity', type=bool, default=True,
                        help='Add gravity')  
    parser.add_argument('--fix', type=bool, default=True,
                        help='Fix a point on the wall')  
    parser.add_argument('--implict', default=True,
                        help='Implict Euler Method, for STVK and Neohookean')
    args = parser.parse_args()
    
    args.ref_node_pos=ref_node_pos

    run_fem_simulation_single(args)

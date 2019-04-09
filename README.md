# project
# Jonah Smith, 156981, Spring 2019
# This project provides code and documentation for the final assignment of MECE 5397
# 
# This repository has 4 subfolders, as follows:
# project/doc/ ~ Contains the final report for this project
# project/src/ ~ Contains the source code (in MATLAB) of the program
# project/bin/ ~ Contains executable files to run the program
# project/tests ~ Contains test cases to be ran by the program
#
# The purpose of this progrm is to find a numerical solution to a simple 2D diffusion equation
#
#      d^2u/dx^2 + d^2u/dy^2 = du/dt
#
# Subject to constraints and boundary conditions:
# a_x < x < b_x;     a_y < y < b_y
# u(x,b_y) = f_b(x)
# u(x,a_y) = g_b(x)
# (du/dx | x=b_x) = 0
# u(a_x,y) = g_b(a_x) + [(y-a_y)/(b_y-a_y)]*(f_b(a_x) - g_b(a_x))
#
# Let  the given constants and functions be defined as the following:
# a_x = a_y = 0
# b_x = b_y  2*pi
# f_b(x) = (b_x - x)^2 * cos(pi*x/b_x)
# g_b(x) = x*(b_x - x)^2
#
#
# Two discretizations will be used and their solutions are compared to one another
# Discretization 1: Explicit
# Discretization 2: Crank-Nicolson (Implicit)
#
# Time integration is carried out to steady state, and grid convergence is satisfied
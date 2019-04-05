# The Kinematics of Containment for N-dimensional Ellipsoids
This article introduces the theory of the Kinematic of Containment and studies in details for the case of ellipsoids in n-dimensional Euclidean space. This theory can be applied to error estimations for parts handling tasks and robotic motion planning in challenging environment. The work has been published in the ASME Journal of Mechanisms and Robotics.

Authors: Sipu Ruan (<ruansp@jhu.edu>), Jianzhong Ding, Qianli Ma, Gregory S. Chirikjian

## Introduction
Knowing the set of allowable motions for a convex body moving inside a slightly larger one is of use in applications such as automated assembly mechanisms, robot motion planning, etc. The theory behind this is called the ``Kinematics of Containment (KC)'', and a simplified model of such a convex body, i.e. an ellipsoid, provides clean descriptions and efficient computational processes. This article, in particular, studies a subset of the allowable motions for an n-dimensional ellipsoid being fully contained in another. The problem is addressed in both algebraic and geometric ways, and two lower bounds of the allowable motions are proposed. Containment checking processes for a specific configuration of the moving ellipsoid and the calculations of the volume of the proposed lower bounds in configuration space (C-space) are introduced. Examples for the proposed lower bounds in 2D and 3D Euclidean space are implemented and the corresponding volumes in C-space are compared with different shapes of the ellipsoids. Practical applications using the proposed theories in motion planning problems and parts-handling mechanisms are then discussed.

## Repository Structure
```bash
├── test
│   ├── "main_aspect_ratio_2d.m": Comparisons of KC volume for 2D w.r.t aspect ratio
│   ├── "main_aspect_ratio_3d.m": Comparisons of KC volume for 3D w.r.t aspect ratio
│   ├── "main_infla_2d.m": Comparisons of KC volume for 2D w.r.t inflation factor
│   ├── "main_infla_3d.m": Comparisons of KC volume for 3D w.r.t inflation factor
│   ├── "main_query_time_2d.m": Comparisons of running time for containment checking for 2D case
│   ├── "main_query_time_3d.m": Comparisons of running time for containment checking for 3D case
│   ├── "main_kcc_plots_2d.m": Visualizations of different KC C-space for 2D case
│   ├── "main_query_visual_cvx_2d.m": Querying results of Convex Lower Bound for 2D case
│   ├── "main_query_visual_cvx_3d.m": Querying results of Convex Lower Bound for 3D case
│   ├── "main_query_visual_geo_2d.m": Querying results of Geometric Lower Bound for 2D case
│   └── "main_query_visual_geo_3d.m": Querying results of Geometric Lower Bound for 3D case
├── src
│   ├── cvx_lower_bound
│   │   ├── "cvxLB_2d.m": Main steps of computing Convex Lower Bound for 2D case
│   │   ├── "cvxLB_3d.m": Main steps of computing Convex Lower Bound for 3D case
│   │   ├── "KC_Extreme.m": Convex optimizations for extreme configuration with max magnitude (2D)
│   │   ├── "KC_Extreme_3d.m": Convex optimizations for extreme configuration with max magnitude (3D)
│   │   ├── "maxAngle.m": Computations for max angle that the smaller ellipsoid can rotate
│   │   ├── "configValidation3D.m": Validations of fully contained configurations for 3D case
│   │   ├── "inpoly.m": Point-in-polyhedron test in n-dimensional Euclidean space
│   │   ├── "insimplex.m": Point-in-simplex test in n-dimensional Euclidean space
│   │   ├── "vol_poly.m": Volume computations for n-D polyhedron
│   │   └── "vol_simplex.m": Volume computations for simplex in n-D Euclidean space
│   ├── geo_lower_bound
│   │   ├── "geoLB_2d.m": Main steps of computing Geometric Lower Bound for 2D case
│   │   ├── "geoLB_3d.m": Main steps of computing Geometric Lower Bound for 3D case
│   │   ├── "cvxShpInMink.m": Computations of polyhedron in Minkowski difference boundary
│   └── └── "ptInShrunkPoly.m": Point query inside polyhedron in shrunk space
├── include
│   ├── "mink_2d.m": Main steps of computing "Actual" KC C-space for 2D case
│   ├── "mink_3d.m": Main steps of computing "Actual" KC C-space for 3D case
│   ├── "extremeDistAxis.m": Extreme distance a sphere can move along semi-axis inside an ellipsoid
│   ├── "Ellipse_plot.m": Plot an ellipsoid for either 2D or 3D case
│   └── "rot2.m": A 2D rotation
├── mat
│   ├── "Hhc_rot_2D.mat": Symbolic matrices of H,h,c in approximated algebraic containment condition
│   └── "Hhc_3D.mat": Symbolic matrices of H,h,c in approximated algebraic containment condition
└── LICENSE
```

## Dependencies
* [Ellipsoidal Toolbox (Version 1.2.0.0)](https://www.mathworks.com/matlabcentral/fileexchange/21936-ellipsoidal-toolbox-et): Implementation of the ellipsoidal calculus and ellipsoidal methods for reachability analysis. (By Alex Kurzhanskiy)

* [Robotics Toolbox (Release 10)](http://petercorke.com/wordpress/toolboxes/robotics-toolbox): Functions that are useful for the study and simulation of classical arm-type robotics, for example such things as kinematics, dynamics, and  trajectory generation. (By Peter Corke)



# 1D CutFEM Poisson Solver: Unfitted Dirichlet Boundary Implementation

## Description

This repository contains a self-contained MATLAB implementation of the Cut Finite Element Method (CutFEM) for a one-dimensional Poisson problem. This code demonstrates how to enforce Dirichlet boundary conditions at locations that do not align with the mesh nodes (unfitted boundaries). This approach is highly relevant for industrial simulations involving moving interfaces or complex geometries where traditional remeshing is computationally expensive.

## Mathematical Problem

The solver addresses the steady-state diffusion equation:

-kappa * (d2u / dx2) = f(x)

The domain is defined from x = 0 to x = L.

Boundary Conditions:

1. Strong Dirichlet conditions are applied at the domain ends (x = 0 and x = L).
2. Internal Dirichlet constraints (u(c1) = 1 and u(c2) = 0) are enforced at arbitrary internal points c1 and c2 using the Symmetric Nitsche Method.
3. A source term (f) is active only within the segment between c1 and c2.

## Technical Features

1. Higher-Order Elements: Supports polynomial degrees p=1 (linear), p=2 (quadratic), and p=3 (cubic).
2. Sub-element Integration: When an element is cut by a boundary point, the code identifies the sub-intervals and performs numerical integration on each piece separately to maintain high accuracy.
3. Symmetric Nitsche Method: A consistent variational framework for enforcing Dirichlet conditions on cut interfaces without the ill-conditioning issues of standard penalty methods.
4. Ghost Penalty Stabilization: To prevent numerical instability when a boundary point is extremely close to a mesh node (the small-cut problem), the code includes ghost penalty stabilization (face-jump penalties).

## Configuration Parameters

The following variables can be modified in the Controls section of the script:

* L: Length of the domain.
* Ne: Number of elements in the mesh.
* pDeg: Polynomial degree (1, 2, or 3).
* kappa: Diffusion coefficient.
* fVal: Source term magnitude.
* c0: Initial position of the left boundary marker.
* shift: Distance to translate the internal segment.
* gammaN: Nitsche penalty parameter (typically between 20 and 100).
* useGhost: Set to true to enable ghost penalty stabilization.
* gammaG: Ghost penalty coefficient.

## How to Use

1. Ensure you have MATLAB installed.
2. Copy the code into a file named cutfem_1d_marked_driver.m.
3. Run the script in the MATLAB Command Window.
4. The solver will generate a plot showing the exact solution versus the FEM approximation.

## Simulation Output

The script generates a visualization that includes:

* Vertical dashed lines representing the element edges.
* Individual colored curves for each element to show the local polynomial solution.
* Specific markers (square and diamond) at the internal points c1 and c2 to verify that the boundary conditions are satisfied exactly at the cut locations.

---

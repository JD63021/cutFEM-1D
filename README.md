# 1D CutFEM Poisson Solver: Unfitted Dirichlet Boundary Implementation

## Description

# 1D CutFEM Poisson Solver with Numerical Integration Visualization

## Description

This repository contains a self-contained MATLAB implementation of the Cut Finite Element Method (CutFEM) for solving a 1D Poisson problem. A specific highlight of this version is the inclusion of diagnostic visualizations that illustrate how numerical integration (Gauss quadrature) is performed on elements that are "cut" by internal boundary points.

This code is ideal for educational purposes or for debugging unfitted mesh methods where standard element-based integration fails due to discontinuities or internal constraints.

## Mathematical Problem

The solver addresses the 1D steady-state diffusion equation:

-kappa * (d2u / dx2) = f(x)

The domain is defined on the interval [0, L].

Boundary Conditions:

1. Strong Dirichlet conditions are applied at the global endpoints x = 0 and x = L.
2. Internal Dirichlet constraints (u = 1 at c1 and u = 0 at c2) are enforced using the Symmetric Nitsche Method. These points do not need to align with the mesh nodes.
3. A source term f is applied only to the region between the internal markers c1 and c2.

## Technical Features

1. Higher-Order Discretization: Supports polynomial degrees p=1, p=2, and p=3.
2. Unfitted Boundary Enforcement: Uses the Symmetric Nitsche Method to impose Dirichlet values at arbitrary coordinates within an element.
3. Ghost Penalty Stabilization: Includes face-jump stabilization to maintain system conditioning when an internal boundary is very close to a mesh node.
4. Sub-element Integration: The code detects "cuts" within an element, divides the element into sub-intervals, and maps Gauss quadrature points to these sub-intervals for exact integration of the stiffness matrix and load vector.

## Visualization of Integration (Gauss Boxes)

Unique to this script are the "Gauss Boxes" plots. When an element is cut by a boundary point, the integration must be split. The script generates:

1. Main FEM View: Shows the final solution, element boundaries, and the accuracy of the internal Dirichlet points.
2. Left Subpiece Figure: A schematic representation of 2-point Gauss quadrature boxes mapped to the left portion of the cut element.
3. Right Subpiece Figure: A schematic representation of Gauss quadrature boxes mapped to the right portion of the cut element.

These boxes help verify that the numerical integration correctly spans the physical sub-segments created by the "cut."

## Configuration Parameters

Users can modify the following in the Controls section:

* L: Domain length.
* Ne: Number of elements (set to 3 by default for clarity in visualization).
* pDeg: Polynomial degree (1, 2, or 3).
* kappa: Diffusion coefficient.
* fVal: Magnitude of the source term.
* c0, shift: Controls the positioning of the internal Dirichlet markers.
* gammaN: Nitsche penalty parameter.
* useGhost: Toggle for ghost penalty stabilization.
* gammaG: Stabilization coefficient.

## How to Run

1. Copy the code into a file named cutfem_1d_marked_driver_boxes.m.
2. Execute the script in MATLAB.
3. Three figures will be produced: the main solution plot and two schematic plots showing the integration logic for the first cut element encountered.

---

**Next Step:** Would you like me to prepare a README for a script involving 2D mesh generation in Gmsh or a specific OpenFOAM case setup?

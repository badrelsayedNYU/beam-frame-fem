# FrameLab - Beam and Frame FEM Stepper

An interactive web tool that walks through the direct stiffness method for **2D and 3D beam/frame problems**, one stage at a time. Built for ME-UY 4214 Finite Element Modeling (NYU Tandon), Bonus Assignment 2.

This follows the browser-based format of the earlier [Truss FEM Stepper](https://github.com/badrelsayedNYU/truss-fem), with beam bending, nodal rotations, and 3D torsion added.

## What it shows

| Stage | Content |
|---|---|
| 1. Model & DOFs | Structural drawing, coordinates, supports, loads, and global DOF numbering |
| 2. Element stiffness | Axial, bending, and torsional contributions; shape functions, coefficients, and local stiffness |
| 3. Transform | Local/global axes, rotation R, transformation T, and global element stiffness |
| 4. Assembly | Local-to-global DOF mapping, individual member contribution, complete K, and inspection of any assembly entry |
| 5. Constraints | Prescribed and free DOFs, complete load vector, reduced Kff, and RHS adjustment for nonzero prescribed movements |
| 6. Solve | Diagonal scaling, Gaussian elimination with partial pivoting, intermediate matrices, row operations, and back substitution |
| 7. Results | Nodal translations/rotations, support reactions, residuals, local member end forces, and magnified deformation |

## Built-in problems

- Textbook Example 5.1: 2D portal frame, four nodes and three members (supplied test case 1).
- 3D cantilever with axial force, bending in two planes, and torque; all six tip responses checked against analytical formulas.
- Three-member 3D bent frame demonstrating spatial orientation and assembly.

Users can edit coordinates, properties, loads, supports, and connections, or add/remove nodes and members without editing code. There is no fixed node/member limit; browser resources determine practical size.

## Verification

| Quantity | Tool | Textbook (rounded) |
|---|---|---|
| Node 2 ux (in) | 0.211363 | 0.211 |
| Node 2 uy (in) | 0.00148133 | 0.00148 |
| Node 2 rz (rad) | -0.00152603 | -0.00153 |
| Node 3 ux (in) | 0.209359 | 0.209 |
| Node 3 uy (in) | -0.00148133 | -0.00148 |
| Node 3 rz (rad) | -0.00148600 | -0.00149 |

Differences reflect textbook rounding of stiffness coefficients. Fifteen mathematical validation groups pass, including analytical 3D predictions, equilibrium, symmetry, arbitrary orientation/roll, energy invariance, prescribed movement, deformation interpolation, split-member equivalence, singularity detection, and reproducibility.

Run with Node.js:

```sh
node tests/fem.test.cjs
```

## Running it

Open the GitHub Pages link in the repository's About section, or download `index.html` and open it in a browser. The root `index.html` contains all HTML, CSS, and JavaScript, with **no installation or external dependencies**. All calculations and explanations work offline. Use the built-in User guide for operation instructions. Save/Open model exports and restores JSON.

Readable source modules are in `src/`. The mathematical `element()` function forms stiffness, transformations, and mappings; `solve()` assembles the equations, applies constraints, solves, and recovers reactions; `curve()` computes the deformed centerline. The inline code in `index.html` is the same implementation, bundled so a single download works independently.

## Assumptions and scope

Linear elastic, small-displacement, prismatic Euler-Bernoulli frame elements with rigid joints and nodal loading. 2D uses 3 DOFs/node; 3D uses 6 DOFs/node. Prescribed nonzero translations/rotations are supported. Use consistent units; rotations are radians and member roll is degrees. Dense storage and solving require O(n²) memory and O(n³) work.

Internal hinges/releases, distributed member loads, tapered sections, shear deformation, nonlinear response, buckling, dynamics, and warping torsion are not implemented. The supplied second test PDF contains multiple special beam exercises without a marked target; those special formulations are not claimed as validated cases.

## AI usage disclosure

OpenAI Codex assisted with implementation, explanations, documentation, and validation checks. The operating tool makes no AI, chatbot, LLM, or external service calls; every result is computed locally by fixed mathematical algorithms. See [AI-DISCLOSURE.txt](AI-DISCLOSURE.txt).

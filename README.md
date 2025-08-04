# Crs2D

A 2D coffee ring effect simulator



### Introduction



The coffee ring effect is a fascinating phenomenon in fluid dynamics and materials science! It occurs when a droplet of liquid containing suspended particles (like coffee) evaporates on a surface, leaving behind a ring-shaped stain rather than a uniform deposit.

Here's how it works: When a coffee droplet sits on a surface, it has a curved edge where the liquid is thinner. This edge, called the contact line, becomes "pinned" to the surface due to surface tension and contamination. As the droplet evaporates, the liquid at the thin edges evaporates faster than at the thicker center. To replace this lost liquid, fluid flows outward from the center toward the edges, carrying suspended particles (coffee grounds, pigments, etc.) along with it.

Since the contact line is pinned and can't move inward, all these particles get deposited at the edge as the liquid evaporates, creating the characteristic ring pattern. The center of the original droplet ends up with fewer particles, so it appears lighter.

This effect isn't limited to coffee - it happens with many colloidal suspensions, including ink, paint, and biological samples. Understanding it has important applications in:



Printing and coating technologies

Medical diagnostics (blood spot analysis)

Manufacturing processes where uniform particle deposition is desired

Art restoration and forensics



Researchers have also found ways to suppress the coffee ring effect when uniform deposition is needed, using techniques like adding surfactants, changing particle shapes, or modifying evaporation conditions.



### Mathematical model



\# Coffee Ring Effect: Governing Equations



\## 1. Droplet Height Evolution



The droplet height $h(r,t)$ evolves due to evaporation:



$$\\frac{\\partial h}{\\partial t} = -J(r,t)$$



where $J(r,t)$ is the local evaporation flux. For a spherical cap droplet with contact angle $\\theta$:



$$J(r) = J\_0 \\left(1 - \\frac{r^2}{R^2}\\right)^{-\\lambda}$$



where:

\- $J\_0$ is the evaporation rate at the center

\- $R$ is the droplet radius

\- $\\lambda = (\\pi - 2\\theta)/(2\\pi - 2\\theta)$ is the singularity parameter

\- For $\\theta < 90°$, $\\lambda > 0$ leading to divergent evaporation at the contact line



\## 2. Fluid Velocity Field



From the continuity equation with evaporation source:



$$\\nabla \\cdot \\vec{v} = -\\frac{J(r)}{h(r)}$$



In axisymmetric geometry, assuming thin film approximation:



$$\\frac{1}{r}\\frac{\\partial}{\\partial r}(r v\_r) = -\\frac{J(r)}{h(r)}$$



Integrating from center ($r=0$) to radius $r$:



$$v\_r(r,t) = \\frac{1}{rh(r)} \\int\_0^r J(r') r' dr'$$



This gives the characteristic \*\*outward radial flow\*\* that transports particles to the edge.



\## 3. Particle Transport Equation



The concentration field $c(r,t)$ follows the advection-diffusion equation:



$$\\frac{\\partial c}{\\partial t} + \\nabla \\cdot (c\\vec{v}) = D\\nabla^2 c + S\_c$$



In axisymmetric coordinates:



$$\\frac{\\partial c}{\\partial t} + \\frac{1}{r}\\frac{\\partial}{\\partial r}(r v\_r c) = \\frac{D}{r}\\frac{\\partial}{\\partial r}\\left(r\\frac{\\partial c}{\\partial r}\\right) + S\_c$$



where:

\- $D$ is the particle diffusion coefficient

\- $S\_c$ is a source term accounting for concentration due to evaporation:



$$S\_c = c \\frac{J(r)}{h(r)}$$



\## 4. Particle Deposition



Particles deposit at the substrate when they reach the contact line or when local conditions favor deposition:



$$\\frac{\\partial \\sigma}{\\partial t} = k\_{dep} \\cdot c(r,0) \\cdot f(h,r)$$



where:

\- $\\sigma(r,t)$ is the surface density of deposited particles

\- $k\_{dep}$ is the deposition rate constant

\- $f(h,r)$ is a function that enhances deposition near the contact line



\## 5. Complete Coupled System



The full system of equations is:



\### Height Evolution:

$$\\frac{\\partial h}{\\partial t} = -J\_0 \\left(1 - \\frac{r^2}{R^2}\\right)^{-\\lambda}$$



\### Radial Velocity:

$$v\_r(r,t) = \\frac{J\_0}{rh(r)} \\int\_0^r \\left(1 - \\frac{r'^2}{R^2}\\right)^{-\\lambda} r' dr'$$



\### Particle Concentration:

$$\\frac{\\partial c}{\\partial t} + \\frac{1}{r}\\frac{\\partial}{\\partial r}(r v\_r c) = \\frac{D}{r}\\frac{\\partial}{\\partial r}\\left(r\\frac{\\partial c}{\\partial r}\\right) + c \\frac{J(r)}{h(r)}$$



\### Particle Deposition:

$$\\frac{\\partial \\sigma}{\\partial t} = k\_{dep} c(r,0) \\exp\\left(-\\frac{h(r)}{h\_0}\\right) \\left(1 + \\alpha e^{-\\beta(R-r)}\\right)$$



\## 6. Key Dimensionless Parameters



\### Péclet Number:

$$Pe = \\frac{U L}{D} = \\frac{v\_r R}{D}$$



\- $Pe >> 1$: Advection dominates → Strong ring formation

\- $Pe << 1$: Diffusion dominates → Uniform deposition



\### Capillary Number:

$$Ca = \\frac{\\mu v\_r}{\\sigma\_{LV}}$$



where $\\mu$ is viscosity and $\\sigma\_{LV}$ is surface tension.



\### Contact Angle:

Determines the singularity strength $\\lambda$ and thus the evaporation profile.



\## 7. Boundary Conditions



\### At center ($r = 0$):

\- Symmetry: $\\frac{\\partial c}{\\partial r} = 0$

\- No radial velocity: $v\_r = 0$



\### At contact line ($r = R$):

\- Pinned contact line: $h(R,t) = h\_{edge}$

\- Enhanced deposition

\- Concentration boundary condition depends on deposition model



\## 8. Analytical Solutions



For small times and thin droplets, approximate solutions exist:



\### Velocity near contact line:

$$v\_r \\approx \\frac{J\_0 R}{2h\_{edge}} \\left(\\frac{R-r}{\\delta}\\right)^{-\\lambda}$$



where $\\delta$ is a small cutoff length.



\### Ring width scaling:

The final ring width scales as:

$$w \\sim \\left(\\frac{D t\_{evap}}{Pe}\\right)^{1/2}$$



where $t\_{evap}$ is the total evaporation time.



\## 9. Extensions and Modifications



\### Marangoni Effects:

Add surface tension gradients:

$$\\frac{\\partial c}{\\partial t} + \\nabla \\cdot (c\\vec{v}) = D\\nabla^2 c + \\nabla \\cdot (c\\vec{v}\_{Ma})$$



\### Particle-Particle Interactions:

Modify diffusion coefficient: $D = D\_0(1 + f(c))$



\### Contact Line Dynamics:

Allow contact line motion: $\\frac{dR}{dt} = g(h, \\theta, J)$



\### 3D Effects:

Include vertical velocity component and full 3D transport.



\## 10. Numerical Considerations



\- The evaporation singularity at $r = R$ requires regularization

\- Adaptive mesh refinement near the contact line

\- Implicit time stepping for stability

\- Conservation of total particle number as a check



This mathematical framework captures the essential physics: \*\*evaporation-driven outward flow overcomes diffusion to transport particles to the contact line, where they deposit preferentially due to geometric constraints.\*\*


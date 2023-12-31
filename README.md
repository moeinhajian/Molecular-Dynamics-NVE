# Molecular-Dynamics-NVE
This repository serves as an illustrative example of an NVE (constant Number of particles, constant Volume, and constant Energy) molecular dynamics (MD) calculation.

## Steps:

1. **Generate Atomic Coordinates:**
   - Generate coordinates for atoms within a cubic simulation box with a periodic boundary condition.
   - Simulation box edge length: 12 (all quantities in reduced Lennard-Jones units).
   - Achieve number density of 0.78.

2. **Initialize Atomic Velocities:**
   - Assign initial velocities to atoms corresponding to a reduced temperature of 1.25.

3. **Lennard-Jones Interaction Potential:**
   - Implement the 6-12 Lennard-Jones interaction potential between atoms.
   - Utilize a potential truncated and shifted at a distance of 2.5.

4. **NVE Molecular Dynamics (MD) Simulation:**
   - Conduct NVE MD simulation, maintaining constant number of particles, volume, and energy.
   - Initially rescale velocity for the first few hundred steps to reach a stabilized state.
   - Continue simulation for equilibration.

5. **Analysis and Comparison:**
   - Examine the radial distribution function (RDF) at the beginning and end of the simulation.
   - Compare RDF profiles to assess structural changes and equilibration.

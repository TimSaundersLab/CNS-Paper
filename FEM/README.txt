Code for mechanical analysis of Ventral Neural Cord.
3D FEM analysis from 2D PIV of spin/confocal microscopy analysis:

AUTHOR: Jose J Muñoz, Universitat Politècnica de Catalunya, 2022.

REFERENCE:
K Karkali, P Tiwari, A Singh, S Tlili2, I Jorba, D Navajas, J.J. Muñoz, T.E. Saunders and E Martin-Blanco. Condensation of the Drosophila Nerve Cord is Oscillatory and depends on Coordinated Mechanical Interactions. Developmental Cell, 2022.


MAIN CONTENTS:
1. Main.m:
-Computes 3D stresses and strains, and creates VTK output for Paraview. 
-INPUT: *.mat file in variable "File"='REPOGRIM' and folder "Set.DispFile"
      Each mat contains structure with:
        t0: reference common time-point in minutes for all embryos 
        dt: time-step size in minutes
        xref: reference x in [um] indicating T1 segment of VNC
        pixels: size of each poixel in [um]
        XEul(:,t),YEul(:,t): Eulerian position (x,y) at time-step t of each measure point
        UEul(:,t),VEul(:,t): Velocity at measure point and time-step t in pixel/time-step 
-OUTPUT: *.mat files with kymo data for analysing stress profiles and peaks
-Usage:
a. For stress/strain  kymograph: Run Main.m with Set.Dimensions=2.
   Generates *.mat files that are sused for stress/strain kymographs. 
b. For 3D model.  Run Main.m with Set.Dimensions=3.


2. PlotFromKymo.m:
-Plots stress-strain kymographs, and strain-stress profiles at different times
-Plots evolution of stress peaks on reference kymograph
-INPUT: *.mat files in step 1


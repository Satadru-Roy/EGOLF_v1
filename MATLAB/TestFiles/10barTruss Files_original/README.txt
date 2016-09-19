A Matlab Implementation of the 10 bar truss problem
===================================================

This folder contatins the functions for the 10 bar truss problem, e.g. 
formulated in
Koth,Berke Venkayya, "Comparison  of Optimality Criteria Algorithms
for Minimum Weight Design of Structures", AIAA 78-469, AIAA Journal, Vol
17(2), 1978

The analysis relies on methods from the textbook "A first course in Finite
Elements" from Belytshko and Fish. The code is available at 
http://bcs.wiley.com/he-bcs/Books?action=index&itemId=0470035803&bcsId=3625
and included in this folder. The following file are used from this package:
- assembly.m
- plottruss.m
- postprocessor.m
- preprocessor.m
- solvedr.m
- truss.m
- trusselem.m
- include_flags.m
This m-files represent a general simple truss analysis.

The following files implement the ten bar truss problem:
- objectiveTenBarTruss
   This function calculates the objective value, which is the mass of the
   ten bar truss
- constraintsTenBarTruss
   This function calculates the normalized constraints values. 
    These are the 10 stresses and the max. displacement of a node
- fitnessTenBarTruss
    This function calculates the fitness with a penalty multiplier 
    formulation.
- ten_bar_truss
   This file contains geometric and load information for the path; do not 
   change this file.
- tenBarTrussFunction 
   This function calculates the mass, displacements and
   stresses in the truss elements. This is really only the analysis and not
   function specific. This function is usefull if the optimizer requires 
   e.g. a penalty function (mass is the objective and stress and 
   displacement and constraints)
- calculateWeight
    This function calculates the weight. It is called by the 
    tenBarTrussFunction
- tenBarTruss
    This function calculates the displacement and the stresses in the bar.
    To see its use, refer to the tenBarTrussFunction
- testVenkayya
    optimizes a purley aluminum truss according to AIAA 78-469 (see above);
    It demonstrates the use of the functions in an optimizations scenario.
    This function uses the aluminumConstraints and the aluminumObjective
    functions

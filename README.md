# Epitaxy 2025
Code for exploring epitaxy using EBSD data

Instructions for using MATLAB script “Dyck2025_Epitaxy.m”

Four sections—labelled “USER INPUT #1”, USER INPUT #2, USER INPUT #3, and USER INPUT #4—require editing. USER INPUTS #1–#3 should be filled out prior to initially running the script. INPUT #4 can be filled out prior to running the script if the user already knows the epitaxial relationships. Line numbers are based on the original script submitted in 2025 and may not represent update versions. 

USER INPUTS #1 (Lines 27–33)
    •	‘**PhaseA**’ – the parent phase of interest  
    •	‘**PhaseB**’ – the daughter phase of interest 
    •	‘**Halfwidth**’ – halfwidth (°) used to construct a continuous function that best describes the data points for the daughter phase. This function can then be used to plot contoured density plots for directions and planes of interest in the daughter phase. The defined half-width broadly impacts how smooth the data is. If the half width is too small, then the reconstructed density function is usually oscillating, and individual sample points are visible sharp peaks. If the half width is too large, then the resulting density function is usually too smooth and does not reproduce the features of the original data. The default value is 10, but the function ‘calcKernel’ can be used to optimize the halfwidth if desired.  
    •	‘**gB**’ – misorientation angle used to define a grain boundary
    •	‘**Know_Epitaxy**’ – if user knows specific parent-daughter relationships and would like to test them statistically set to 1 and fill in ‘USER INPUT #4’ prior to running script. If not, leave at default (0). 

USER INPUTS #2 (Lines 34–56)
    •	Import EBSD data. To generate a code to import your EBSD Data type into the command line ‘import_wizard’ and press enter.
    •	Click on the tab that says ‘ebsd’ and click the ‘+’ button on the far right. 
    •	Navigate to and select desired EBSD data, then press open. 
    •	Click through, ensuring that you input the correct ‘Specimen Coordinate System’ and ‘MTEX Plotting Convention’ for your EBSD data/SEM set up. 
    •	A script should pop up with the code required to import the EBSD data, copy and paste this code into the section entitled ‘USER INPUTS #2’.
    
USER INPUTS #3 (Lines 66–73)
    •	‘**h_ADir**’ – directions within the parent phase of interest, replace {h,k,l} with miller or miller-bravais indices of interest.
    •	‘**h_APlan**’ – planes within the parent phase of interest, replace {h,k,l} with miller or miller-bravais indices of interest. 
    •	‘**h_BDir**’ – directions within the daughter phase of interest, replace {h,k,l} with miller or miller-bravais indices of interest. 
    •	‘**h_BPlan**’ – planes within the daughter phase of interest, replace {h,k,l} with miller or miller-bravais indices of interest. 
    
USER INPUTS #4 (Lines 175–185)

    •	‘**PhaseA_Type**’ – state if the epitaxial relationship of interest is a ‘plane’ or a ‘direction’ in the parent grains.
    •	‘**PhaseB_Type**’ – state if the epitaxial relationship of interest is a ‘plane’ or a ‘direction’ in the daughter grains.  
    •	‘**PhaseA_Int**’ – plane/direction of interest in parent phase 
    •	‘**PhaseB_Int**’ – plane/direction of interest in daughter phase 
    •	‘**mis**’ – maximum misorientation between relationships to be considered in statistical analysis. Default is 10°.  
    
    OUTPUT FIGURES
    •	Figure 1: Phase map displaying grain pairs of interest 
    •	Figures 2–5: Pole figures showing the average orientations of parent and daughter grains prior to any rotation
    •	Figures 6–9: Pole figures showing the average orientations of parent and daughter grains after rotation into the parent reference frame (to be used to identify epitaxial/topotaxial relationships).
    
Any issues with data processing and scripts contact Rellie Goddard (rellie.goddard@gmail.com)

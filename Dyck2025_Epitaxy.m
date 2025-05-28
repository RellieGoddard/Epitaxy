%% February 2025 %%
% Written by Rellie M. Goddard (rellie.goddard@gmail.com)

% We present a generalised script for analysing parent-daughter
% relationships using EBSD data. 

% This script first identifies neighbouring parent-daughter pairs. The script then rotates all of 
% the parent grains, that is 'Phase A', by the inverse of their respective Euler angles, such that all Phase A 
% grains are orientated based on the default orientation of the EBSD map,
% that is [c]/[001] aligning with the ‘z’ direction of the fixed co-ordinate system. For each parent-daughter pair, 
% We apply the same rotation to the daughter grain, that is 'Phase B', as we did to the parent grain. 
% Pole figures can then be compared to identify any crystallographic alignments/relationship. 

% For a full description of EBSD analysis see supplementary methods 

% Note there are FOUR USER INPUT sections. USER INPUTS #1-#3 need to be
% completed prior to running the script.
% USER INPUTS #4 is optional and can only be filled out should the user know the epitaxial relationships  

% Note: PLOTTING CONVENTION does not work currently with MTEX version 6.0 beta2

close all; clear all

% Comment out if you prefer images popping up (but there are a lot...)
set(0,'DefaultFigureWindowStyle','docked')

%% USER INPUTS #1 %%
PhaseA = ''; % Parent Phase
PhaseB = ''; % Daughter Phase
halfwidth = ; % In degrees
gB = ; % Grain boundary misorientation, in degrees. 
Know_Epitaxy = [0]; % Set to 1 and fill in USER INPUTS #4 if there are known parent-daughter relationships

%% USER INPUTS #2 %%
% Specify Crystal and Specimen Symmetries %

% Crystal Symmetry
CS = {...};

% plotting convention, set based on SEM/EBSD collection software
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

%% Specify File Names %%

% Path to files
pname = '';

% File name
fname = [pname '\.ctf'];

% Import the Data %
% Create an EBSD variable 
ebsd = EBSD.load(fname,CS,'interface','ctf',...
'convertEuler2SpatialReferenceFrame');

%% CLEAN THE DATA %%
% Edit this in line with your personal preference for cleaning data
BAD_MAD = ebsd.mad > 1;
ebsd(BAD_MAD) = [];
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',gB*degree, 'boundary', 'tight');
ebsd(grains(grains.grainSize<=3)) = [];
ebsd = fill(ebsd);
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',gB*degree, 'boundary', 'tight');

%% USER INPUTS #3 %%
% Directions and planes of interest
% Replace 'h,k,l' with Miller/Miller-Bravais indices of interest
h_ADir = Miller({h,k,l},{h,k,l},ebsd(PhaseA).CS,'direction'); % Directions of interest, Phase A
h_APlan = Miller({h,k,l},{h,k,l},ebsd(PhaseA).CS,'plane'); % Planes of interest, Phase A
h_BDir = Miller({h,k,l},{h,k,l},ebsd(PhaseB).CS,'direction'); % Direction of interest, Phase B
h_BPlan = Miller({h,k,l},{h,k,l},ebsd(PhaseB).CS,'plane'); % Planes of interest, Phase B

%% IDENTIFY ALL GRAINS THAT HAVE THE DESIRED PARENT-DAUGHTER PHASE BOUNDARY %%
grain_id = grains.boundary(PhaseA, PhaseB).grainId;

% Plot grains of interest 
figure
plot(grains(grain_id))

% Find unique boundary pairs
grain_id_unique = unique(grain_id, 'rows');

%% INITAL ORIENTATIONS OF PARENT/DAUGHTER GRAINS %%

% Create a subgroup just with grains that have a parent-daughter (P-D) boundary
% Note this is a n x 2 grain array, with each row containing information on
% one parent and one daughter grain. Should, for example, a parent grain
% have two possible daughters, information on that parent grain will be stored,
% in two places in the sub structure (and vice versa).  

sub = grains(grain_id_unique);

% Separate the two phases, grains are stored in the same order as in ‘sub’
sub_PA_pre = sub(PhaseA);
sub_PB_pre = sub(PhaseB);

Ori_sub_PA_pre = sub_PA_pre.meanOrientation;
Ori_sub_PB_pre = sub_PB_pre.meanOrientation;
Sub_PB_ODF_pre = calcDensity(Ori_sub_PB_pre, 'halfwidth', halfwidth*degree);


%% PLOT ORIGINAL ORIENTATION %%  
% Note these are not true pole figures as they count grains with multiple
% P-D relationship twice.
% NOTE: Change contour levels to fit data 

figure,
plotPDF(Ori_sub_PA_pre, h_ADir, 'anitpodal', 'lower')
hold on;

figure,
plotPDF(Ori_sub_PA_pre, h_APlan, 'anitpodal', 'lower')
hold on;

figure,
plotPDF(Sub_PB_ODF_pre, h_BDir, 'anitpodal', 'lower')
hold on;
plotPDF(Sub_PB_ODF_pre,h_BDir,'antipodal','lower', 'contour',[1:2:10],'linewidth',0.5,'linecolor','k')
mtexColorMap parula
mtexColorbar

figure,
plotPDF(Sub_PB_ODF_pre, h_BPlan, 'anitpodal', 'lower')
hold on;
plotPDF(Sub_PB_ODF_pre,h_BPlan,'antipodal','lower', 'contour',[1:2:10],'linewidth',0.5,'linecolor','k')
mtexColorMap parula
mtexColorbar

%% ROTATE PARENT GRAINS INTO ONE REFERENCE FRAME, APPLY THE SAME ROTATION TO DAUGHTER GRAINS%%
% Replicate the two datasets (to preserve the OG orientation)
sub_PA_post = sub_PA_pre;
sub_PB_post = sub_PB_pre;

for i = 1:1:length(grain_id_unique)

    % Rotate parent by the inverse of its Euler angle 
    sub_PA_post(i).meanOrientation = inv(sub_PA_pre(i).meanOrientation)*sub_PA_pre(i).meanOrientation;
   
    % Apply same rotation to daughter grain
    sub_PB_post(i).meanOrientation = inv(sub_PA_pre(i).meanOrientation)*sub_PB_pre(i).meanOrientation;

end

%% CALULATE MEAN ORIENTATION %% 
Ori_sub_PA_post = sub_PA_post.meanOrientation;
Ori_sub_PB_post = sub_PB_post.meanOrientation;

Sub_PB_ODF = calcDensity(Ori_sub_PB_post, 'halfwidth', halfwidth*degree);

%% PLOT NEW ORIENTATIONS TO IDENTIFY EPITAXY %%
fontsize(8,"points")
figure,
plotPDF(Ori_sub_PA_post, h_ADir, 'anitpodal', 'lower')
hold on;

figure,
plotPDF(Ori_sub_PA_post, h_APlan, 'anitpodal', 'lower')
hold on;

figure,
plotPDF(Sub_PB_ODF, h_BDir, 'anitpodal', 'lower', 'ShowText','off')
hold on;
plotPDF(Sub_PB_ODF,h_BDir,'antipodal','lower', 'contour',[1:2:50],'linewidth',0.5,'linecolor','k', 'ShowText','off')
mtexColorMap parula
mtexColorbar

figure,
plotPDF(Sub_PB_ODF, h_BPlan, 'anitpodal', 'lower', 'ShowText','off')
hold on;
plotPDF(Sub_PB_ODF,h_BPlan,'antipodal','lower', 'contour',[1:2:50],'linewidth',0.5,'linecolor','k', 'ShowText','off')
mtexColorMap parula
mtexColorbar

%% USER INPUTS #4 %%
% STATISTICS OF IDENTIFIED P-D RELATIONSHIPS %
PhaseA_Type = ''; % State if the 'parent' in the P/D relationship of interest is a 'plane' or 'direction'
PhaseB_Type = ''; % State if the 'daughter' in the P/D relationship of interest is a 'plane' or 'direction'

% Specify planes/directions of interest in parent (PhaseA) and daughter
% (PhaseB) in Miller or Miller-Bravis indices 
PhaseA_Int = [h,k,l];
PhaseB_Int = [h,k,l];
mis = 10; % Maximum deviation/misorientation (in degrees) between P/D grains to be considered within error of matching the predefined epitaxial relationship. 

%% CREATE RANDOM DISTRIBUTION FOR COMPARISION %% 
o_PhaseA = orientation.rand(length(grain_id_unique), ebsd(PhaseA).CS);
o_PhaseB = orientation.rand(length(grain_id_unique), ebsd(PhaseB).CS);

%% CALULATE MISORIENTATION & PLOT HISTOGRAM %%

if Know_Epitaxy == 1
    Mis_ori = zeros(length(grain_id_unique), 1);
    Mis_ori_random = zeros(length(grain_id_unique),1);

    for a = 1:1:length(grain_id_unique)
        if length(PhaseA_Int) == 3 & length(PhaseB_Int) == 3
            Vec_A = sub_PA_post(a).meanOrientation*symmetrise(Miller(PhaseA_Int(1), PhaseA_Int(2), PhaseA_Int(3),ebsd(PhaseA).CS,PhaseA_Type));
            Vec_B =  sub_PB_post(a).meanOrientation*symmetrise(Miller(PhaseB_Int(1), PhaseB_Int(2), PhaseB_Int(3), ebsd(PhaseB).CS,PhaseB_Type));
            temp = (180/pi)*angle(unique(Vec_A), unique(Vec_B'));
            temp_2 = min(temp(:));
            Mis_ori(a,1) = temp_2(1);

            Vec_A_rand = o_PhaseA(a)*symmetrise(Miller(PhaseA_Int(1), PhaseA_Int(2), PhaseA_Int(3),ebsd(PhaseA).CS,PhaseA_Type));
            Vec_B_rand =  o_PhaseB(a)*symmetrise(Miller(PhaseB_Int(1), PhaseB_Int(2), PhaseB_Int(3), ebsd(PhaseB).CS,PhaseB_Type));
            temp_rand = (180/pi)*angle(unique(Vec_A_rand), unique(Vec_B_rand'));
            temp_2_rand = min(temp_rand(:));
            Mis_ori_random(a,1) = temp_2_rand(1); 


        elseif length(PhaseA_Int) == 4 & length(PhaseB_Int) == 3
            Vec_A = sub_PA_post(a).meanOrientation*symmetrise(Miller(PhaseA_Int(1), PhaseA_Int(2), PhaseA_Int(3), PhaseA_Int(4), ebsd(PhaseA).CS,PhaseA_Type));
            Vec_B =  sub_PB_post(a).meanOrientation*symmetrise(Miller(PhaseB_Int(1), PhaseB_Int(2), PhaseB_Int(3), ebsd(PhaseB).CS,PhaseB_Type));
            temp = (180/pi)*angle(unique(Vec_A), unique(Vec_B'));
            temp_2 = min(temp(:));
            Mis_ori(a,1) = temp_2(1); 

            Vec_A_rand = o_PhaseA(a)*symmetrise(Miller(PhaseA_Int(1), PhaseA_Int(2), PhaseA_Int(3),PhaseA_Int(4), ebsd(PhaseA).CS,PhaseA_Type));
            Vec_B_rand =  o_PhaseB(a)*symmetrise(Miller(PhaseB_Int(1), PhaseB_Int(2), PhaseB_Int(3), ebsd(PhaseB).CS,PhaseB_Type));
            temp_rand = (180/pi)*angle(unique(Vec_A_rand), unique(Vec_B_rand'));
            temp_2_rand = min(temp_rand(:));
            Mis_ori_random(a,1) = temp_2_rand(1); 


        elseif length(PhaseA_Int) == 3 & length(PhaseB_Int) == 4
            Vec_A = sub_PA_post(a).meanOrientation*symmetrise(Miller(PhaseA_Int(1), PhaseA_Int(2), PhaseA_Int(3),ebsd(PhaseA).CS,PhaseA_Type));
            Vec_B =  sub_PB_post(a).meanOrientation*symmetrise(Miller(PhaseB_Int(1), PhaseB_Int(2), PhaseB_Int(3), PhaseB_Int(4), ebsd(PhaseB).CS,PhaseB_Type));
            temp = (180/pi)*angle(unique(Vec_A), unique(Vec_B'));
            temp_2 = min(temp(:));
            Mis_ori(a,1) = temp_2(1);

            Vec_A_rand = o_PhaseA(a)*symmetrise(Miller(PhaseA_Int(1), PhaseA_Int(2), PhaseA_Int(3),ebsd(PhaseA).CS,PhaseA_Type));
            Vec_B_rand =  o_PhaseB(a)*symmetrise(Miller(PhaseB_Int(1), PhaseB_Int(2), PhaseB_Int(3),PhaseB_Int(4), ebsd(PhaseB).CS,PhaseB_Type));
            temp_rand = (180/pi)*angle(unique(Vec_A_rand), unique(Vec_B_rand'));
            temp_2_rand = min(temp_rand(:));
            Mis_ori_random(a,1) = temp_2_rand(1); 

        elseif length(PhaseA_Int) == 4 & length(PhaseB_Int) == 4
            Vec_A = sub_PA_post(a).meanOrientation*symmetrise(Miller(PhaseA_Int(1), PhaseA_Int(2), PhaseA_Int(3),PhaseA_Int(4), ebsd(PhaseA).CS,PhaseA_Type));
            Vec_B =  sub_PB_post(a).meanOrientation*symmetrise(Miller(PhaseB_Int(1), PhaseB_Int(2), PhaseB_Int(3), PhaseB_Int(4), ebsd(PhaseB).CS,PhaseB_Type));
            temp = (180/pi)*angle(unique(Vec_A), unique(Vec_B'));
            temp_2 = min(temp(:));
            Mis_ori(a,1) = temp_2(1); 

            Vec_A_rand = o_PhaseA(a)*symmetrise(Miller(PhaseA_Int(1), PhaseA_Int(2), PhaseA_Int(3), PhaseA_Int(4),ebsd(PhaseA).CS,PhaseA_Type));
            Vec_B_rand =  o_PhaseB(a)*symmetrise(Miller(PhaseB_Int(1), PhaseB_Int(2), PhaseB_Int(3), PhaseB_Int(4),ebsd(PhaseB).CS,PhaseB_Type));
            temp_rand = (180/pi)*angle(unique(Vec_A_rand), unique(Vec_B_rand'));
            temp_2_rand = min(temp_rand(:));

            Mis_ori_random(a,1) = temp_2_rand(1); 

        end
    end

    % PERCENTAGE OF DAUGHTERS WITH AT LEAST ONE EPITAXYAL BOUNDARY %
    PhaseB_grains_ID = sub_PB_post.id;
    PhaseA_grains_ID = sub_PA_post.id;

    % Find the number of unique daughter grains
    unique_PhaseB = unique(PhaseB_grains_ID);
    unique_PhaseA = unique(PhaseA_grains_ID);

    PhaseB_IDs_Within = [];

    for i = 1:1:length(Mis_ori)
        if Mis_ori(i)<=mis
            PhaseB_IDs_Within = [PhaseB_IDs_Within, PhaseB_grains_ID(i)];
        end
    end


    % PRINT PERCENTAGE OF DAUGHTERS WITH SPECIFIED DEGREES OF OFFSET %
    DAUGHTER_Percentage = round((length(unique(PhaseB_IDs_Within))/length(unique_PhaseB))*100,0);
    BOUNDARY_Percentage = round((length(find(Mis_ori<=mis))/length(Mis_ori))*100,0);
    RAND_Percentage = round((length(find(Mis_ori_random<=mis))/length(Mis_ori_random))*100,0);
end


%% PRINT PERCENTAGE OF DAUGHTERS WITH SPECIFIED DEGREES OF OFFSET %% 

display(DAUGHTER_Percentage)
display(BOUNDARY_Percentage)
display(RAND_Percentage)




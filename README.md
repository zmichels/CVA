# CVA and PGA

Code for Crystallographic Vorticity Axis analysis and Principal Geodesic Analysis


The provided code snippets and functions can be run in MATLAB and are written to
leverage functionality from the open source toolbox, MTEX, v5.X+ and (which
provides tools for importing EBSD data and constructing GrainSet objects). Some
details about the CVA method (and appropriate references) are provided in
comments within the functions.


Please contact Z.Michels (zachary.michels@gmail.com) with any questions or for consultation and feel free to modify these functions as allowed by the MIT open source license guidelines.
Additionally, older versions of the code are available for use with earlier versions of MTEX upon request; however, these have been removed from the repository and are no
longer being maintained. Occasionally users may find that updating to new MTEX versions results in an error with functions here -- this can happen when the MTEX functions change in new versions. In such a case, please return to this repository to check for an update that will work with latest MTEX (I check compatibility with every MTEX release, and update the scripts, if needed, immediately).



REFERENCES:
Zachary D. Michels, Seth C. Kruckenberg, Joshua R. Davis, and Basil
Tikoff, Determining vorticity axes from grain-scale dispersion of
crystallographic orientations, Geology, G36868.1, first published on July 17,
2015, doi:10.1130/G36868.1



# Example application or workflow


After you have imported an EBSD dataset and constructed grains, perhaps try the following workflow example as a starting point.

1. Isolate a set of grains that are appropriate for CVA analysis by creating an index of grains that match conditions suitable for CVA analysis.

        a. grainSize >= 3  --  grains with 3 or more solutions -- grains(grains.grainSize>=3)
    
        b. GOS >= 1-degree --  grains with an orientation spread of 1 degree or more -- grains(grains.GOS>=1*degree)   
    
        c. indexed grains  --  only analyze grains that contain orientation solutions of an indexed phase (grains('indexed') or grains(grains.phase>0)



Example:

condition = grains(grains.grainSize>=3 & grains.GOS>=1*degree & grains.phase>0);

grains = grains(condition);


%% compute orientation dispersion tensor and CVA results for each grain

[gCVA, bulk] = grainsCVA(grains,ebsd)


% NOTE: When you run the analysis, the Matlab terminal will update with "percent done" message to help gage the time left. 



%% Exploring the CVA results

% the output from the CVA analysis will include a grain set that has been appended with CVA-related data. Previous versions of the scripts output the data separately, but this new approach provides for easier workflow, in my opinion.



% Filter results to exclude axes computed from grains with [potentially] insignificant intragranular distortion / lattice curvature.
% perhaps 1-degree of grain-orientation-spread (GOS) is a conservative bare minimum requirement for a grain/kernel to be considered.

gCVA = gCVA(gCVA.GOS>=1*degree)




%% To visualize all of the rotation axes from all of the grains:

figure,

plot(gCVA.CVA,'antipodal','lower','smooth','halfwidth',10*degree)

mtexColorbar('Title','M.U.D.')

% add the "bulk" vorticity axis -- highest MUD using a KDE with 10-degree halfwidth

hold on

plot(bulk,'antipodal','lower','Marker','^','MarkerSize',15,'MarkerEdgeColor','w','MarkerFaceColor','k')




%% Visualize results for grains of a specific phase (ex: 'forsterite')

figure,

plot(gCVA('forsterite').CVA,'antipodal','lower','smooth','halfwidth',10*degree)

mtexColorbar('Title','M.U.D.')






%% NEW: plot mean dispersion tensor

meanT = mean(gCVA.ODT);


figure,

plot(meanT,'antipodal','lower')




# Contact and support

Please feel free to write to Zach Michels zachary.michels@gmail.com) with any questions or concerns, or for any assistance.

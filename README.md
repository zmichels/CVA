# CVA
Code for Crystallographic Vorticity Axis analysis

The provided code snippets and functions can be run in MATLAB and are written to
leverage functionality from the open source toolbox, MTEX, v5.X+ and (which
provides tools for importing EBSD data and constructing GrainSet objects). Some
details about the CVA method (and appropriate references) are provided in
comments within the functions.

Please contact Z.Michels (mich0201@umn.edu) with any questions or for consultation and feel free to modify these functions as allowed by the MIT open source license guidelines.
Additionally, older versions of the code are available for use with earlier versions of MTEX upon request; however, these have been removed from the repository and are no
longer being maintained. Occasionally users may find that updating to new MTEX versions results in an error with functions here -- this can happen when the MTEX functions change in new versions. In such a case, please return to this repository to check for an update that will work with latest MTEX (I check compatibility with every MTEX release, and update the scripts, if needed, immediately).


REFERENCES:
Zachary D. Michels, Seth C. Kruckenberg, Joshua R. Davis, and Basil
Tikoff, Determining vorticity axes from grain-scale dispersion of
crystallographic orientations, Geology, G36868.1, first published on July 17,
2015, doi:10.1130/G36868.1



# Example application or workflow
After you have imported an ebsd dataset and constructed grains, perhaps try the following workflow example as a starting point.

%% Isolate a set of grains that are appropriate for CVA analysis
% create an index of grains that match parameters that will allow for CVA
% analysis:
% grainSize >= 3    --  grains with 3 or more solutions
% GOS               --  grain orientation spread of 1 degree or more
% phase > 0         --  only grains with indexed data (not-indexed = 0).

condition = grains.grainSize>=3&grains.GOS>=1*degree & grains.phase>0;

% select the grains
gCVA = grains(condition);

%% compute orientation dispersion tensor for each grain
[gCVA,bulk] = calcGrainsDispersion(gCVA,ebsd(gCVA));

% NOTE: When you run the analysis, the Matlab terminal will update with "percent done" message to help gage the time left. 

%% Exploring the CVA results
% the output from the CVA analysis will include a grain set that has been appended with CVA-related data. Previous versions of the scripts output the data separately, but this new approach provides for easier workflow, in my opinion.

% to visualize all of the rotation axes from all of the grains:
figure,
plot(gCVA.CVA,'antipodal','lower','smooth','halfwidth',10*degree)
mtexColorbar('Title','M.U.D.')

% add the "bulk" vorticity axis -- highest MUD using a KDE with 10-degree halfwidth
hold on
plot(bulk,'antipodal','lower','Marker','^','MarkerSize',15,'MarkerEdgeColor','w','MarkerFaceColor','k')


% visualize results for grains of a specific phase (ex: 'forsterite')
figure,
plot(gCVA('forsterite').CVA,'antipodal','lower','smooth','halfwidth',10*degree)
mtexColorbar('Title','M.U.D.')


# Weighting or filtering CVA results??
Some folks have noted that when applying CVA to samples which do not exhibit much intragranular distortion of the lattice (and hence not much orientation dispersion of the lattice), the results of the CVA analysis are "noisy" due to the fact that the CVA analysis will assign a best-fit rotation for any amount of dispersion in the orientations (even if the orientation spread is potentially insignificantly small). One proposed method for teasing out a more significant signal from the results may be to filter or weight the results prior to visualizing or further analyzing. A colleague suggested using a ratio of the axes magnitudes of the [best-fit] dispersion tensor to ignore results for which the ratio is "not significant". Determining that value is arbitrary. In the following example snippet, I use the value of 1.4 as a threshold for the ratio between the magnitudes of the primary and secondary principal axes of each grain-scale dispersion tensor. I borrowed this value from my experience analyzing shape-preferred orientations of grains, in which I typically ignore grain-shape ellipses that have an aspect ratio less than 1.4. By imposing a lower limit on the CVA results, you exclude results that do not have a strongly "preferred" best-fit to their orientation dispersion. Increasing the lower limit, and/or plotting subsets of the data within different ratio ranges, may be one interesting way to explore details of CVA results. For many datasets, using a 1.4 ratio cutoff does not seem to filter out many results and does not seem to have a perceptible influence on the plot.  However, for some samples, it has a *dramatic* effect on the plot, and it may be worth exploring such results in more detail.


% Plot a filtered version of the CVA data
figure,
plot(gCVA(gCVA.mag1./gCVA.mag2>=1.4).CVA,'antipodal','lower','smooth','halfwidth',10*degree)
mtexColorbar('Title','M.U.D.')



# Contact and support
Please feel free to write to Zach Michels (mich0201@umn.edu) with any questions or concerns, or for any assistance.











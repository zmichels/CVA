function [rotCVA] = rotToCVAref(gCVA)
%%
% dispersion axes as orientation
cvaOr = (orientation.map(vector3d.X,gCVA.CVA,vector3d.Z,gCVA.eV3,crystalSymmetry('mmm')));

% ODF of dispersion axes/orientations
cvaODF = calcDensity(cvaOr,'halfwidth',10*degree);

% mode of CVA
cvaMode = calcModes(cvaODF);

%% compute rotation
rotCVA = rotation.map(cvaMode*vector3d.X,vector3d.Z,cvaMode*vector3d.Y,vector3d.Y);
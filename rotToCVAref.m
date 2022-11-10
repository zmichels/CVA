function [rotCVA] = rotToCVAref(gCVA)
%% mean orientation dispersion tensor
mODT = mean(gCVA.ODT);

[v,~] = eig(mODT);


% % dispersion axes as orientation
% cvaOr = (orientation.map(vector3d.X,gCVA.CVA,vector3d.Z,gCVA.eV3,crystalSymmetry('mmm')));
% 
% % ODF of dispersion axes/orientations
% cvaODF = calcDensity(cvaOr,'halfwidth',10*degree);
% 
% % mode of CVA
% cvaMode = calcModes(cvaODF);

%% compute rotation
rotCVA = rotation.map(v(3),vector3d.Z,v(1),vector3d.Y);
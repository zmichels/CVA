function [xtalCVA] = invCVA(gPhase)

%% note only one phase permitted
meanO = orientation(gPhase.meanRotation,gPhase.CS);

%% xtalCVA (inversion of CVA with mean orientation)
xtalCVA.CVA = meanO.\gPhase.CVA(:);
xtalCVA.eV1 = meanO.\gPhase.eV1(:);
xtalCVA.eV2 = meanO.\gPhase.eV2(:);
xtalCVA.eV3 = meanO.\gPhase.eV3(:);
xtalCVA.ODT = meanO.\gPhase.ODT(:);


%% plot it
% setup new mtex plot
mtexFig = mtexFigure('position',[0 0 1000 1000]);
titleOpt = {'visible','on','color','k'};

% firt plot is inverted CVA (dominant dispersion/rotation axis)
plot(xtalCVA.CVA,'antipodal','lower','smooth','halfwidth',10*degree,'fundamentalRegion')
mtexTitle('invCVA (eV1)',titleOpt{:})

% second plot is inverted intermediate eigenvector of the dispersion tensor
nextAxis
plot(xtalCVA.eV2,'antipodal','lower','smooth','halfwidth',10*degree,'fundamentalRegion')
mtexTitle('inv eV2',titleOpt{:})

% third plot is inverted axes of smallest eignevcevtor of disperision tensor
nextAxis
plot(xtalCVA.eV3,'antipodal','lower','smooth','halfwidth',10*degree,'fundamentalRegion')
mtexTitle('inv eV3',titleOpt{:})





end
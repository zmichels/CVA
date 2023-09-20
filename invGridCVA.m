function [xtalGridCVA] = invGridCVA(eCVAPhase)


%% note only one phase permitted
o = orientation(eCVAPhase.orientations,eCVAPhase.CS);

%% xtalCVA (inversion of CVA with mean orientation)
xtalGridCVA.CVA = o.\eCVAPhase.CVA(:);
xtalGridCVA.eV1 = o.\eCVAPhase.eV1(:);
xtalGridCVA.eV2 = o.\eCVAPhase.eV2(:);
xtalGridCVA.eV3 = o.\eCVAPhase.eV3(:);
xtalGridCVA.ODT = o.\eCVAPhase.ODT(:);


%% plot it

% Fourth plot is inverted mean ODT
nextAxis
plot(mean(xtalGridCVA.ODT),'antipodal','lower','smooth','halfwidth',10*degree,'fundamentalRegion')
mtexTitle('inv ODT')

cond = eCVAPhase.mag1 > quantile(eCVAPhase.mag1,.95);
nextAxis
plot(xtalGridCVA.CVA(cond),'antipodal','lower','smooth','halfwidth',10*degree,'fundamentalRegion')
mtexTitle('inv CVA')



end
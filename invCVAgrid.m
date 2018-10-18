%% inverse calculation to get CVA vectors in xtal reference frame
xtalGrainCVA = vector3d();
eCVAphase = eCVA(phase);
xtalGraineV2 = vector3d();
xtalGrainCVAor = orientation();
xtalGraineV3 = vector3d();


for i = 1:length(eCVAphase)
    xtalGraineV1(i) = inv(eCVAphase(i).orientations) * eCVAphase.CVA(i);
    xtalGraineV2(i) = inv(eCVAphase(i).orientations) * eCVAphase.eV2(i);
    xtalGraineV3(i) = inv(eCVAphase(i).orientations) * eCVAphase.eV3(i);
end


%% symmetrize
[xtalSymeV1,~,~] = symmetrise(xtalGraineV1,eCVAphase.CS);
[xtalSymeV3,~,~] = symmetrise(xtalGraineV3,eCVAphase.CS);
[xtalSymeV2,~,~] = symmetrise(xtalGraineV2,eCVAphase.CS);


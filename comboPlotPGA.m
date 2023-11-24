function [] = comboPlotPGA(gCVA)

%
figure
m = newMtexFigure;
plot(gCVA.eV1,'antipodal','lower','smooth','halfwidth',10*degree)
mtexTitle('eV1')


nextAxis
plot(gCVA.eV2,'antipodal','lower','smooth','halfwidth',10*degree)
mtexTitle('eV2')

nextAxis
plot(gCVA.eV3,'antipodal','lower','smooth','halfwidth',10*degree)
mtexTitle('eV3')

setColorRange('equal');
cb = mtexColorbar('Title','M.U.D.');

if cb.Limits(2)>2
    setColorRange([1 floor(cb.Limits(2))])
else 
    setColorRange([1 ceil(cb.Limits(2))])
end

m.parent.Name = 'Eigenvectors of Dispersion Tensors';
hold off
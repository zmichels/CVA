warning off % because figure saving/exporting gives warnings lately

%% Check for variables used in the script

% if no phase is specified, assume forsterite
if exist('phase') == 0
    phase = 'fo';
end


% if no sample name is designated, use the name 'sample'
if exist('sampleName') == 0
    sampleName = 'sample';
end



%% ******************************************
%*******************************************%
%    CVA axes to crystal reference frame    %
%*******************************************%
%% inverse calculation to get CVA vectors in xtal reference frame
xtalGrainCVA = vector3d();
gCVAphase = gCVA(phase);
xtalGraineV2 = vector3d();
xtalGrainCVAor = orientation();
xtalGraineV3 = vector3d();
o = gCVAphase.meanOrientation(:);
cva = gCVAphase.CVA(:);
eV2 = gCVAphase.eV2(:);
eV3 = gCVAphase.eV3(:);

xtalGraineV1 = o.\cva;
xtalGraineV2 = o.\eV2;
xtalGraineV3 = o.\eV3;

D=orientation('map',xvector,xtalGraineV1,yvector,xtalGraineV2);


%% symmetrize
[xtalSym] = symmetrise(D,gCVAphase.CS);

%% plot
cs = xtalSym.CS;
h = [Miller(1,0,0,'direction',cs),Miller(0,1,0,'direction',cs),Miller(0,0,1,'direction',cs)];
f = figure;
plot(vector3d(xtalSym*xvector),'antipodal','lower','smooth','halfwidth',10*degree)
annotate(h(3),'label',{'[001]'},'BackgroundColor','w')
annotate(h(1),'label',{'[100]'},'BackgroundColor','w')
annotate(h(2),'label',{'[010]'},'BackgroundColor','w')
f.Position = [0,0,400,400];
pos = f.Position;
mtexColorbar('title','R1 axes: M.U.D.')
saveas(gcf,sprintf('%s R1 axes.png',sampleName));

f = figure;
plot(vector3d(xtalSym*yvector),'antipodal','lower','smooth','halfwidth',10*degree)
annotate(zvector,'label',{'[001]'},'BackgroundColor','w')
annotate(xvector,'label',{'[100]'},'BackgroundColor','w')
annotate(yvector,'label',{'[010]'},'BackgroundColor','w')
f.Position = [0,0,400,400];
pos = f.Position;
mtexColorbar('title','R2 axes: M.U.D.')
saveas(gcf,sprintf('%s R2 axes.png',sampleName));

f = figure;
plot(vector3d(xtalSym*zvector),'antipodal','lower','smooth','halfwidth',10*degree)
annotate(zvector,'label',{'[001]'},'BackgroundColor','w')
annotate(xvector,'label',{'[100]'},'BackgroundColor','w')
annotate(yvector,'label',{'[010]'},'BackgroundColor','w')
f.Position = [0,0,400,400];
pos = f.Position;
mtexColorbar('title','R3 axes:  M.U.D.')
saveas(gcf,sprintf('%s R3 axes.png',sampleName));


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
xtalCVA = vector3d();
eCVAphase = eCVA(phase);
xtaleV2 = vector3d();
xtalGrainCVAor = orientation();
xtaleV3 = vector3d();
o = eCVAphase.orientations(:);
cva = eCVAphase.CVA(:);
eV2 = eCVAphase.eV2(:);
eV3 = eCVAphase.eV3(:);

xtaleV1 = o.\cva;
xtaleV2 = o.\eV2;
xtaleV3 = o.\eV3;

D=orientation('map',xvector,xtaleV1,yvector,xtaleV2);


% 
% for i = 1:length(eCVAphase)
%     xtaleV1(i) = inv(eCVAphase(i).orientations) * eCVAphase.CVA(i);
%     xtaleV2(i) = inv(eCVAphase(i).orientations) * eCVAphase.eV2(i);
%     xtaleV3(i) = inv(eCVAphase(i).orientations) * eCVAphase.eV3(i);
% end


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

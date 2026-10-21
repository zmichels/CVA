function tests=testGrainsCVAMtex61
tests=functiontests(localfunctions);
end
function testThreePhaseLoopAndBatched(t)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
new=fullfile(root,'matlab','mtex61');
t.assertTrue(isfile(fullfile(new,'grainsCVA.m')),'Maintained grainsCVA must exist.');
old=pwd; cleanup=onCleanup(@()cd(old)); cd(tempdir); addpath(new,'-begin'); clear grainsCVA PGA;
ebsd=mtexdata('forsterite');
[grains,ebsd.grainId]=calcGrains(ebsd('indexed'),'angle',10*degree);
g=grains(grains.numPixel>=3 & grains.GOS>0.5*degree);
% Bound the fixture while retaining every indexed phase.
ids=[];
for p=unique(g.phaseId).'
 gp=g(g.phaseId==p); ids=[ids;gp.id(1:min(8,length(gp)))];
end
g=g('id',ids); ebsd=ebsd(g);
[a,~,da]=grainsCVA(g,ebsd,'method','loop','computeBulkVector',false,'storeTensor',true);
[b,~,db]=grainsCVA(g,ebsd,'method','batched','computeBulkVector',false,'storeTensor',true);
t.verifyEqual(a.id,b.id); t.verifyEqual(a.phaseId,b.phaseId);
t.verifyEqual(a.id,sort(ids));
t.verifyEqual(numel(unique(a.phaseId)),3);
maxTensor=0; maxEigen=0; maxAxis=0;
for i=1:length(a)
 ta=matrix(a.ODT(i)); tb=matrix(b.ODT(i));
 maxTensor=max(maxTensor,norm(ta-tb,'fro')/max(norm(ta,'fro'),1e-20));
 ma=[a.mag1(i),a.mag2(i),a.mag3(i)]; mb=[b.mag1(i),b.mag2(i),b.mag3(i)];
 maxEigen=max(maxEigen,max(abs(ma-mb)./max(abs(ma),1e-20)));
 va=matrix([a.eV1(i),a.eV2(i),a.eV3(i)]); vb=matrix([b.eV1(i),b.eV2(i),b.eV3(i)]);
 maxAxis=max(maxAxis,max(acos(min(1,abs(sum(va.*vb,1))))/degree));
end
t.verifyLessThan(maxTensor,1e-4); t.verifyLessThan(maxEigen,1e-4); t.verifyLessThan(maxAxis,0.01);
t.verifyEqual(da.meanMatrices,db.meanMatrices,'AbsTol',1e-12);
% The unmodified public grain path is the independent selection/KDE oracle.
addpath(root,'-begin'); clear PGA grainsCVA;
[legacy,legacyBulk]=grainsCVA(g,ebsd);
addpath(new,'-begin'); clear PGA grainsCVA;
[maintained,maintainedBulk]=grainsCVA(g,ebsd,'method','batched');
t.verifyEqual(maintained.id,legacy.id); t.verifyEqual(maintained.phaseId,legacy.phaseId);
t.verifyEqual(maintained.mag1,legacy.mag1,'AbsTol',1e-12);
t.verifyEqual(maintained.mag2,legacy.mag2,'AbsTol',1e-12);
t.verifyEqual(maintained.mag3,legacy.mag3,'AbsTol',1e-12);
t.verifyEqual(matrix(maintainedBulk),matrix(legacyBulk),'AbsTol',1e-12);
for j=1:length(legacy)
 t.verifyEqual(matrix(maintained.ODT(j)),matrix(legacy.ODT(j)),'AbsTol',1e-12);
end

fprintf('GRAIN_PARITY n=%d phases=%d tensor_rel=%.12g eigen_rel=%.12g axis_deg=%.12g\n',length(a),numel(unique(a.phaseId)),maxTensor,maxEigen,maxAxis);
for repeat=1:3
 tic; grainsCVA(g,ebsd,'method','loop','computeBulkVector',false,'storeTensor',false); loop=toc;
 tic; grainsCVA(g,ebsd,'method','batched','computeBulkVector',false,'storeTensor',false); batch=toc;
 fprintf('GRAIN_TIMING repeat=%d loop=%.6f batched=%.6f\n',repeat,loop,batch);
end
end

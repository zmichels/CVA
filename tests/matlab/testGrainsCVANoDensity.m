function tests = testGrainsCVANoDensity
%TESTGRAINSCVANODENSITY Clean-exit three-phase batch/loop regression.
tests = functiontests(localfunctions);
end

function testThreePhaseExactMeansAndTensors(t)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
maintained = fullfile(root,'matlab','mtex61');
oldpath = path; oldpwd = pwd;
cleanup = onCleanup(@() restoreState(oldpath,oldpwd));
cd(tempdir); addpath(maintained,'-begin'); clear grainsCVA PGA
t.assertEqual(which('grainsCVA'),fullfile(maintained,'grainsCVA.m'));

ebsd = mtexdata('forsterite');
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',10*degree);
eligible = grains(grains.numPixel>=3 & grains.GOS>0.5*degree);
ids = [];
for phase = unique(eligible.phaseId).'
    inPhase = eligible(eligible.phaseId==phase);
    ids = [ids;inPhase.id(1:min(8,length(inPhase)))]; %#ok<AGROW>
end
selected = eligible('id',ids); ebsd = ebsd(selected);
[loop,~,loopInfo] = grainsCVA(selected,ebsd,'method','loop','computeBulkVector',false);
[batch,~,batchInfo] = grainsCVA(selected,ebsd,'method','batched','computeBulkVector',false);

t.verifyEqual(length(batch),24);
t.verifyEqual(numel(unique(batch.phaseId)),3);
t.verifyEqual(batch.id,loop.id);
t.verifyEqual(batch.phaseId,loop.phaseId);
t.verifyEqual(batchInfo.meanMatrices,loopInfo.meanMatrices,'AbsTol',1e-12);
for i=1:length(batch)
    A=matrix(loop.ODT(i)); B=matrix(batch.ODT(i));
    t.verifyLessThan(norm(A-B,'fro')/max(norm(A,'fro'),1e-20),1e-12);
    t.verifyEqual([batch.mag1(i),batch.mag2(i),batch.mag3(i)], ...
        [loop.mag1(i),loop.mag2(i),loop.mag3(i)],'AbsTol',1e-12);
    t.verifyLessThan(acos(min(1,abs(dot(batch.eV1(i),loop.eV1(i))))),0.01*degree);
end
end

function restoreState(oldpath,oldpwd)
path(oldpath); cd(oldpwd);
end

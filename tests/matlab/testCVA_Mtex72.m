function tests = testCVA_Mtex72
%TESTCVA_MTEX72 Bounded MTEX 7.2 compatibility and phase-aware CVA checks.
tests = functiontests(localfunctions);
end

function setupOnce(t)
t.TestData.previousFolder = pwd;
t.TestData.root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
t.TestData.maintained = fullfile(t.TestData.root,'matlab','mtex61');
t.TestData.pathAdded = false;
t.assertTrue(startsWith(getMTEXpref('version'),'mtex-7.2'), ...
    'This suite requires MATLAB MTEX 7.2.');
cd(tempdir);
addpath(t.TestData.maintained,'-begin');
t.TestData.pathAdded = true;
clear PGA grainsCVA gridCVA
t.assertEqual(which('PGA'),fullfile(t.TestData.maintained,'PGA.m'));
end

function teardownOnce(t)
if t.TestData.pathAdded
    cd(t.TestData.previousFolder);
    rmpath(t.TestData.maintained);
end
end

function testPGAAgainstMtex61Fixtures(t)
raw = readmatrix(fullfile(t.TestData.root,'tests','fixtures','pga_rotations.csv'));
fixture = jsondecode(fileread(fullfile(t.TestData.root,'tests','fixtures','pga_expected.json')));
for k = 1:numel(fixture.cases)
    c = fixture.cases(k);
    rows = raw(raw(:,1) == k,2:10);
    matrices = zeros(3,3,size(rows,1));
    for j = 1:size(rows,1)
        matrices(:,:,j) = reshape(rows(j,:),3,3).';
    end
    ori = orientation.byMatrix(matrices,crystalSymmetry(c.crystal_symmetry),specimenSymmetry('1'));
    [~,m,T] = PGA(ori);
    [~,mExplicit,TExplicit] = PGA(ori,mean(ori));
    t.verifyLessThan(norm(matrix(T)-c.tensor,'fro'),1e-10);
    t.verifyLessThan(max(abs(m(:)-c.eigenvalues(:))),1e-10);
    t.verifyLessThan(norm(matrix(T)-matrix(TExplicit),'fro'),1e-12);
    t.verifyLessThan(max(abs(m(:)-mExplicit(:))),1e-12);
end
end

function testThreePhaseGrainsAndGrid(t)
[x,y] = meshgrid(0:8,0:8);
n = numel(x);
symmetries = {'notIndexed',crystalSymmetry('mmm'), ...
    crystalSymmetry('m-3m'),crystalSymmetry('222')};
phase = 2*ones(n,1);
phase(x(:)>=3) = 3;
phase(x(:)>=6) = 4;
rotations = rotation.byEuler((10+x(:)*.4+y(:)*.1)*degree, ...
    (20+y(:)*.3)*degree,(30+x(:).*y(:)*.07)*degree);
ebsd = EBSD(vector3d(x(:),y(:),zeros(n,1)),rotations,phase, ...
    symmetries,struct(),'phaseMap',[0 1 2 3]);
[grains,ebsd] = calcGrains(ebsd('indexed'),'angle',10*degree);
[loop,~] = grainsCVA(grains,ebsd,'method','loop','computeBulkVector',false);
[batch,~] = grainsCVA(grains,ebsd,'method','batched','computeBulkVector',false);
t.verifyEqual(length(batch),3);
t.verifyEqual(numel(unique(batch.phaseId)),3);
t.verifyEqual(batch.id,loop.id);
for k = 1:length(batch)
    t.verifyLessThan(norm(matrix(batch.ODT(k))-matrix(loop.ODT(k)),'fro'),1e-10);
    t.verifyLessThan(max(abs([batch.mag1(k),batch.mag2(k),batch.mag3(k)]- ...
        [loop.mag1(k),loop.mag2(k),loop.mag3(k)])),1e-10);
    for name = {'eV1','eV2','eV3'}
        a = batch.prop.(name{1});
        b = loop.prop.(name{1});
        t.verifyLessThan(acos(min(1,abs(dot(a(k),b(k))))),1e-5*degree);
        t.verifyLessThanOrEqual(a(k).z,1e-12);
    end
end
[grid,~,info] = gridCVA(ebsd,'computeBulkVector',false, ...
    'showProgress',false,'batchSize',2);
t.verifyEqual(info.numCandidateKernels,49);
t.verifyEqual(length(grid),49);
t.verifyTrue(all(isfinite(grid.mag1)));
end

function testGrainPGAAgainstMtex61Fixtures(t)
fixture = jsondecode(fileread(fullfile(t.TestData.root,'tests','fixtures','grains_expected.json')));
t.verifyEqual(numel(fixture.grains),24);
for k = 1:numel(fixture.grains)
    g = fixture.grains(k);
    matrices = permute(g.rotations,[2 3 1]);
    ori = orientation.byMatrix(matrices,crystalSymmetry(g.crystal_symmetry),specimenSymmetry('1'));
    [~,m,T] = PGA(ori);
    t.verifyLessThan(norm(matrix(mean(ori))-g.mean,'fro'),1e-9);
    t.verifyLessThan(norm(matrix(T)-g.tensor,'fro'),1e-10);
    t.verifyLessThan(max(abs(m(:)-g.eigenvalues(:))),1e-10);
end
end

function testGridAgainstMtex61Fixtures(t)
for mode = {'mean','center'}
    if strcmp(mode{1},'mean')
        filename = 'grid_expected.json';
    else
        filename = 'grid_center_expected.json';
    end
    fixture = jsondecode(fileread(fullfile(t.TestData.root,'tests','fixtures',filename)));
    input = fixture.input;
    symmetries = {'notIndexed'};
    for j = 1:numel(input.phases)
        symmetries{input.phases(j).id} = crystalSymmetry(input.phases(j).symmetry);
    end
    rotations = rotation.byMatrix(permute(input.rotations,[2 3 1]));
    ebsd = EBSD(vector3d(input.x(:),input.y(:),zeros(numel(input.x),1)), ...
        rotations(:),input.phase_ids(:),symmetries,struct(), ...
        'phaseMap',0:numel(symmetries)-1);
    ebsd.id = input.id(:);
    ebsd.prop.grainId = input.grain_ids(:);
    for sameGrain = [false,true]
        if sameGrain, expected = fixture.same_grain;
        else, expected = fixture.default; end
        [grid,~,info] = gridCVA(ebsd,'reference',mode{1}, ...
            'sameGrain',sameGrain,'computeBulkVector',false, ...
            'showProgress',false,'batchSize',2);
        t.verifyEqual(grid.id,[expected.id].');
        t.verifyEqual(grid.phaseId,[expected.phase_id].');
        for k = 1:length(grid)
            e = expected(k);
            center = find(info.centerId == e.id);
            members = info.windowId(center,:);
            t.verifyEqual(members(members>0).',e.member_ids(:));
            t.verifyEqual(logical(info.windowMask(center,:)), ...
                logical(e.window_mask(:).'));
            t.verifyLessThan(norm(matrix(grid.ODT(k))-e.tensor,'fro'),1e-10);
            t.verifyLessThan(max(abs([grid.mag1(k),grid.mag2(k),grid.mag3(k)]- ...
                e.eigenvalues(:).')),1e-10);
            t.verifyLessThan(abs(grid.kos(k)-e.kos),1e-8);
            t.verifyLessThan(norm(matrix(grid.meanRotation(k))-e.mean,'fro'),1e-8);
        end
    end
end
end

function tests = testGridCVAMtex61
tests = functiontests(localfunctions);
end
function setupOnce(t)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
t.TestData.root = root;
t.TestData.cwd=pwd; cd(fullfile(root,'tests','matlab')); addpath(root);
addpath(fullfile(root,'matlab','mtex61'),'-begin'); clear gridCVA PGA
end
function teardownOnce(t)
cd(t.TestData.cwd);
end
function testMaintainedPathExists(t)
t.verifyEqual(exist(fullfile(t.TestData.root,'matlab','mtex61','gridCVA.m'),'file'),2);
end
function testEmptyIndexedSupport(t)
e=syntheticGrid;
inputs={e([]),e}; inputs{2}.phaseId(:)=1;
for j=1:numel(inputs)
 [r,b,info]=gridCVA(inputs{j},'showProgress',false);
 verifyEmptyResult(t,r,b,info);
end
boundary=struct('boundary',struct('ebsdId',e.id));
[r,b,info]=gridCVA(e,boundary,'showProgress',false);
verifyEmptyResult(t,r,b,info);
end
function verifyEmptyResult(t,r,b,info)
t.verifyClass(r,'EBSD'); t.verifyEmpty(r.id); t.verifyTrue(isnan(b));
t.verifyEqual(info.numCandidateKernels,0); t.verifyEqual(info.numKernels,0);
t.verifySize(info.windowId,[0 9]);
for name={'CVA','eV1','eV2','eV3','mag1','mag2','mag3','kos','kax','meanRotation','ODT','kernelSize','centerPhase','centerGrainId'}
 t.verifyTrue(isfield(r.prop,name{1}),name{1});
end
end
function testIndependentSelection(t)
e=syntheticGrid;
expectedIds=[8;9;10;11;14;16;17;20;21;22;23;26;27;28;29];
for same=[false,true]
 [r,~,info]=gridCVA(e,'sameGrain',same,'computeBulkVector',false,'showProgress',false);
 t.verifyEqual(r.id,expectedIds);
 [centers,windows]=independentWindows(e,same);
 t.verifyEqual(info.centerId,centers); t.verifyEqual(info.windowId,windows);
 t.verifyEqual(info.windowMask,windows>0);
 % Center 16 lies against the phase boundary; opposite-phase column is masked.
 j=find(centers==16); t.verifyTrue(any(r.id==16));
 if same, t.verifyEqual(windows(j,:),[0 0 0 10 16 0 11 17 0]);
 else, t.verifyEqual(windows(j,:),[9 0 0 10 16 0 11 17 0]); end
end
end
function [centers,windows]=independentWindows(e,same)
% Select by physical square-grid offsets, independent of implementation IDs.
centers=e.id(e.x>min(e.x)&e.x<max(e.x)&e.y>min(e.y)&e.y<max(e.y)&e.phaseId>1);
windows=zeros(numel(centers),9); dx=[-1 0 1 -1 0 1 -1 0 1]; dy=[-1 -1 -1 0 0 0 1 1 1];
for j=1:numel(centers)
 c=find(e.id==centers(j));
 for k=1:9
  m=find(e.x==e.x(c)+dx(k)&e.y==e.y(c)+dy(k));
  if isempty(m)||e.phaseId(m)~=e.phaseId(c)||e.phaseId(m)==1,continue;end
  if same&&e.grainId(m)~=e.grainId(c),continue;end
  windows(j,k)=e.id(m);
 end
end
end
function testMixedPhaseWindowOracle(t)
e = syntheticGrid;
[r,b,info] = gridCVA(e,'computeBulkVector',false,'showProgress',false,'batchSize',2);
t.verifyTrue(isnan(b));
t.verifyGreaterThan(length(r),0);
t.verifyTrue(any(info.kernelSize < 9 & info.kernelSize >= 3));
for k = 1:length(r)
 j=find(info.centerId==r.id(k)); ids=info.windowId(j,:); ids=ids(ids>0);
 o=e('id',ids).orientations; ref=mean(o); v=log(o,ref,SO3TangentSpace.leftVector);
 M=v*v'; actual=double(r.ODT(k));
 t.verifyEqual(actual,M,'AbsTol',1e-12);
 t.verifyLessThan(angle(rotation(r.meanRotation(k)),rotation(ref)),1e-7);
 t.verifyEqual(r.kos(k),max(angle(o,ref)),'AbsTol',1e-12);
 t.verifyLessThan(angle(r.kax(k),mean(axis(o,ref))),1e-7);
 t.verifyEqual(r.centerPhase(k),e('id',r.id(k)).phaseId);
end
end
function testSameGrainAndChunkInvariance(t)
e=syntheticGrid;
[a,~,ia]=gridCVA(e,'sameGrain',true,'computeBulkVector',false,'showProgress',false,'batchSize',1);
[b,~,ib]=gridCVA(e,'sameGrain',true,'computeBulkVector',false,'showProgress',false,'batchSize',100);
t.verifyEqual(a.id,b.id); t.verifyEqual(a.mag1,b.mag1,'AbsTol',1e-12);
t.verifyEqual(ia.windowId,ib.windowId);
for j=1:numel(ia.centerId)
 ids=ia.windowId(j,:); ids=ids(ids>0);
 t.verifyTrue(all(e('id',ids).grainId==e('id',ia.centerId(j)).grainId));
end
end
function testCenterReferenceNamedAlternative(t)
e=syntheticGrid;
[r,~,info]=gridCVA(e,'reference','center','computeBulkVector',false,'showProgress',false);
t.verifyEqual(info.reference,'center');
for k=1:length(r)
 t.verifyLessThan(angle(rotation(r.meanRotation(k)),rotation(e('id',r.id(k)).orientations)),1e-7);
end
end
function testBoundedOptionsAndLowSpread(t)
e=syntheticGrid;
[r,~,info]=gridCVA(e,'maxKernels',4,'storeTensor',false,'computeBulkVector',false,'showProgress',false);
t.verifyEqual(info.numCandidateKernels,4); t.verifyFalse(isfield(r.prop,'ODT'));
t.verifyTrue(all(r.CVA.z<=0));
e.rotations=repmat(rotation.id,length(e),1);
r=gridCVA(e,'computeBulkVector',false,'showProgress',false);
t.verifyEqual(length(r),0);
end
function testOptionalBoundaryExclusion(t)
e=syntheticGrid;
boundary=struct('boundary',struct('ebsdId',[8;9]));
[~,~,info]=gridCVA(e,boundary,'computeBulkVector',false,'showProgress',false);
t.verifyFalse(any(ismember(info.windowId(:),[8;9])));
e.prop=rmfield(e.prop,'grainId');
t.verifyError(@()gridCVA(e,'sameGrain',true),'CVA:grid:MissingGrainId');
end
function testSharedCenterFixture(t)
file=fullfile(t.TestData.root,'tests','fixtures','grid_center_expected.json');
t.assertEqual(exist(file,'file'),2);
f=jsondecode(fileread(file)); base=jsondecode(fileread(fullfile(t.TestData.root,'tests','fixtures','grid_expected.json')));
t.verifyEqual(f.input,base.input);
e=syntheticGrid;
for same=[false,true]
 if same, kernels=f.same_grain; else, kernels=f.default; end
 [r,~,info]=gridCVA(e,'reference','center','sameGrain',same,'computeBulkVector',false,'showProgress',false);
 t.verifyEqual(r.id,[kernels.id].');
 [centers,windows]=independentWindows(e,same);
 t.verifyEqual(info.windowId,windows); t.verifyEqual(info.centerId,centers);
 for k=1:numel(kernels)
  q=kernels(k); j=find(centers==q.id); members=windows(j,:); members=members(members>0);
  t.verifyEqual(q.member_ids(:),members(:));
  o=e('id',members).orientations; ref=e('id',q.id).orientations;
  v=log(o,ref,SO3TangentSpace.leftVector); M=v*v';
  t.verifyEqual(q.tensor,M,'AbsTol',1e-12);
  t.verifyEqual(q.mean,matrix(ref),'AbsTol',1e-12);
  t.verifyEqual(q.kos,max(angle(o,ref)),'AbsTol',1e-12);
  a=mean(axis(o,ref)); t.verifyEqual(q.kax(:),[a.x;a.y;a.z],'AbsTol',1e-12);
  t.verifyEqual(q.eigenvalues(:),[r.mag1(k);r.mag2(k);r.mag3(k)],'AbsTol',1e-12);
  axes=[r.eV1(k),r.eV2(k),r.eV3(k)]; A=[axes.x(:),axes.y(:),axes.z(:)].';
  t.verifyEqual(abs(sum(q.eigenvectors.*A,1)),ones(1,3),'AbsTol',1e-10);
  t.verifyEqual(q.tensor*q.eigenvectors,q.eigenvectors*diag(q.eigenvalues),'AbsTol',1e-12);
 end
end
end
function testForsteriteLegacyParity(t)
e=EBSD.load('/Users/Z/Documents/MATLAB/mtex-6.1.1/data/ebsd/Forsterite.ctf');
g=gridify(e); e=g(1:20,1:20);
root=t.TestData.root; maintained=fullfile(root,'matlab','mtex61');
rmpath(maintained); clear gridCVA PGA
cleanup=onCleanup(@() addpath(maintained,'-begin'));
w=warning; wc=onCleanup(@() warning(w));
legacy=gridCVA(e);
addpath(maintained,'-begin'); clear gridCVA PGA
[new,~,info]=gridCVA(e,'computeBulkVector',false,'showProgress',false);
t.verifyEqual(new.id,legacy.id);
t.verifyEqual(new.phaseId,legacy.phaseId);
t.verifyEqual(new.kos,legacy.kos,'AbsTol',1e-10);
t.verifyLessThan(max(angle(new.meanRotation,legacy.meanRotation)),0.01*degree);
t.verifyLessThan(max(acos(min(1,abs(dot(new.CVA(:),legacy.CVA(:)))))),0.01*degree);
t.verifyEqual(new.mag1(:),legacy.mag1(:),'RelTol',1e-4,'AbsTol',1e-12);
t.verifyEqual(double(new.ODT),double(legacy.ODT),'RelTol',1e-4,'AbsTol',1e-12);
t.verifyLessThan(max(angle(new.kax,legacy.kax)),0.01*degree);
t.verifyTrue(info.samePhase);
end
function e=syntheticGrid
[x,y]=meshgrid(0:5,0:5); n=numel(x);
cs={ 'notIndexed',crystalSymmetry('mmm'),crystalSymmetry('m-3m')};
phase=2*ones(n,1); phase(x(:)>=3)=3; phase(15)=1;
o=rotation.byEuler((10+x(:)*.4+y(:)*.1)*degree,(20+y(:)*.3)*degree,(30+x(:).*y(:)*.07)*degree);
e=EBSD(vector3d(x(:),y(:),0),o,phase,cs,struct(),'phaseMap',[0 1 2]);
e.prop.grainId=1+(y(:)>=3)+2*(x(:)>=3);
end

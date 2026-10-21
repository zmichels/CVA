function tests = testCVAParity
tests = functiontests(localfunctions);
end
function testExplicitReferenceAndLegacy(t)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
new=fullfile(root,'matlab','mtex61');
t.verifyTrue(isfile(fullfile(new,'PGA.m')),'The maintained PGA path must exist.');
if ~isfile(fullfile(new,'PGA.m')), return; end
for symmetry={'1','222'}
 cs=crystalSymmetry(symmetry{1});
 base=orientation.byEuler(37*degree,63*degree,21*degree,cs);
 v=vector3d([1 0 0 1 2 -1],[0 1 0 2 -1 1],[0 0 1 -1 1 2]);
 for scale=[1,1e-7]
  ori=rotation.byAxisAngle(normalize(v),scale*[-3 2 1 4 -2 0.5]*degree).*base;
  ref=mean(ori);
  oldpwd=pwd; cleanup=onCleanup(@()cd(oldpwd)); cd(tempdir);
  addpath(root,'-begin'); clear PGA;
  t.assertEqual(which('PGA'),fullfile(root,'PGA.m'));
  [eo,mo,To]=PGA(ori);
  addpath(new,'-begin'); clear PGA;
  t.assertEqual(which('PGA'),fullfile(new,'PGA.m'));
  [en,mn,Tn]=PGA(ori,ref);
  t.verifyLessThanOrEqual(norm(matrix(Tn)-matrix(To),'fro'),1e-12);
  t.verifyEqual(mn,mo,'AbsTol',1e-12);
  if scale==1
   t.verifyLessThan(max(acos(min(1,abs(sum(matrix(en).*matrix(eo),1))))),1e-6);
  end
  t.verifyTrue(all(diff(mn)<=0));
  clear cleanup;
 end
end
end
function testCommittedFixtures(t)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
new=fullfile(root,'matlab','mtex61');
t.verifyTrue(isfile(fullfile(new,'PGA.m')));
if ~isfile(fullfile(new,'PGA.m')), return; end
oldpwd=pwd; cleanup=onCleanup(@()cd(oldpwd)); cd(tempdir);
addpath(new,'-begin'); clear PGA;
t.assertEqual(which('PGA'),fullfile(new,'PGA.m'));
raw=readmatrix(fullfile(root,'tests','fixtures','pga_rotations.csv'));
f=jsondecode(fileread(fullfile(root,'tests','fixtures','pga_expected.json')));
for k=1:numel(f.cases)
 c=f.cases(k); rows=raw(raw(:,1)==k,2:10); M=zeros(3,3,size(rows,1));
 for j=1:size(rows,1), M(:,:,j)=reshape(rows(j,:),3,3).'; end
 ori=orientation.byMatrix(M,crystalSymmetry(c.crystal_symmetry),specimenSymmetry('1'));
 [~,m,T]=PGA(ori);
 t.verifyEqual(matrix(mean(ori)),c.mean,'AbsTol',1e-12);
 t.verifyEqual(matrix(T),c.tensor,'AbsTol',1e-12);
 t.verifyEqual(m(:),c.eigenvalues(:),'AbsTol',1e-12);
end
end

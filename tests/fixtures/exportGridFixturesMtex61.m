function exportGridFixturesMtex61
%EXPORTGRIDFIXTURESMTEX61 Rebuild compact mean/center MATLAB oracle fixtures.
% From any folder, after startup_mtex:
% addpath('/absolute/path/to/CVA/tests/fixtures'); exportGridFixturesMtex61
% Requires MTEX 6.1.1. Uses direct per-window MTEX mean/log/tensor eig calls,
% not the maintained gridCVA window masks or batched calculation. The input
% is fixed by grid_expected.json; both output modes use exactly that input.
% The public legacy gridCVA comparison remains the separate bounded test.
folder=fileparts(mfilename('fullpath'));
original=jsondecode(fileread(fullfile(folder,'grid_expected.json')));
in=original.input;
cs={'notIndexed'};
for j=1:numel(in.phases), cs{in.phases(j).id}=crystalSymmetry(in.phases(j).symmetry); end
rot=rotation.byMatrix(permute(in.rotations,[2 3 1])); rot=rot(:);
e=EBSD(vector3d(in.x(:),in.y(:),0),rot,in.phase_ids(:),cs,struct(),'phaseMap',0:numel(cs)-1);
e.id=in.id(:); e.prop.grainId=in.grain_ids(:);
for mode={'mean','center'}
 refMode=mode{1};
 metadata=struct('matlab',version,'mtex','6.1.1','reference',refMode, ...
  'tangent_frame','left specimen','tensor_units','radians squared unnormalized', ...
  'matrix_layout','N by 3 by 3','oracle','independent coordinate windows and per-window MTEX calls');
 fixture=struct('metadata',metadata,'input',in);
 for same=[false,true]
  centers=find(e.x>min(e.x)&e.x<max(e.x)&e.y>min(e.y)&e.y<max(e.y)&e.isIndexed);
  records={}; dx=[-1 0 1 -1 0 1 -1 0 1]; dy=[-1 -1 -1 0 0 0 1 1 1];
  for c=centers(:).'
   ids=zeros(1,9);
   for k=1:9
    j=find(e.x==e.x(c)+dx(k)&e.y==e.y(c)+dy(k));
    if isempty(j)||~e.isIndexed(j)||e.phaseId(j)~=e.phaseId(c),continue;end
    if same&&e.grainId(j)~=e.grainId(c),continue;end
    ids(k)=e.id(j);
   end
   members=ids(ids>0); if numel(members)<3,continue;end
   o=e('id',members).orientations; o=o(:);
   if strcmp(refMode,'mean'),ref=mean(o);else,ref=e(c).orientations;end
   kos=max(angle(o,ref)); kax=mean(axis(o,ref));
   if kos<=.01*degree||isnan(kax),continue;end
   tangent=log(o,ref,SO3TangentSpace.leftVector);
   M=tangent*tangent'; [v,d]=eig(tensor(M,'rank',2));
   v=v([3 2 1]); d=d([3 2 1]); v(v.z>0)=-v(v.z>0);
   records{end+1,1}=struct('id',e.id(c),'phase_id',e.phaseId(c), ...
    'window_ids',ids,'window_mask',ids>0,'member_ids',members, ...
    'mean',matrix(ref),'tensor',M,'eigenvalues',d(:).', ...
    'eigenvectors',[v.x(:),v.y(:),v.z(:)].', ...
    'kos',kos,'kax',[kax.x,kax.y,kax.z]); %#ok<AGROW>
  end
  if same,fixture.same_grain=records;else,fixture.default=records;end
 end
 if strcmp(refMode,'mean'), name='grid_expected.json';else,name='grid_center_expected.json';end
 fid=fopen(fullfile(folder,name),'w');
 if fid<0,error('CVA:fixture:WriteFailed','Cannot write %s',name);end
 cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(fixture,PrettyPrint=true)); clear cleanup
 fprintf('Exported %s: %d default, %d same-grain kernels\n',name,numel(fixture.default),numel(fixture.same_grain));
end
end

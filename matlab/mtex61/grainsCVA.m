function [gCVA,bv,diagnostics]=grainsCVA(gCVA,ebsd,varargin)
%GRAINSCVA Grain CVA for MTEX 6.1.1 with exact per-grain MTEX means.
% method: 'batched' (default) or 'loop'. storeTensor defaults true;
% computeBulkVector defaults true (legacy antipodal 10 degree KDE).
p=inputParser;
p.addParameter('method','batched',@(x)any(strcmpi(x,{'loop','batched'})));
p.addParameter('computeBulkVector',true,@islogical);
p.addParameter('storeTensor',true,@islogical);
p.addParameter('showProgress',false,@islogical); p.parse(varargin{:}); opt=p.Results;
if isempty(ebsd.grainId), error('CVA:MissingGrainId','There is no ebsd.grainId. Run calcGrains first.'); end
gCVA=gCVA('indexed');
gCVA=gCVA(gCVA.numPixel>=3 & gCVA.GOS>0.01*degree);
ebsd=ebsd(gCVA); ebsd=ebsd('indexed');
[gid,~,group]=unique(ebsd.grainId); gCVA=gCVA('id',gid);
n=numel(gid); vals=nan(3,n); axes=nan(3,3,n); tensors=zeros(3,3,n); means=zeros(3,3,n);
for phase=unique(gCVA.phaseId).'
 grainIndices=find(gCVA.phaseId==phase); pixelIndices=find(ebsd.phaseId==phase);
 ori=ebsd(pixelIndices).orientations; localGroup=group(pixelIndices);
 refs=cell(numel(grainIndices),1);
 for j=1:numel(grainIndices)
  i=grainIndices(j); refs{j}=mean(ori(localGroup==i)); means(:,:,i)=matrix(refs{j});
 end
 if strcmpi(opt.method,'loop')
  for j=1:numel(grainIndices)
   i=grainIndices(j); [v,m,T]=PGA(ori(localGroup==i),refs{j});
   axes(:,:,i)=matrix(v); vals(:,i)=m(:); tensors(:,:,i)=matrix(T);
  end
 else
  % Construct references in this phase only; preserve crystal symmetry.
  ref=ori;
  for j=1:numel(grainIndices), ref(localGroup==grainIndices(j))=refs{j}; end
  tangents=matrix(log(ori,ref,SO3TangentSpace.leftVector));
  components=[1 1;2 2;3 3;1 2;1 3;2 3];
  for c=1:6
   x=components(c,1); y=components(c,2);
   sums=accumarray(localGroup,tangents(x,:).'.*tangents(y,:).',[n 1]);
   tensors(x,y,grainIndices)=sums(grainIndices); tensors(y,x,grainIndices)=sums(grainIndices);
  end
  for i=grainIndices(:).'
   [v,m]=eig(tensor(tensors(:,:,i),'rank',2)); [m,order]=sort(m,'descend');
   axes(:,:,i)=matrix(v(order)); vals(:,i)=m(:);
  end
 end
 if opt.showProgress, fprintf('CVA phase %d: %d grains\n',phase,numel(grainIndices)); end
end
v1=vector3d(reshape(axes(:,1,:),3,n)); v2=vector3d(reshape(axes(:,2,:),3,n)); v3=vector3d(reshape(axes(:,3,:),3,n));
v1(v1.z>0)=-v1(v1.z>0); v2(v2.z>0)=-v2(v2.z>0); v3(v3.z>0)=-v3(v3.z>0);
gCVA.prop.CVA=v1; gCVA.prop.eV1=v1; gCVA.prop.eV2=v2; gCVA.prop.eV3=v3;
gCVA.prop.mag1=vals(1,:); gCVA.prop.mag2=vals(2,:); gCVA.prop.mag3=vals(3,:);
if opt.storeTensor, gCVA.prop.ODT=tensor(tensors,'rank',2); end
bv=vector3d.empty;
if opt.computeBulkVector && n>0
 r=plotS2Grid('resolution',0.25*degree,'antipodal');
 kde=calcDensity([v1 -v1],r,'antipodal','halfwidth',10*degree); [~,I]=max(kde);
 bv=[r(I),-r(I)]; bv(bv.z>0)=[];
end
diagnostics=struct('method',opt.method,'grainIds',gid,'phaseIds',gCVA.phaseId,...
 'meanMatrices',means,'meanMethod','exact MTEX mean per grain','tangentFrame','specimen leftVector',...
 'tensorUnits','radians squared, unnormalised sum');
end

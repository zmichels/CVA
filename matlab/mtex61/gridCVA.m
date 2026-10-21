function [eCVA,bv,info] = gridCVA(e,varargin)
%GRIDCVA Mean-referenced, center-phase 3-by-3 CVA for MTEX 6.1.1.
% [eCVA,bv,info] = gridCVA(ebsd,[grains],Name,Value)
% Optional grains removes boundary-adjacent pixels, as in the public legacy
% function. sameGrain instead masks neighbors by center grain membership.
% reference='center' is a scientifically different explicit alternative.
% batchSize bounds post-mean tangent/tensor batches; all means use MTEX mean.
% ODT is an unnormalized sum of left tangent outer products, in radians^2.

grains=[];
if ~isempty(varargin) && ~(ischar(varargin{1}) || isstring(varargin{1}))
    grains=varargin{1}; varargin(1)=[];
end
p=inputParser;
addParameter(p,'sameGrain',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'reference','mean');
addParameter(p,'storeTensor',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'computeBulkVector',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'showProgress',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'maxKernels',inf,@(x)isnumeric(x)&&isscalar(x)&&x>0&&(isinf(x)||fix(x)==x));
addParameter(p,'batchSize',1024,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x);
parse(p,varargin{:}); opt=p.Results;
opt.reference=validatestring(opt.reference,{'mean','center'});
if ~isempty(grains) || opt.sameGrain
    if ~isfield(e.prop,'grainId') || isempty(e.grainId), error('CVA:grid:MissingGrainId','Run calcGrains first.'); end
end
if ~isempty(grains)
    boundaryIds=grains.boundary.ebsdId(:); boundaryIds(boundaryIds==0)=[];
    e('id',intersect(e.id,boundaryIds))=[];
end
if ~any(e.isIndexed(:))
    [eCVA,bv,info]=emptyGridResult(e,opt);
    return
end
g=gridify(e('indexed'));
g(~g.isIndexed).phase=nan;
g=gridify(g);
[nr,nc]=size(g.id);
centers=g.id(2:max(1,nr-1),2:max(1,nc-1));
mask=g.isIndexed(2:max(1,nr-1),2:max(1,nc-1));
centers=centers(mask & centers>0); centers=centers(:);
if numel(centers)>opt.maxKernels
    centers=centers(unique(round(linspace(1,numel(centers),opt.maxKernels))));
end
n=numel(centers);
% Legacy column-major center and neighbor order, in grid IDs.
windows=centers+[-nr-1,-1,nr-1,-nr,0,nr,-nr+1,1,nr+1];
phases=reshape(g('id',windows).phaseId,size(windows));
centerPhase=g('id',centers).phaseId; centerPhase=centerPhase(:);
keep=phases==centerPhase & reshape(g('id',windows).isIndexed,size(windows));
centerGrain=nan(n,1); windowGrain=nan(size(windows));
if isfield(g.prop,'grainId') && ~isempty(g.grainId)
    centerGrain=g('id',centers).grainId; centerGrain=centerGrain(:);
    windowGrain=reshape(g('id',windows).grainId,size(windows));
end
if opt.sameGrain, keep=keep & windowGrain==centerGrain; end
windows(~keep)=0; kernelSize=sum(keep,2);
refs=orientation.nan(n,1); kos=nan(n,1); kax=vector3d.nan(n,1);
valid=false(n,1);
for k=1:n
    if kernelSize(k)<3, continue; end
    o=g('id',windows(k,keep(k,:))).orientations;
    if strcmp(opt.reference,'mean'), ref=mean(o);
    else, ref=g('id',centers(k)).orientations; end
    spread=max(angle(o,ref));
    if spread<=.01*degree, continue; end
    refs(k)=ref; kos(k)=spread; kax(k)=mean(axis(o,ref)); valid(k)=true;
end
axes=vector3d.nan(3,n); mags=nan(3,n);
if opt.storeTensor, tensors=repmat(tensor(nan(3),'rank',2),n,1); end
for phase=unique(centerPhase(valid)).'
    selected=find(valid & centerPhase==phase);
    for first=1:opt.batchSize:numel(selected)
        ix=selected(first:min(first+opt.batchSize-1,end));
        W=windows(ix,:); [row,col]=find(W>0); row=row(:); col=col(:);
        member=W(sub2ind(size(W),row,col));
        o=g('id',member(:)).orientations; o=o(:);
        ref=orientation(rotation(refs(ix(row))),o.CS,o.SS); ref=ref(:);
        tangent=log(o,ref,SO3TangentSpace.leftVector);
        xyz=[tangent.x(:),tangent.y(:),tangent.z(:)];
        components=zeros(numel(ix),6); pairs=[1 1;1 2;1 3;2 2;2 3;3 3];
        for j=1:6
            components(:,j)=accumarray(row,xyz(:,pairs(j,1)).*xyz(:,pairs(j,2)),[numel(ix),1]);
        end
        M=zeros(3,3,numel(ix));
        for j=1:6
            M(pairs(j,1),pairs(j,2),:)=reshape(components(:,j),1,1,[]);
            M(pairs(j,2),pairs(j,1),:)=reshape(components(:,j),1,1,[]);
        end
        T=tensor(M,'rank',2); [v,d]=eig(T);
        axes(:,ix)=v([3 2 1],:); mags(:,ix)=d([3 2 1],:);
        if opt.storeTensor, tensors(ix)=T; end
        if opt.showProgress, fprintf('gridCVA phase %g: %d/%d kernels\n',phase,min(first+opt.batchSize-1,numel(selected)),numel(selected)); end
    end
end
axes(axes.z>0)=-axes(axes.z>0);
eCVA=g('id',centers);
eCVA.prop.CVA=axes(1,:).'; eCVA.prop.eV1=axes(1,:).';
eCVA.prop.eV2=axes(2,:).'; eCVA.prop.eV3=axes(3,:).';
eCVA.prop.mag1=mags(1,:).'; eCVA.prop.mag2=mags(2,:).'; eCVA.prop.mag3=mags(3,:).';
eCVA.prop.kos=kos; eCVA.prop.kax=kax; eCVA.prop.meanRotation=refs;
eCVA.prop.kernelSize=kernelSize; eCVA.prop.centerPhase=centerPhase;
eCVA.prop.centerGrainId=centerGrain;
if opt.storeTensor, eCVA.prop.ODT=tensors; end
accepted=valid & ~isnan(kax(:)) & ~isnan(axes(1,:).') & norm(axes(1,:).')>0;
eCVA=eCVA(accepted); eCVA=eCVA(~isnan(eCVA.rotations));
bv=vector3d.nan(1,1);
if opt.computeBulkVector && ~isempty(eCVA)
    r=plotS2Grid('resolution',.25*degree,'antipodal');
    density=calcDensity([eCVA.CVA -eCVA.CVA],r,'antipodal','halfwidth',10*degree);
    [~,j]=max(density); bv=[r(j),-r(j)]; bv(bv.z>0)=[];
end
info=struct('samePhase',true,'sameGrain',opt.sameGrain,'reference',opt.reference, ...
    'storedTensor',opt.storeTensor,'computedBulkVector',opt.computeBulkVector, ...
    'batchSize',opt.batchSize,'numCandidateKernels',n,'numKernels',length(eCVA), ...
    'centerId',centers,'centerPhase',centerPhase,'centerGrainId',centerGrain, ...
    'windowId',windows,'windowMask',keep,'windowPhase',phases, ...
    'windowGrainId',windowGrain,'kernelSize',kernelSize);
end

function [eCVA,bv,info]=emptyGridResult(e,opt)
% No indexed support is a normal empty result, before MTEX gridify.
eCVA=e([]); empty=zeros(0,1); v=vector3d.nan(0,1);
for name={'CVA','eV1','eV2','eV3','kax'}, eCVA.prop.(name{1})=v; end
for name={'mag1','mag2','mag3','kos','kernelSize','centerPhase','centerGrainId'}
    eCVA.prop.(name{1})=empty;
end
eCVA.prop.meanRotation=orientation.nan(0,1);
if opt.storeTensor, eCVA.prop.ODT=tensor(zeros(3,3,0),'rank',2); end
bv=vector3d.nan(1,1);
info=struct('samePhase',true,'sameGrain',opt.sameGrain,'reference',opt.reference, ...
    'storedTensor',opt.storeTensor,'computedBulkVector',opt.computeBulkVector, ...
    'batchSize',opt.batchSize,'numCandidateKernels',0,'numKernels',0, ...
    'centerId',empty,'centerPhase',empty,'centerGrainId',empty, ...
    'windowId',zeros(0,9),'windowMask',false(0,9),'windowPhase',zeros(0,9), ...
    'windowGrainId',zeros(0,9),'kernelSize',empty);
end

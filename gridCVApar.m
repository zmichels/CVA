% This function will perform principal geodesic analysis (PGA) of [3 x 3]
% kernels for a single phase. This analysis is similar to the grain-scale
% crystallographic vorticity axis (CVA) analysis, except that instead of
% using sets of orientations in whole grains, it uses a [3 x 3] window, and
% it is only applied to a single phase at a time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% Zachary D. Michels, Seth C. Kruckenberg, Joshua R. Davis, and Basil Tikoff
% Determining vorticity axes from grain-scale dispersion of
% crystallographic orientations Geology, G36868.1, first published on July
% 17, 2015, doi:10.1130/G36868.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% input:
%        ebsd:
%
% output:
%        eCVA (EBSD variable appended with PGA results, such as the
%        following):
%
%            eV:      PGA eigenvectors as vector3d  ––  eV(1,:) = cva
%          mags:      Eigenvalue of eigenvectors    ––  mags(1,:) for cva's
%            bv:      Preferred cva direction as vector3d
%           eId:      ids of ebsd at center of 3x3 window
%
%
% example usage:
%
% % load some data

% mtexdata forsterite
%
% % compute the kernel CVA analysis (use either of the following syntax)
% [eCVA,bv] = gridCVA(ebsd)            % all points
% [eCVA,bv] = gridCVA(ebsd, grains)    % excluding points adjacent to
%                                             grain boundaries
%
%
% % plot results and best-fit/preferred cva vector
% figure,
% plot(eCVA.CVA,'antipodal','lower','smooth')
% hold on
% plot(bv,'antipodal','lower','Marker','^','MarkerSize',10,...
%     'MarkerFaceColor','k','MarkerEdgeColor','w')


%%
function [eCVA,bv] = gridCVApar(ePhase,varargin)


%%
narginchk(1,2)
nargin;
varargin;
grains = [varargin{1}];

if nargin > 1
    % check if grainID exists
    if isempty(ePhase.grainId)
        
        error('There is no ebsd.grainId. Run calcGrains first.')
        
    end
    
    % get ebsd IDs at grain on either side of grain boundaries
    eId = grains.boundary.ebsdId(:);
    
    % remove any zeros
    eId(eId==0)=[];
    
    % remove the identified pixels from the dataset
    ePhase('id',intersect(ePhase.id,eId)) = [];
end

%%
% gridify
egrid = ePhase('indexed').gridify;

% ebsd ids in grid/matrix
ids = egrid.id;

% size of matrix
[a1,b1] = size(ids);

% window width
w = 3;

% indices/ids of center points
row = 2:1:max(a1)-1;
col = 2:1:max(b1)-1;

inds = ids(2:a1-1,2:b1-1);
eId = inds(:);

phases = unique(egrid('indexed').phase);

%% initialize window
num = length(eId);

win1 = zeros(w,w,num);

for s = 1:num
    % s = 1;
    c = inds(s);
    
    win1(:,:,s) =   [   c-a1-1     c-1    c+a1-1;
        
                        c-a1       c      c+a1;
                     
                        c-a1+1     c+1    c+a1+1   ];
                 
    
end



%% assign rotations
oRot = egrid(win1).rotations;
pID = egrid(win1).phase;
CSList = egrid.CSList;
mineralList = egrid.mineralList;

%% pre-allocate
eV = [vector3d.nan(1,num); vector3d.nan(1,num); vector3d.nan(1,num)];
mags = nan(3,num);
kos = nan(size(eId));
kax = vector3d.nan(1,num);
meanRotation = orientation.nan(num,1);
T = repmat(tensor(nan(3,3),'rank',2),[num,1]);


%% parpool setup
if isempty(gcp('nocreate'))
availableGPUs = gpuDeviceCount("available");
pool = parpool('local',availableGPUs);
end

%% analysis loop
% for keeping track of progress in for loop:
div=round(num/20);
count=div;

fprintf('\n%i kernels\n',num)
% fprintf('\n%i%% done\n',0)

WaitMessage = parfor_wait(num,'Waitbar',true);

parfor n = 1:num
    WaitMessage.Send;
    pInd = pID(:,:,n)==pID(2,2,n)&pID(:,:,n)>0;
    rots = oRot(pInd(1,:),pInd(:,1),n);
    rots = rots(~isnan(rots(:)));

    if length(rots)>2 && max(angle(rots,mean(rots)))>.01*degree
        
        [eV(:,n),mags(:,n),T(n,:)] = PGA(rots);
     
        % kernel mean orientation
        meanRotation(n,:) = mean(rots);
        % kernel orientation spread (KOS - like mis2mean for kernel)
        kos(n,:) = max(angle(rots,mean(rots)));
        % kernel mean KOS axis
        kax(n,:) = mean(axis(rots,mean(rots)));
        
    end
%     % Keep track of for loop progress and print to consoloe screen:
%         perc=round(n/num*100);
%         if n==count            
%             fprintf('\n%i%% done...\n',perc)
%             count=count+div;
%         end

end
WaitMessage.Destroy
% project to lower hemisphere
eV(eV.z>0)=-eV(eV.z>0);


%% append ebsd variable
eCVA = egrid(eId);
eCVA.prop.CVA           = eV(1,:);
eCVA.prop.eV1           = eV(1,:);
eCVA.prop.eV2           = eV(2,:);
eCVA.prop.eV3           = eV(3,:);
eCVA.prop.mag1          = mags(1,:);
eCVA.prop.mag2          = mags(2,:);
eCVA.prop.mag3          = mags(3,:);
eCVA.prop.kos           = kos;
eCVA.prop.kax           = kax;
eCVA.prop.meanRotation  = meanRotation;
eCVA.prop.cvaTensors    = T;



%% Handle results
% identify null solutions
cond = (norm(eCVA.CVA)==0 | isnan(eCVA.mag1) | isnan(eCVA.CVA) | isnan(eCVA.kax));

% apply condition
eCVA(cond) = [];


%% Kernel Density Estimation to get a best fit "bulk" vorticity vector.
% Define a kernel density estimation with specified halfwidth. MTEX default
% uses the de la Vallee Poussin kernel:

r = plotS2Grid('resolution',0.25*degree,'antipodal');
kde = calcDensity([eCVA.CVA -eCVA.CVA],r,'antipodal','halfwidth',10*degree);
[~,I]=max(kde);

% get vector and negated vector (antipodal) of best-fit axis:
bv=[r(I),-r(I)];
bv(bv.z>0) = [];

%% close / delete parpool
delete(pool)
end

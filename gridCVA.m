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
%           ODT:      Orientation dispersion tensor result from PGA
%            eV:      ODT eigenvectors as vector3d  ––  eV(1,:) = cva
%          mags:      Eigenvalues of eigenvectors    ––  mags(1,:) for cva's




% example usage:

% % load some data

% mtexdata forsterite

% % compute the kernel CVA analysis (use either of the following syntax)
% [eCVA,bv] = gridCVA(ebsd)            % all points
% [eCVA,bv] = gridCVA(ebsd, grains)    % excluding points adjacent to



% % plot results and best-fit/preferred cva vector
% figure,
% plot(eCVA.CVA,'antipodal','lower','smooth')
% hold on
% plot(bv,'antipodal','lower','Marker','^','MarkerSize',10,...
%     'MarkerFaceColor','k','MarkerEdgeColor','w')


%%
function [eCVA,bv] = gridCVA(ebsd,varargin)

warning off
%%
narginchk(1,2)
nargin;
if nargin > 1
    varargin;
    grains = [varargin{1}];


    % check if grainID exists
    if isempty(ebsd.grainId)

        error('There is no ebsd.grainId. Run calcGrains first.')

    end

    % get ebsd IDs at grain on either side of grain boundaries
    eIdB = grains.boundary.ebsdId(:);

    % remove any zeros
    eIdB(eIdB==0)=[];

    % remove the identified pixels from the dataset
    ebsd('id',intersect(ebsd.id,eIdB)) = [];
end

%%
% gridify
egrid = ebsd.gridify;

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

phases = unique(egrid.phase);

%% initialize window
num = length(eId);

win1 = zeros(num,w*w);

for s = 1:num
    % s = 1;
    c = inds(s);

    win1(s,:) =   [   c-a1-1     c-1    c+a1-1    c-a1       c      c+a1  c-a1+1     c+1    c+a1+1   ];


end


%% assign rotations
% nan not-indexed
egrid(~(egrid.isIndexed)).phase = nan;
% phase IDs
pID = egrid(win1).phase;
% Crystal Symmetry list
CSList = egrid.CSList;
% mieral name list
mineralList = egrid.mineralList;
% phase numbers
p = egrid('id',eId).phase;
% for logical query
pind1 = pID(:)==reshape(repmat(p,1,w*w),[w*w*num,1]);
% logical index by phase of points in the kernel window with a phase that
% matches the central point
pind = reshape(pind1,size(pID));
% get the indexed phase numbers
pp = pind.*pID;
% ids of points in the windows
winId = egrid('id',win1).id;
% use logical index to leave only ones of same phase
winId(~pind) = 0;


%% pre-allocate
eV = [vector3d.nan(1,num); vector3d.nan(1,num); vector3d.nan(1,num)];
mags = nan(3,num);
kos = nan(size(eId));
kax = vector3d.nan(1,num);
meanOrientation = orientation.nan(num,1);
T = repmat(tensor(nan(3,3),'rank',2),[num,1]);


%% analysis loop
% for keeping track of progress in for loop:
div=round(num/20);
count=div;

fprintf('\n%i kernels\n',num)
fprintf('\n%i%% done\n',0)

for n = 1:num
    
    
    if sum(~isnan(egrid(winId(n,winId(n,:)>0))))>2

        % orientations of same phase in the kernel
        o = egrid(winId(n,winId(n,:)>0)).orientations;

        if length(o)>2 && max(angle(o,mean(o)))>.01*degree

            [eV(:,n),mags(:,n),T(n)] = PGA(o);

            % kernel mean orientation
            meanOrientation(n) = mean(o);
            % kernel orientation spread (KOS - like mis2mean for kernel)
            kos(n) = max(angle(o,mean(o)));
            % kernel mean KOS axis
            kax(n) = mean(axis(o,mean(o)));

        end

    end
    % Keep track of for loop progress and print to consoloe screen:
    perc=round(n/num*100);
    if n==count
        fprintf('\n%i%% done...\n',perc)
        count=count+div;
    end
end

% project to lower hemisphere
eV(eV.z>0)=-eV(eV.z>0);


%% append ebsd variable
eCVA = egrid(eId);
eCVA.prop.CVA = eV(1,:);
eCVA.prop.eV1 = eV(1,:);
eCVA.prop.eV2 = eV(2,:);
eCVA.prop.eV3 = eV(3,:);
eCVA.prop.mag1 = mags(1,:);
eCVA.prop.mag2 = mags(2,:);
eCVA.prop.mag3 = mags(3,:);
eCVA.prop.kos = kos;
eCVA.prop.kax = kax;
eCVA.prop.meanRotation = meanOrientation;
eCVA.prop.ODT = T;



%% Handle results

% identify null solutions
cond1 = (norm(eCVA.CVA)==0 | isnan(eCVA.mag1) | isnan(eCVA.CVA) | isnan(eCVA.kax));

% apply condition
eCVA = eCVA(~cond1);


eCVA = eCVA(~isnan(eCVA.rotations));



%% Kernel Density Estimation to get a best fit "bulk" vorticity vector.
% Define a kernel density estimation with specified halfwidth. MTEX default
% uses the de la Vallee Poussin kernel:

r = plotS2Grid('resolution',0.25*degree,'antipodal');
kde = calcDensity([eCVA.CVA -eCVA.CVA],r,'antipodal','halfwidth',10*degree);
[~,I]=max(kde);

% get vector and negated vector (antipodal) of best-fit axis:
bv=[r(I),-r(I)];
bv(bv.z>0) = [];



end

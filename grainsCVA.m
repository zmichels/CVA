% This function will calculate CVA for grains. CVA is based ona principal
% geodesic analysis (PGA) which identifies the best-fit rotation to
% describe the dispersion of orientations in a set of orientations.
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
%        grains:        grainSet
%          ebsd:        ebsd (with grainID)
% 
% 
% output:
%          gCVA:        grainSet appended with CVA resuls for each grain
%
%
% example usage:
%
% % load some data
% mtexdata forsterite
%
% % compute grains
% [grains,ebsd.grainId]=calcGrains(ebsd)
%
% % compute CVA for grains
% [eV,mags,bv,eId] = gridCVA(ebsd('f'))
%
% % plot results and best-fit/preferred cva vector
% figure,
% plot(gCVA.CVA,'antipodal','lower','smooth')
% hold on
% plot(bv,'antipodal','lower','Marker','^','MarkerSize',10,...
%     'MarkerFaceColor','k','MarkerEdgeColor','w')


function [gCVA, bv] = grainsCVA(gCVA, ebsd)

% check if grainID exists
if isempty(ebsd.grainId)
     
    error('There is no ebsd.grainId. Run calcGrains first.')
     
end

%% Setup: select grains and initialize variables for sake of loops
gCVA = gCVA('indexed');
gCVA = gCVA(gCVA.grainSize>=3&gCVA.GOS>0.01*degree);
ebsd = ebsd(gCVA);


ebsd(ebsd.phase<1) = [];

[gid,~,eindex] = unique(ebsd.grainId);


% number of grains:
num=length(gid);

% for keeping track of progress in for loop:
div=round(num/10);
count=div;

% pre-allocate
eVg = [vector3d(nan(3,num));vector3d(nan(3,num));vector3d(nan(3,num))];
mags = nan(3,num);

%% anlysis loop
fprintf('\n\n%i grains\n\n',num)

for n = 1:length(gid)
    
    [eVg(:,n),mags(:,n)] = PGA(ebsd(eindex==n).orientations);
    
    % Keep track of for loop progress and print to consoloe screen:
    perc=round(n/num*100);
        if n==count            
            fprintf('\n\n%i percent done...\n\n',perc)
            count=count+div;
        end
    
end

%% project to lower hemisphere
eVg(eVg.z>0)=-eVg(eVg.z>0);


%% Kernel Density Estimation to get a best fit "bulk" vorticity vector.
% Define a kernel density estimation with specified halfwidth. MTEX default
% uses the de la Vallee Poussin kernel:

r = plotS2Grid('resolution',0.25*degree,'antipodal');
kde = calcDensity([eVg(1,:) -eVg(1,:)],r,'antipodal','halfwidth',10*degree);
[~,I]=max(kde);

% get vector and negated vector (antipodal) of best-fit axis:
bv=[r(I),-r(I)];
bv(bv.z>0) = [];


%% Append the grainset
gCVA.prop.CVA = eVg(1,:);
gCVA.prop.eV1 = eVg(1,:);
gCVA.prop.eV2 = eVg(2,:);
gCVA.prop.eV3 = eVg(3,:);
gCVA.prop.mag1 = mags(1,:);
gCVA.prop.mag2 = mags(2,:);
gCVA.prop.mag3 = mags(3,:);

%% CVA inversion to crystal reference frame
phases = unique(gCVA.phase);

for j = 1:length(phases)
    pInd = gCVA.phase==phases(j);
gCVA(pInd).prop.inv_CVA = round(gCVA(pInd).meanOrientation.\gCVA(pInd).CVA');
gCVA(pInd).prop.inv_eV2 = round(gCVA(pInd).meanOrientation.\gCVA(pInd).eV2');
gCVA(pInd).prop.inv_eV3 = round(gCVA(pInd).meanOrientation.\gCVA(pInd).eV3');

end
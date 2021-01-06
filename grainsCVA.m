function [gCVA, bulkVort] = grainsCVA(gCVA, eCVA)

% calculate CVA for grains
% input:
%      ebsd (with grainID)
% 
% 
%   output:
%      gCVA:          grainSet appended with CVA resuls for each grain
%
%
%   usage:
%   [grains,ebsd.grainId]=calcGrains(ebsd)
%   [eVg, gid] = grainsPGA(ebsd)


%   check if grainID exists
if isempty(eCVA.grainId)
     
    error('There is no ebsd.grainId. Run calcGrains first.')
     
end

%% Setup: select grains and initialize variables for sake of loops
gCVA = gCVA('indexed');
gCVA = gCVA(gCVA.grainSize>=3&gCVA.GOS>0.01*degree);
eCVA = eCVA(gCVA);


eCVA(eCVA.phase<1) = [];

[gid,~,eindex] = unique(eCVA.grainId);


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
    
    [eVg(:,n),mags(:,n)] = PGA(eCVA(eindex==n).orientations);
    
    % Keep track of for loop progress and print to consoloe screen:
    perc=round(n/num*100);
        if n==count            
            fprintf('\n\n%i percent done...\n\n',perc)
            count=count+div;
        end
    
end

%% project to lower hemisphere
eVg(eVg.z>0)=-eVg(eVg.z>0);


eV1 = eVg(1,:);
%% Try Kernel Density Estimation to get a best fit "bulk" vorticity vector.
% Define a kernel density estimation with specified halfwidth. MTEX defaul

% uses the de la Vallee Poussin kernel:
r = plotS2Grid('resolution',0.25*degree,'antipodal');
kde = calcDensity([eV1 -eV1],r,'antipodal','halfwidth',10*degree);
[~,I]=max(kde);

% get vector and negated vector (antipodal) of best-fit axis:
bulkVort=[r(I),-r(I)];
bulkVort(bulkVort.z>0) = [];


%% Append the grainset
gCVA.prop.CVA = eV1;
gCVA.prop.eV1 = eV1;
gCVA.prop.eV2 = eVg(2,:);
gCVA.prop.eV3 = eVg(3,:);
gCVA.prop.mag1 = mags(1,:);
gCVA.prop.mag2 = mags(2,:);
gCVA.prop.mag3 = mags(3,:);


end
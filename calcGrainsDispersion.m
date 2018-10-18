function [ gCVA , bulkVort ]  =  calcGrainsDispersion ( gCVA , eCVA )
% written by Zachary D. Michels; February, 2015
% edited 2/23/2018 for compatibility with MTEX 5.0 

%% References
% Zachary D. Michels, Seth C. Kruckenberg, Joshua R. Davis, and Basil Tikoff
% Determining vorticity axes from grain-scale dispersion of
% crystallographic orientations Geology, G36868.1, first published on July
% 17, 2015, doi:10.1130/G36868.1


%% NOTES TO USERS and REFERENCE INFORMATION

% !*** IMPORTANT ***! This function requires MTEX v5.x Some alteration may
% be required for use with subsequent versions of MTEX. However, MOST of
% the calculations do not require MTEX specific functions (mainly MTEX is
% used for passing intragranular orientation solutions to the custom
% scripts), and the code should be easily adapted. For consultation
% regarding this analysis, please contact Z. Michels through GitHub
% repository:(https://github.com/zmichels/CVA.git)

% NOTE: All grain objects passed to this function must contain more than 1
% indexed orientation solution... a good rule of thumb is >3 solutions. The
% comments immediately below give an example of how this might be
% accomplished using MTEX 4.0 to make a new set of grains, 'g', which
% satisfy these conditions:

% ** EXAMPLE **
% [grains,ebsd.grainId,ebsd.mis2mean] = 
% calcGrains(ebsd,'angle',10*degree,'keepNotIndexed');
%
% condition = grains.grainSize >= 3 & grains.phase>0 & grains.GOS>0.001*degree;
%
% g=grains(condition);


% ************************************************************************
%% References
% Zachary D. Michels, Seth C. Kruckenberg, Joshua R. Davis, and Basil Tikoff
% Determining vorticity axes from grain-scale dispersion of
% crystallographic orientations Geology, G36868.1, first published on July
% 17, 2015, doi:10.1130/G36868.1

% ************************************************************************



%% Setup: select grains and initialize variables for sake of loops
gCVA = gCVA('indexed');
gCVA = gCVA(gCVA.grainSize>=3&gCVA.GOS>0.01*degree);
eCVA = eCVA(gCVA);

% number of grains:
numg=size(gCVA,1)

% for keeping track of progress in for loop:
div=round(numg/10);
count=div;

% initialize variables (MTEX vector3d object for vorticity vectors)
vorts = vector3d();
D = rotation();
eV1 = vector3d();
eV2 = vector3d();
eV3 = vector3d();
mags=[];



%% Analysis Loop:
% This loop passes data from one grain at a time to the intragranular
% (single grain) dispersion script. The intragranular dispersion function
% applies a principal geodesic analysis and returns a best-fit axis that
% describes rotational crystallographic distortion within each grain. A
% vector for each grain-scale vorticity axis is stored in the variable
% 'vorts'.

[vorts(1),D(1),eV1(1), eV2(1), eV3(1),mags(:,1),T]=intragranularDispersion(eCVA(gCVA(1)));


for i = 2:numg
    % Pass one grain object from the GrainSet to the function for
    % intragranular (single-grain) dispersion analysis using Principal
    % Geodesic Analysis
    
    [vorts(i),D(i),eV1(i), eV2(i), eV3(i),mags(:,i),t]=intragranularDispersion(eCVA(gCVA(i)));
    T(i) = t;
    mags=real(mags);
    
    % Keep track of for loop progress and print to consoloe screen:
    perc=round(i/numg*100);
        if i==count            
            fprintf('\n\n%i percent done...\n\n',perc)
            count=count+div;
        end
end


%% Remove nonsensical vectors (i.e., [0,0,0]) 
% for sake of kernel density estimation (below):

[x,y,z]=double(vorts);
IzeroX = find(~real(x))';
IzeroY = find(~real(y))';
IzeroZ = find(~real(z))';
intXYZ=intersect(intersect(IzeroX,IzeroY,'rows'),IzeroZ);
vorts(intXYZ)=[];


vorts(vorts.z>0) = -vorts(vorts.z>0);


%% Try Kernel Density Estimation to get a best fit "bulk" vorticity vector.
% Define a kernel density estimation with specified halfwidth. MTEX defaul



% uses the de la Vallee Poussin kernel:
r = plotS2Grid('resolution',0.25*degree,'antipodal');
kde = calcDensity([vorts -vorts],r,'antipodal','halfwidth',10*degree);
[~,I]=max(kde);

% get vector and negated vector (antipodal) of best-fit axis:
bulkVort=r(I);


%% Append the grainset
gCVA.prop.CVA = vorts;
gCVA.prop.D = D;
gCVA.prop.eV1 = eV1;
gCVA.prop.eV2 = eV2;
gCVA.prop.eV3 = eV3;
gCVA.prop.mag1 = mags(1,:);
gCVA.prop.mag2 = mags(2,:);
gCVA.prop.mag3 = mags(3,:);
gCVA.prop.T = T;

% angle CVA to zvector
gCVA.prop.ang2z = angle(gCVA.CVA,zvector,'antipodal')./degree;

gCVA.prop.ind = ones(size(gCVA.CVA));

% Filter
gCVA = gCVA((gCVA.mag1./gCVA.mag2)>2);


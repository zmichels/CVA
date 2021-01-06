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
%        ePhase:
% 
% output:
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
% % compute the kernel CVA analysis
% [eV,mags,bv,eId] = gridCVA(ebsd('f'))
%
% % plot results and best-fit/preferred cva vector
% figure,
% plot(eV(1,:),'antipodal','lower','smooth')
% hold on
% plot(bv,'antipodal','lower','Marker','^','MarkerSize',10,...
%     'MarkerFaceColor','k','MarkerEdgeColor','w')


%%
function [eV,mags,bv,eId] = gridCVA(ePhase)
%%

% gridify
ePhase = ePhase.gridify;

% ebsd ids in grid/matrix
ids = ePhase.id;

% size of matrix
[a1,b1] = size(ids);

% window width
w = 3;

% make sure the matrix of ids is divisible by w
a = w*floor(a1/w);
b = w*floor(b1/w);
ids = ids(1:a,1:b);

% break up the matrix into cells w x w cells
win1 = mat2cell(ids,[w*ones(a/w,1)'],[w*ones(b/w,1)']);



% list of cells
wins = win1(:);

% number of grains:
num=length(wins);


% for keeping track of progress in for loop:
div=round(num/10);
count=div;

% pre-allocate
eV = [vector3d(nan(3,num));vector3d(nan(3,num));vector3d(nan(3,num))];
mags = nan(3,num);
eId = nan(1,num);

%% analysis loop
fprintf('\n\n%i kernels\n\n',num)

for n = 1:length(wins)
    m = wins{n};
    if sum(ePhase(m).phase)>0
        o = ePhase(m).orientations;
        o = o(~isnan(o));
        
        if length(o)>=3
            [eV(:,n),mags(:,n)] = PGA(o);
            eId(n) = m(ceil(numel(m)/2));
        end
    end
    
    % Keep track of for loop progress and print to consoloe screen:
    perc=round(n/num*100);
        if n==count            
            fprintf('\n\n%i percent done...\n\n',perc)
            count=count+div;
        end
end

% project to lower hemisphere
eV(eV.z>0)=-eV(eV.z>0);

%% Handle results
% % assign 'cva' from eignevector results
% eV1 = eV(1,:);
% 
% % and mags
% mag = mags(1,:);

% identify null solutions
cond = (norm(eV(1,:))==0 | isnan(mags(1,:)));

% apply condition
eV(:,cond) = [];
eId(cond) = [];
mags(:,cond) = [];

%% Kernel Density Estimation to get a best fit "bulk" vorticity vector.
% Define a kernel density estimation with specified halfwidth. MTEX default
% uses the de la Vallee Poussin kernel:

r = plotS2Grid('resolution',0.25*degree,'antipodal');
kde = calcDensity([eV(1,:) -eV(1,:)],r,'antipodal','halfwidth',10*degree);
[~,I]=max(kde);

% get vector and negated vector (antipodal) of best-fit axis:
bv=[r(I),-r(I)];
bv(bv.z>0) = [];

end

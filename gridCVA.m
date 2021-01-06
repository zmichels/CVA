%% GRID CVA

% This function will perform a sliding-window principal geodesic analysis
% (CVA) in which the grain boundaries are ignored ? so that only
% intragranular orientations are considered. This function has been tested
% to work with MTEX v5.3

% NOTE: The sliding-window approach will take longer than a whole-grain
% approach due to the increase in analyses that are required. The script
% includes estimates for remaining time during updates, but these are often
% incorrect and usually a significant underestimate.


function [cva,mag,bv,eId] = gridCVA(ePhase)
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
% assign 'cva' from eignevector results
cva = eV(1,:);

% and mags
mag = mags(1,:);

% identify null solutions
cond = (norm(cva)==0 | isnan(mag));

% apply condition
cva(cond) = [];
eId(cond) = [];
mag(cond) = [];

%% Try Kernel Density Estimation to get a best fit "bulk" vorticity vector.
% Define a kernel density estimation with specified halfwidth. MTEX defaul

% uses the de la Vallee Poussin kernel:
r = plotS2Grid('resolution',0.25*degree,'antipodal');
kde = calcDensity([cva -cva],r,'antipodal','halfwidth',10*degree);
[~,I]=max(kde);

% get vector and negated vector (antipodal) of best-fit axis:
bv=[r(I),-r(I)];
bv(bv.z>0) = [];

end

%% GRID CVA

% This function will perform a sliding-window principal geodesic analysis
% (CVA) in which the grain boundaries are ignored ? so that only
% intragranular orientations are considered. This function has been tested
% to work with MTEX v5.3

% NOTE: The sliding-window approach will take longer than a whole-grain
% approach due to the increase in analyses that are required. The script
% includes estimates for remaining time during updates, but these are often
% incorrect and usually a significant underestimate.


function [eCVA] = gridCVA(e,pad,grain)
%%
gb = grain.boundary;
id1 = gb.ebsdId(:,1);
id2 = gb.ebsdId(:,2);

uIds = unique([id1(:);id2(:)]);
e.prop.newId = 1:length(e);

e(intersect(e.newId,uIds)) = [];


e = e.gridify;

lims = size(e)-1;
ex = unique(e.x);
ey = unique(e.y);

stepx = abs(ex(2)-ex(1));
stepy = abs(ey(2)-ey(1));

ecellArea = 3*stepx*3*stepy;

% Set overlap
% pad = 3;

% initialize variables (MTEX vector3d object for vorticity vectors)
CVAgrid = vector3d();
Dgrid = rotation();
eV1grid = vector3d();
eV2grid = vector3d();
eV3grid = vector3d();
magsgrid=[];

stepsInX = size(e,1);
stepsInY = size(e,2);
%%
% Need to select orientations of same phase, and make sure selected
% measurements comprise >3 orientations
i = 1;
numg = stepsInX*stepsInY;
inc = 10;
% for keeping track of progress in for loop:
div=round(numg/inc);
count=div;
n = 1;

centx = nan(size(1+pad:stepsInX-pad));
centy = nan(size(1+pad:stepsInY-pad));

% eGrid = e(1+pad:stepsInX-pad,1+pad:stepsInY-pad);
%%
tic
tocInt=toc;
toc
for j = length(centx)
    for k = length(centy)
        eWin = e([j-pad:j+pad],[k-pad:k+pad]);
        if size(eWin(eWin.phase>0),1) >= 3
            
            [eV,mags] = PGA(eWin(eWin.phase>0).orientations);
            if mags(1,i)>0
            centx(i)= e(j,k).x;
            centy(i) = e(j,k).y;
            magsgrid=real(magsgrid);
            id(i) = e(j,k).id;
            i = i+1;
            
            end
        end
        
        % Keep track of for loop progress and print to consoloe screen:
        perc=round(n/numg*100);
        
        if n==count
            fprintf('\n\n%i percent done ...\n',perc)
            toc
            fprintf('\n     %0.2f minutes total',roundn(toc/60,-2))
            fprintf('\n     %0.2f minutes since last',roundn(toc-tocInt,-2)/60)
            fprintf('\n     %0.2f hours remaining (%0.2f minutes)\n\n\n',roundn(toc-tocInt,-2)/60/60*((100-perc)/inc),roundn(toc-tocInt,-2)/60*((100-perc)/inc))
            tocInt=toc;
            
            count=count+div;
        end
        n = n+1;
    end
    
    
    % project to lower hemisphere
    eV(eV.z>0)=-eV(eV.z>0);


    eV1 = eV(1,:);


    eCVA = e(id);
    eCVA.prop.CVA = eV1;
    eCVA.prop.eV1 = eV1;
    eCVA.prop.eV2 = eV(2,:);
    eCVA.prop.eV3 = eV(3,:);
    eCVA.prop.mag1 = magsgrid(1,:);
    eCVA.prop.mag2 = magsgrid(2,:);
    eCVA.prop.mag3 = magsgrid(3,:);
    eCVA.prop.centx = centx;
    eCVA.prop.centy = centy;
    
end

%% GRID CVA

function [CVAgrid, Dgrid, Tgrid, eV1grid, eV2grid, eV3grid, magsgrid, centx, centy, eCVA] = gridCVA(e,pad)
%%
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
Tgrid = tensor([0 0 0; 0 0 0; 0 0 0],'rank',2);

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


% eGrid = e(1+pad:stepsInX-pad,1+pad:stepsInY-pad);
%%
tic
tocInt=toc;
toc
for j = 1+pad:stepsInX-pad
    for k = 1+pad:stepsInY-pad
        eWin = e([j-pad:j+pad],[k-pad:k+pad]);
        if size(eWin(eWin.phase>0),1) >= 3
            
            [CVAgrid(i),Dgrid(i),eV1grid(i), eV2grid(i), eV3grid(i),magsgrid(:,i),t]=intragranularDispersion(eWin(eWin.phase>0));
            if magsgrid(1,i)>0
            centx(i)= e(j,k).x;
            centy(i) = e(j,k).y;
            Tgrid(i) = t;
            magsgrid=real(magsgrid);
            id(i) = e(j,k).id;
            i = i+1;
            
            end
            
            
          
        else
            %             CVA(i) = vector3d(NaN, NaN, NaN);
            %             eV1(i) = vector3d(NaN, NaN, NaN);
            %             eV2(i) = vector3d(NaN, NaN, NaN);
            %             eV3(i) = vector3d(NaN, NaN, NaN);
            %             D(i) = rotation('Euler',0, 0, 0);
            %             T(i) = tensor([NaN NaN NaN; NaN NaN NaN; NaN NaN NaN]);
            %             centx(i)= e(j,k).x;
            %             centy(i) = e(j,k).y;
            %             mags(:,i) = [0; 0; 0];
            %             i = i+1;
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
    eCVA = e(id);
    eCVA.prop.CVA = CVAgrid;
    eCVA.prop.eV1 = eV1grid;
    eCVA.prop.eV2 = eV2grid;
    eCVA.prop.eV3 = eV3grid;
    eCVA.prop.D = Dgrid;
    eCVA.prop.T = Tgrid;
    eCVA.prop.mag1 = magsgrid(1,:);
    eCVA.prop.mag2 = magsgrid(2,:);
    eCVA.prop.mag3 = magsgrid(3,:);
    eCVA.prop.centx = centx;
    eCVA.prop.centy = centy;
    
    % remove any spurious results
    eCVA = eCVA(eCVA.mag1>0&eCVA.mag2>0&eCVA.mag3>0);
end

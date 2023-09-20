function [gNew] = newGrainsPGA(eNew)


%%
[gNew,eNew.grainId] = calcGrains(eNew('indexed'),'threshold',[1 10]*degree,'maxDist',2);

%%
gNew = gNew(phase);
pairs = gNew.neighbors;

% number of grains:
num=length(gNew);

% for keeping track of progress in for loop:
div=round(num/10);
count=div;

%%
eV = vector3d.nan(3,num);
mags = nan(3,num);
T = tensor.nan(num,'rank',2);
%% 
for k = 191291:num

    neighbs = unique([pairs(pairs(:,1) == k,2); pairs(pairs(:,2) == k,1)]);
    neighbs = neighbs(neighbs <= num);
    if length(neighbs)>2
    oTemp = gNew(neighbs).meanOrientation;
    [eV(:,k),mags(:,k),T(k)] = PGA(oTemp);
    end

    % Keep track of for loop progress and print to consoloe screen:
    perc=round(k/num*100);
        if k==count            
            fprintf('\n\n%i percent done...\n\n',perc)
            count=count+div;
        end


end



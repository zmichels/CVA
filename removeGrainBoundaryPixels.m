%% remove pixels at grain boundaries for intragranular analyses

function [eIntra] = removeBoundaryPixels(ebsd,grains)

%%
% check if grainID exists
if isempty(ebsd.grainId)
     
    error('There is no ebsd.grainId. Run calcGrains first.')

end

% get ebsd IDs at grain on either side of grain boundaries
eId = grains.boundary.ebsdId(:);

% remove any zeros
eId(eId==0)=[];

% remove the identified pixels from the dataset
ebsd('id',eId) = [];

% assign to new variable
eIntra = ebsd;
%%
end
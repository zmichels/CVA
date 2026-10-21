function matlab_mtex61()
% Bounded three-phase example using the maintained MATLAB CVA functions.
% MTEX 6.1.1 must already be installed on the MATLAB path.

repo = fileparts(fileparts(mfilename('fullpath')));
previousFolder = pwd;
restoreFolder = onCleanup(@() cd(previousFolder));
cd(tempdir); % Avoid root legacy functions taking precedence over addpath.
addpath(fullfile(repo,'matlab','mtex61'),'-begin');
assert(strcmp(which('PGA'),fullfile(repo,'matlab','mtex61','PGA.m')), ...
  'The maintained PGA.m is not first on the MATLAB path.');

ebsd = mtexdata('forsterite','silent');
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',10*degree);
selected = grains(grains.GOS > 0.5*degree);
selected = selected(1:min(24,length(selected)));

[grainResult,~] = grainsCVA(selected,ebsd, ...
  'method','batched','computeBulkVector',false);
[gridResult,~] = gridCVA(ebsd, ...
  'reference','mean','maxKernels',500,'computeBulkVector',false);

fprintf('CVA example: %d grains, %d windows, %d phases\n', ...
  length(grainResult),length(gridResult),numel(unique(grainResult.phaseId)));
fprintf('Leading grain CVA axis: [%.6f %.6f %.6f]\n', ...
  grainResult.CVA(1).x,grainResult.CVA(1).y,grainResult.CVA(1).z);
end

function result = matlab_mtex72_grainset(inputEbsd, makePlots)
%MATLAB_MTEX72_GRAINSET Two-pass grain CVA and full-eigenframe example.
% MTEX 7.2 must already be started. With no input, the example uses the
% small, multiphase Forsterite EBSD map supplied by mtexdata.
%
%   result = matlab_mtex72_grainset();
%   result = matlab_mtex72_grainset(myEbsd);
%   result = matlab_mtex72_grainset(myEbsd,false); % compute without plots
%
% This function can be called from any folder. For follow-up interactive
% CVA calls, work outside the repository root, which holds legacy functions
% with the same names. The maintained functions are in matlab/mtex61/ and
% have also been tested with MTEX 7.2.

if nargin < 2, makePlots = true; end
if nargin < 1, inputEbsd = []; end
assert(islogical(makePlots) && isscalar(makePlots), ...
  'makePlots must be a logical scalar.');
assert(startsWith(getMTEXpref('version'),'mtex-7.'), ...
  'Start MATLAB MTEX 7 before running this example.');

repo = fileparts(fileparts(mfilename('fullpath')));
maintained = fullfile(repo,'matlab','mtex61');
previousFolder = pwd;
restoreFolder = onCleanup(@() cd(previousFolder)); %#ok<NASGU>
cd(tempdir); % The repository root contains historical functions of the same name.
addpath(maintained,'-begin');
for name = {'PGA','grainsCVA','gridCVA'}
  expected = fullfile(maintained,[name{1} '.m']);
  assert(strcmp(which(name{1}),expected), ...
    'Expected maintained CVA function at %s; found %s.',expected,which(name{1}));
end

if isempty(inputEbsd), inputEbsd = mtexdata('small','silent'); end

%% Prepare the selected grainset without losing grain IDs
% Two passes reproduce the select-then-reconstruct workflow. The second
% calcGrains output is an EBSD object, not an array of grain IDs in MTEX 7.
ebsdCVA = inputEbsd('indexed');
[grains,ebsdCVA] = calcGrains(ebsdCVA,'threshold',[10 1]*degree);
grains = grains(grains.numPixel >= 3);
ebsdCVA = ebsdCVA(grains); % Use the EBSD object with grainId from calcGrains.
[grains,ebsdCVA] = calcGrains(ebsdCVA,'threshold',[10 1]*degree);

%% Grain and local 3-by-3 kernel analyses
% Disabling the optional bulk-axis KDE does not change the tensor or axes.
[gCVA,~] = grainsCVA(grains,ebsdCVA, ...
  'method','batched','computeBulkVector',false);
[eCVA,~] = gridCVA(ebsdCVA,grains, ...
  'reference','mean','computeBulkVector',false,'showProgress',false);

% CVA == eV1. eV2 and eV3 are the intermediate and least eigenvalue axes.
% Each result also carries mag1/2/3 and an unnormalized ODT in radians^2.
fprintf('MTEX 7 grainset example: %d reconstructed grains, %d grain CVAs, %d grid CVAs\n', ...
  length(grains),length(gCVA),length(eCVA));

figures = struct('singleGrain',[],'grainDensity',[],'gridDensity',[]);
if makePlots && ~isempty(gCVA)
  %% One grain: display its tensor and complete eigenframe
  g = gCVA(1);
  figures.singleGrain = figure('Name','One-grain ODT eigenframe','Color','w');
  plot(g.ODT,'antipodal','lower','figSize','small');
  hold on
  plot(g.eV1,'antipodal','lower','Marker','^', ...
    'MarkerFaceColor',[0.95 0.95 0.95],'MarkerEdgeColor','k','MarkerSize',12);
  plot(g.eV2,'antipodal','lower','Marker','s', ...
    'MarkerFaceColor',[0.95 0.95 0.95],'MarkerEdgeColor','k','MarkerSize',12);
  plot(g.eV3,'antipodal','lower','Marker','o', ...
    'MarkerFaceColor',[0.95 0.95 0.95],'MarkerEdgeColor','k','MarkerSize',12);
  mtexTitle(sprintf('Grain ID %d ODT eigenframe',g.id));

  %% Populations: one axis per grain versus one axis per accepted kernel
  figures.grainDensity = plotEigenframeDensity(gCVA,'Grain',length(gCVA));
  if ~isempty(eCVA)
    figures.gridDensity = plotEigenframeDensity(eCVA,'Grid kernel',length(eCVA));
  end
end

result = struct('grains',grains,'ebsd',ebsdCVA, ...
  'gCVA',gCVA,'eCVA',eCVA,'figures',figures);
end

function fig = plotEigenframeDensity(result,label,n)
% Three specimen-frame, antipodal lower-hemisphere density panels.
figure('Name',sprintf('%s eigenframe (%d axes)',label,n),'Color','w');
newMtexFigure('layout',[1 3],'figSize','large');
fig = gcf;
fig.Name = sprintf('%s eigenframe (%d axes)',label,n);
v = {result.eV1,result.eV2,result.eV3};
titles = {sprintf('eV1 (CVA), n = %d',n),'eV2 (intermediate)','eV3 (least)'};
for k = 1:3
  if k > 1, nextAxis; end
  plot(v{k},'antipodal','lower','smooth','halfwidth',10*degree);
  mtexTitle(titles{k});
end
setColorRange('equal');
mtexColorbar('title','M.U.D.');
end

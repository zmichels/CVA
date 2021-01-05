function [eV,mags] = PGA(ori)

%% P.G.A. – Principal Geodesic Analysis

%% References
% Zachary D. Michels, Seth C. Kruckenberg, Joshua R. Davis, and Basil Tikoff
% Determining vorticity axes from grain-scale dispersion of
% crystallographic orientations Geology, G36868.1, first published on July
% 17, 2015, doi:10.1130/G36868.1


%% ABOUT:
% This function performs a principal geodesic analysis on crystallographic
% orientations extracted from a single grain to identify a grain-scale
% vorticity vector associated with intragranular crystallographic
% dispersion. 

% input:
%                       ori     -   set of orientations

% output:
%                        eV     -   eigenvectors of PGA result
%                      mags     -   magnitudes of the eigenvectors
    

%% INTRODUCTION:
%   It has been observed that intragranular crystallographic dispersion
% axes can reflect the orientation of bulk vorticity in geologically
% deformed crystalline materials (Bestmann & Prior, 2003; Reddy & Buchan,
% 2005). We show that the coordinates of such rotational axes can be fit
% numerically for each grain in a sample and assembled to determine a
% preferred sample-scale vorticity axis.
%   The general approach in the following intragranular analysis can be
% applied to any set of concentrated orientation data. However, the
% following code is designed to accommodate electron backscatter
% diffraction data processed through the MTEX texture analysis toolbox for
% Matlab. Specifically, this function requires input of a single MTEX grain
% object (usually passed to the function from a larger MTEX GrainSet).
% Grain objects passed to this function must contain multiple orientation
% solutions for successful analysis.

%% 

% compute the projection of orientations onto the three dimensional
% tangential space centered about the mean
t = log(ori, mean(ori),'left');

% covariance matrix
cvm = t*t';


%%
% eigenvector analysis of the rotation covariance matrix:
[Vec,Diagonal] = eig(cvm);

% eigen-values : convert diagonal matrix to vector of eigen values:
value = diag(Diagonal);

% Sort eigen-values in descending order:
[val,index] = sort(value,'descend');

% sort and assign eigenvectors by corresponding eigenvalues
V1 = (Vec(:,index(1)));
V2 = (Vec(:,index(2)));
V3 = (Vec(:,index(3)));


% eigenvalue to axes magnitudes:
mags = val;

% sorted eigenvector directions:
dirs = [V1, V2, V3];

% sorted eigenvector vector3d objects for output:
eV = normalize(vector3d(mags'.*dirs));

end
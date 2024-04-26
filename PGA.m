function [eV,mags,T] = PGA(ori)

%%
ori = ori(:);


%% P.G.A. – Principal Geodesic Analysis

%% References
% Zachary D. Michels, Seth C. Kruckenberg, Joshua R. Davis, and Basil
% Tikoff Determining vorticity axes from grain-scale dispersion of
% crystallographic orientations Geology, G36868.1, first published on July
% 17, 2015, doi:10.1130/G36868.1

% Fletcher, P. Thomas, et al. "Principal geodesic analysis for the study of
% nonlinear statistics of shape." IEEE transactions on medical imaging 23.8
% (2004): 995-1005.


%% ABOUT:
% This function performs a principal geodesic analysis on crystallographic
% orientations extracted from a single grain to identify a grain-scale
% vorticity vector associated with intragranular crystallographic
% dispersion.

% input:
%                       ori     -   set of orientations

% output:
%                        eV     -   eigenvectors of dispersion tensor
%                      mags     -   magnitudes of the eigenvectors
%                         T     -   Dispersion tensor

% Prinicipal Geodesic Analysis is used to conduct a principal component
% analysis on sets of orientations/rotations. The result of the analysis
% is a covariance/dispersion matrix/tensor.


%%

oriRef = mean(ori);

% compute the projection of orientations onto the three dimensional
% tangential space centered about the mean
ver = extract(getMTEXpref('version'),digitsPattern);

if str2double(ver{2}) < 11
    t = log(ori, oriRef,'left');
else
    t = log(ori, oriRef,SO3TangentSpace.leftVector);
end

% covariance matrix as tensor
T = tensor(t*t,'rank',2);

% eigenvector of the tensor
[eV, mags] = eig(T);


%% sort largest first for compatibility with previous versions
[~,ind] = sort (mags,'descend');

eV = eV(ind);
mags = mags(ind);

end
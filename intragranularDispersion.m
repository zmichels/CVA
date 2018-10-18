function [dispersion_axis,D,eV1, eV2, eV3,mags,T]=intragranularDispersion(e)

%% References
% Zachary D. Michels, Seth C. Kruckenberg, Joshua R. Davis, and Basil Tikoff
% Determining vorticity axes from grain-scale dispersion of
% crystallographic orientations Geology, G36868.1, first published on July
% 17, 2015, doi:10.1130/G36868.1



% This function performs a principal geodesic analysis on crystallographic
% orientations extracted from a single grain to identify a grain-scale
% vorticity vector associated with intragranular crystallographic
% dispersion. This function is intended to be called within the function
% sample_scale_vorticity.m, which will perform the intragranular analysis
% on many grains from an MTEX GrainSet object to calculate the vector
% coordinates of a preferred bulk vorticity axis.

% input:
%                         e     -   EBSD data from a single grain from an
%                                   individual MTEX grain object comprising
%                                   multpile orientation solutions.

% output:
%           dispersion_axis     -   a vector3d object (MTEX v4.0 ) matching
%                                   the grain-scale crystallographic
%                                   dispersion/vorticity axis (CVA).
    

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


%% METHODS:
%   The following procudeure for Principal Geodesic Analysis includes...
%       - extraction of orientations from the grain object
%       - conversion of orientations to 3x3 rotation matrices
%       - conversion of 3x3 rotation matrices to (axis,angle) format
%       - calculation of rotation covariance matrix
%       - eigenvector analysis of the rotation covariance matrix


%% Extract orientations and mean orientation from grain object
% all intragranular orientations:
o=e.orientations;

%%
% mean orientation:
meano=mean(o);


%% Convert orientations to rotation matrices
% MTEX provides a function 'matrix' to convert orientations to rotation
% matrices. This produces a 3x3 matrix for each crystallographic
% orientation in which the columns of the matrix correspond to the primary
% direction of the crystallographic unit cell (equivalent to Miller Indices
% [100], [010], [001]). 

r=(matrix(o));

% The format of our analysis requires that the unit cell vectors be in the
% rows of the matrix rather than columns. So we transpose them here.
for i=1:size(r,3)
    r(:,:,i)=r(:,:,i)';
end

% Convert the mean orientation to rotation matrix and transpose. Uses MTEX
% function 'matrix' to convert quaternion to direction cosine matrix.
rotMean=matrix(meano)';


%% Covariance matrix, approximated in the tangent space at a given rotation
% Appropriate only if the sample is tightly clustered near the center,
% which is typically the case for intragranular crystallographic
% orientations. Requires a list of 3x3 rotation matrices (constructed
% above) and a 3x3 mean rotation matrix (calculated from the group of
% orientations above). The result of this section of code is to populate
% the covariance matrix variable 'cvm'.

L=length(r);

for i=1:L
        
    % compute statistical crossproduct of one rotation with the mean
    % rotation:
    cp=(r(:,:,i)'*rotMean);
    
    % convert to (axis, angle) format:
    cosine=(trace(cp)-1)/2;
    if cosine>=1
        a=0;
        u=[0 0 -1];      
        
    else
        u=sqrt(max([0;1+(cp(1,1)-1)/(1-cosine)]));
        
        if cp(3,2)<cp(2,3)
            u=-u;
        end
        
        v=sqrt(max([0;1+(cp(2,2)-1)/(1-cosine)]));
        
        if cp(1,3)<cp(3,1)
            v=-v;
        end
        
        w=sqrt(max([0;1+(cp(3,3)-1)/(1-cosine)]));
        
        if cp(2,1)<cp(1,2)
            w=-w;
        end
        
        nrm=sqrt(u^2 + v^2 + w^2);
%         nrm = 1;
        a=acos(max([-1;cosine]));
        u=[u/nrm v/nrm w/nrm];
        
    end
    
    % Matrix logarithm to produce infinitesimal from finite rotation:
    rL=a.*[0 u(3) -u(2); -u(3) 0 u(1); u(2) -u(1) 0];
    
    % Tangent vector from rotation:
    tv(i,:)=[rL(3,2) rL(1,3) rL(2,1)];
    
    % outer product of tangent vector with itself
    ms(:,:,i)=[tv(i,:)]'*[tv(i,:)];
end

% Initialize variable for covariance matrix
cvm = zeros(size(ms(:,:,1),1), size(ms(:,:,1),2));
sz=length(ms);
for i=1:sz
    cvm = cvm + ms(:,:,i);
end

% THIS IS THE COVARIANCE MATRIX
cvm = cvm / sz;
T = tensor(cvm);
%%

%% Rotation Principal Components (Principal Geodesic Analysis)

% eigenvector analysis of the rotation covariance matrix:
[Vec,Diagonal] = eig(cvm);

% eigen-values : convert diagonal matrix to vector of eigen values:
value = diag(Diagonal);

% Sort eigen-values in descending order:
[val,index] = sort(value,'descend');

% sort and assign eigenvectors by corresponding eigenvalues
V1 = real(Vec(:,index(1)));
V2 = real(Vec(:,index(2)));
V3 = real(Vec(:,index(3)));


% eigenvalue to axes magnitudes:
mags = val;



% eigenvectors:
dirs = [V1, V2, V3];
eV1 = vector3d(V1);
eV2 = vector3d(V2);
eV3 = vector3d(V3);

%% Select the 1st principal direction for the intragranular dispersion axis
dispersion_axis=vector3d((mags(1)*dirs(:,1))');

D=orientation('map',xvector,eV1,yvector,eV2);

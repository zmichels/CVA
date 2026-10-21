function [eV,mags,T] = PGA(ori,oriRef)
%PGA Specimen-frame principal geodesic analysis for MTEX 6.1.1 and 7.2.
% The unnormalised dispersion tensor is a sum in radians squared.
ori=ori(:);
if nargin<2 || isempty(oriRef), oriRef=mean(ori); end
t=log(ori,oriRef,SO3TangentSpace.leftVector);
xyz=matrix(t);
T=tensor(xyz*xyz','rank',2);
[eV,mags]=eig(T);
[mags,order]=sort(mags,'descend');
eV=eV(order);
end

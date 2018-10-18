function [vmax] = maxv(v,deg)

r = plotS2Grid('resolution',0.25*degree,'antipodal');
kde = calcDensity([v -v],r,'antipodal','halfwidth',deg*degree);
[~,I]=max(kde);

% get vector and negated vector (antipodal) of best-fit axis:
vmax=[r(I),-r(I)];

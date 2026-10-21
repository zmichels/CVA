function exportDensityFixture()
% Export a fixed-axis MTEX 6.1.1 bulk-vector reference for Python parity.
% Run after startup_mtex from the repository root.

axes = normalize(vector3d( ...
  [0.20, 0.25, 0.30, 0.23], ...
  [0.60, 0.58, 0.62, 0.56], ...
  [-0.77, -0.77, -0.72, -0.80]));
grid = plotS2Grid('resolution',0.25*degree,'antipodal');
density = calcDensity([axes,-axes],grid,'antipodal','halfwidth',10*degree);
[maximum,index] = max(density);
bulk = [grid(index),-grid(index)];
bulk(bulk.z>0) = [];

data = struct();
data.mtex_version = getMTEXpref('version');
data.axis_xyz = [axes.x(:), axes.y(:), axes.z(:)];
data.peak_xyz = [bulk.x, bulk.y, bulk.z];
data.peak_density = maximum;
data.grid_resolution_degrees = 0.25;
data.halfwidth_degrees = 10;
data.antipodal = true;

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outputPath = fullfile(root,'tests','fixtures','density_expected.json');
file = fopen(outputPath,'w');
assert(file>0,'Could not open %s',outputPath);
cleanup = onCleanup(@() fclose(file));
fprintf(file,'%s\n',jsonencode(data,'PrettyPrint',true));
fprintf('DENSITY_FIXTURE_PEAK=(%.9g,%.9g,%.9g) MAX=%.9g\n', ...
  bulk.x,bulk.y,bulk.z,maximum);
end

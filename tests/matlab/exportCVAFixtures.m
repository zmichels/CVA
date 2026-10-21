function exportCVAFixtures(destination)
%EXPORTCVAFIXTURES Regenerate fixed reference data with root legacy PGA.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin<1, destination=fullfile(root,'tests','fixtures'); end
oldpath=path; oldpwd=pwd; cleanup=onCleanup(@()restore(oldpath,oldpwd));
cd(tempdir); addpath(root,'-begin'); clear PGA;
rows=[]; cases=struct([]);
for k=1:4
 if k==2, sym='222'; else, sym='1'; end
 base=orientation.byEuler(37*degree,63*degree,21*degree,crystalSymmetry(sym));
 v=normalize(vector3d([1 0 0 1 2 -1],[0 1 0 2 -1 1],[0 0 1 -1 1 2]));
 a=[-3 2 1 4 -2 0.5]*degree;
 names={'rotated_triclinic','rotated_orthorhombic','near_zero','near_degenerate'};
 if k==3, a=a*1e-7; end
 if k==4
  v=vector3d([1 -1 0 0 0 0],[0 0 1 -1 0 0],[0 0 0 0 1 -1]);
  a=ones(1,6)*degree;
 end
 ori=rotation.byAxisAngle(v,a).*base; ref=mean(ori);
 [e,m,T]=PGA(ori); E=matrix(e);
 for j=1:3
  [~,i]=max(abs(E(:,j))); if E(i,j)<0, E(:,j)=-E(:,j); end
 end
 C=struct('name',names{k},'crystal_symmetry',sym,'mean',matrix(ref),...
 'tensor',matrix(T),'eigenvalues',m(:).','eigenvectors',E,...
 'tangent_vectors',matrix(log(ori,ref,SO3TangentSpace.leftVector)));
 if k==1, cases=C; else, cases(k)=C; end
 M=matrix(ori);
 for j=1:length(ori), rows(end+1,:)=[k reshape(M(:,:,j).',1,9)]; end
end
metadata=struct('reference','root PGA.m, MATLAB R2025b / MTEX 6.1.1',...
 'units','tangents radians; tensor and eigenvalues radians squared',...
 'tangent_frame','specimen, SO3TangentSpace.leftVector','normalization','sum, not mean',...
 'csv_layout','case_id,m11,m12,m13,m21,m22,m23,m31,m32,m33',...
 'specimen_symmetry','1','axis_convention','columns; axial; largest absolute component positive');
writematrix(rows,fullfile(destination,'pga_rotations.csv'));
fid=fopen(fullfile(destination,'pga_expected.json'),'w');
fprintf(fid,'%s\n',jsonencode(struct('metadata',metadata,'cases',cases),PrettyPrint=true)); fclose(fid);
exportGrains(destination);
end
function restore(p,d)
path(p); cd(d); clear PGA;
end

function exportGrains(destination)
ebsd=mtexdata('forsterite');
[grains,ebsd.grainId]=calcGrains(ebsd('indexed'),'angle',10*degree);
g=grains(grains.numPixel>=3 & grains.GOS>0.5*degree);
ids=[];
for phase=unique(g.phaseId).'
 gp=g(g.phaseId==phase); ids=[ids;gp.id(1:min(8,length(gp)))];
end
g=g('id',ids); ebsd=ebsd(g); samples=struct([]);
for i=1:length(g)
 ori=ebsd(ebsd.grainId==g.id(i)).orientations;
 [v,m,T]=PGA(ori);
 c=struct('grain_id',g.id(i),'phase_id',g.phaseId(i),...
 'crystal_symmetry',char(ori.CS.pointGroup),'mean',matrix(mean(ori)),...
 'rotations',permute(matrix(ori),[3 1 2]),'tensor',matrix(T),...
 'eigenvalues',m(:).','eigenvectors',matrix(v));
 if i==1, samples=c; else, samples(i)=c; end
end
metadata=struct('source','MTEX 6.1.1 forsterite bundled dataset',...
 'selection','calcGrains indexed, 10 degree threshold; numPixel>=3 and GOS>0.5 degrees; first 8 grain IDs per phase',...
 'rotation_layout','N by 3 by 3','reference','root legacy PGA',...
 'units','radians squared for tensors and eigenvalues','tangent_frame','specimen leftVector');
fid=fopen(fullfile(destination,'grains_expected.json'),'w');
fprintf(fid,'%s\n',jsonencode(struct('metadata',metadata,'grains',samples),PrettyPrint=true)); fclose(fid);
end

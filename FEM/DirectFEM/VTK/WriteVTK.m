    function [ ] = WriteVTK( C,E,Set,S,Sv,x,u,uExp,utot,t)
%% WriteVTK writes resuls in vtk format
% INPUT:
% C(e,:) connectivity of hexahedral element e
% E(e,:) = strains of element e: Exx,Eyy, Eyy, Exy,Exy,Exz,Eyz
% S(e,:) = Elastic stress of element e: Sxx,Syy, Syy,Sxy,Sxy,Sxz,Syz
% Sv(e,:) = viscous stress of element e: Sxx,Syy, Syy,Sxy,Sxy,Sxz,Syz
% x(i,:) : coordinates of node i
% u(i,:) : displacement of node i
% t      : time step
% OUTPUT:
% Creates folder OutputVTK with VRK files in it
%
CD=pwd();
if isfield(Set,'Data')
    NewFolder=strcat('OutputVTK_',Set.Data);
else
    NewFolder=strcat('OutputVTK');
end
FileName='MeshDisp';
if exist(NewFolder,'dir')
    cd(NewFolder);
    if t==0 % Only delete if no time is specified, since files for previous times are preserved
        delete('*.*');
    end
    cd(CD);
else
    mkdir(NewFolder);
end
FileName=strcat(NewFolder,Esc,FileName,sprintf('%i',t),'.vtk');
VTKfid=fopen(FileName,'w');
WriteVTKMesh(C,VTKfid,x);
First=true; % First set of results
WriteVTKResNodal('Displacements',VTKfid,u,First);
WriteVTKResNodal('DispExp',VTKfid,uExp,~First);
WriteVTKResNodal('DispTot',VTKfid,utot,~First);
% Compute averaged strains/stresses
nelem=size(S,1);
nnod=size(C,2);
nstrs=[1 1 3 3 3 3 3 6];
Ea=zeros(nelem,nstrs(nnod));
Sa=Ea;
Sav=Ea;
for e=1:nelem
    Ea(e,:)=sum(E(e,:,:))/size(E,2); % Average strain
    Sa(e,:)=sum(S(e,:,:))/size(S,2); % Average stress
    Sav(e,:)=sum(Sv(e,:,:))/size(S,2); % Average stress
end
WriteVTKResElemental('Strains',VTKfid,Ea,First);
WriteVTKResElemental('StressE',VTKfid,Sa,~First);
WriteVTKResElemental('StressV',VTKfid,Sav,~First);
if strcmp(Set.Rheology(1:6),'Kelvin')==1
    WriteVTKResElemental('StressTot',VTKfid,Sav+Sa,~First);
else
    WriteVTKResElemental('StressTot',VTKfid,Sa,~First);
end

fclose(VTKfid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write Results to VTK files
function WriteVTKResNodal(s,VTKfid,Res,First)
% Write in VTKFile the mesh part of the results
% INPUT:
% s = text with results name
% Res=list nof nodal reslts
% PointData=true, write "POINT_DATA" (first time results for points are written)
npoints=size(Res,1);
ncomp=size(Res,2);
Vector=ncomp>1;
if ncomp<3
    Res=[Res zeros(npoints,3-ncomp)];
end
% Write results
if First % First results
    fprintf(VTKfid,'POINT_DATA %i\n',npoints);
end
if Vector
    fprintf(VTKfid,'VECTORS  %s float\n',s);
    for n=1:npoints
        fprintf(VTKfid,'%e ',Res(n,:));
        fprintf(VTKfid,'\n');
    end
else
    fprintf(VTKfid,'SCALARS  %s float 1\n',s);
    fprintf(VTKfid,'LOOKUP_TABLE default\n');
    for n=1:npoints
        fprintf(VTKfid,'%e \n',Res(n));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteVTKResElemental(s,VTKfid,Res,First)
% Write in VTKFile the mesh part of the results
% INPUT:
% s = text with results name
% Res=list nof nodal reslts
% PointData=true, write "POINT_DATA" (first time results for points are written)
nelem=size(Res,1);
ncomp=size(Res,2);
Tensor=ncomp>1;
if ncomp<6
    Res=[Res(:,1:2) zeros(nelem,1) Res(:,3) zeros(nelem,2)];
end
% Write results
if First % First results
    fprintf(VTKfid,'CELL_DATA %i\n',nelem);
end
if Tensor
    ind=[1 4 5
        4 2 6
        5 6 3];
    fprintf(VTKfid,'TENSORS  %s float\n',s);
    for n=1:nelem
        for i=1:3
            fprintf(VTKfid,'%e ',Res(n,ind(i,:)));
            fprintf(VTKfid,'\n');
        end
        fprintf(VTKfid,'\n');
    end
else
    fprintf(VTKfid,'SCALARS  %s float 1\n',s);
    fprintf(VTKfid,'LOOKUP_TABLE default\n');
    for n=1:npoints
        fprintf(VTKfid,'%e \n',Res(n));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write Mesh to VTK files
function WriteVTKMesh(C,VTKfid,X)
% Write in VTKFile the mesh part of the results
% INPUT:
% X = nodal coordinates
% C = elemental connectivity
npoints=size(X,1);
nnodes=size(C,2);
if nnodes==8
    type=12;
elseif nnodes==4
    type=9;
else
    warning('Wrong numberof nodes %i',nnodes)
    type=0;
end
ncells=size(C,1);
dim=size(X,2);
if dim==2
    X=[X zeros(npoints,1)];
end
% Write coordinates
fprintf(VTKfid,'# vtk DataFile Version 3.98\n');
fprintf(VTKfid,'CNS_Direct_vtk\n');
fprintf(VTKfid,'ASCII\n');
fprintf(VTKfid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(VTKfid,'POINTS %i double\n',npoints);
for n=1:npoints
    fprintf(VTKfid,'%e %e %e \n',X(n,:));
end
% Write connectivity
fprintf(VTKfid,'CELLS %i %i\n',ncells,ncells*(nnodes+1));
for e=1:ncells
    fprintf(VTKfid,'%i ',nnodes,C(e,:)-1);
    fprintf(VTKfid,'\n');
end
% Write type of element
fprintf(VTKfid,'CELL_TYPES %i\n',ncells);
for e=1:ncells
    fprintf(VTKfid,'%i\n',type);
end
end

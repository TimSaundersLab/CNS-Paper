Pwd=pwd;
esc=Esc();
if ~exist(strcat(Pwd,esc,'Main.m'),'file')
    error('Please run program from folder with file Main.m in it');
end
addpath(strcat(Pwd,esc,'DirectFEM',Esc,'FEM'));
addpath(strcat(Pwd,esc,'DirectFEM',Esc,'Geometry'));
addpath(strcat(Pwd,esc,'DirectFEM',Esc,'VTK'));

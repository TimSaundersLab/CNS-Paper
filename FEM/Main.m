% Code for Stress Analysis of CNS from 2D velocity data
% Reads data in Set.DispFile. Datafile *.mat in the folder contains: *File
% inside folder contains:
%   XEul(:,i) = list of x positions in pixels at time i where velocity is read (Eulerian coorindates, fixed in space)
%   YEul(:,i) = list of y positions in pixels at time i where velocity is read (Eulerian coorindates, fixed in space)
%   UEul(:,i) = list of x velocities in pixel/time-step at time i (Eulerian coorindates, fixed in space)
%   VEul(:,i) = list of y velocities in pixel/time-step at time i (Eulerian coorindates, fixed in space)
%   Time: stime-steps
%   xref: position of reference point in um. PIV and Mesh are corrected by -xref for having the same 0 value
%   t0= time of initial time-step in minutes.
%   dt= time-step size in minutes.
% pixels=size of pixel in um.
% xref=reference point for analysis
% INSTRUCTIONS:
% Set in DispFile the file to read the velocity field.
% Run Main.m
% This runs DirectFEM\MainFEM with a forward FEM analysis with rheology
% model in Set.Rheology.
clearvars;
Set.MeshFile='CNS_Curved100x_3D.dat'; % Used in 3D: '1element.dat'; % Mesh file crated with GID
Set.Long=3; % =3: Longest direction in datafile is z. This coordinate is set as first (x) in CNS_Direct.
Set.CenterY=true; % true=correct location of PIV points such that they are centred at Y=0
Set.Dimensions=3; % Select 2 or 3 dimensions
Set.dx=2; % Size of element in um. Used in 2D. Options are 1, 2, 4 (lower ersolution)
% Sham
File=strcat('ShamMatFiles','results_velocity3');
Set.DispFile=strcat('ShamMatFiles',Esc,'results_velocity3'); % 'results_1elementIncremental', 'results_1elementNIncremental';
Set.ConstrainY=false; % do not constrain long Y direction, just minimal direction
% Tim files
Data={'WTZIP1','WTZIP2','ELAVZIP2','ELAVZIP4','REPOZIP1','REPOZIP4','TUBULIN','REPOGRIM'};
Set.ts=410; % Time in minues where analysis will start, and skipping previous increments. If =0, starts from the first increment.
            % Files writen in VTK are numbered 1,2,.. from this time onwards.
xmax=192; % Maximum length of analysis. End coordinate is resultats.xref+xmin. If NaN, maximum length is whole length             
for i=1:length(Data)
    Set.Data=Data{i};
    Set.DispFile=strcat('TimMatFiles',Esc,Set.Data,Esc,'resultats'); % WTZIP1, ELAVZIP2, REPOZIP1
    Set.ConstrainY=true; % do not constrain long Y direction, just minimal direction
    load(Set.DispFile);
    Set.ReductionX=resultats.pixels; % Convert from pixesl to um. Factor being applied to (Xlag,Ylag) coordinates and (Ulag,Vlag).
    Set.TransX=-resultats.xref; % Translation in um of PIV and Mesh (the latter after trimming). Translates this amount resultats.Xlag after scaling.
    Set.TrimPIV=[NaN NaN -25 25 ]; % trim data to [xmin xmax ymin ymax]. Triming is applied after centering if CenterY=true
    Set.AccumulateU=true; % Experiental applied displacements are accumulated in elastic problems.
    Set.PIVtol=0.75; % PIV displacment applied on nodes closer to dx(i)*PIVtol from measured position
    xminmax=[min(resultats.XEul(:,1)) max(resultats.XEul(:,1))]*resultats.pixels;
    Set.TrimMesh=[resultats.xref min(xmax+resultats.xref,floor(0.9*xminmax(2))) NaN NaN NaN NaN]; % trim mesh to [xmin xmax ymin ymax zmin zmax]. [80 NaN NaN NaN NaN NaN] Remove elements/nodes on first 80 units.
    Set.E=75; % in Pa. Default 1000
    Set.v=0.3; % Default 0.
    Set.eta=7*Set.E; % Charactersitic time ~7. Viscous coefficient
    Set.Rheology='Maxwell'; % Options: Kelvin, Maxwell, or elastic
    Set.MatFile=strcat('KimoCNS_',Set.Data,'.mat');
    Set.WriteVTK=true;
    Set.WriteVTKdt=1;
    if isfield(resultats,'dt')
        Set.WriteVTKdt=ceil(10/resultats.dt); % Increment in time-steps
    end
    % Call to MainFEM
    addpath(strcat(pwd,Esc,'DirectFEM'));
    tic;
    if isfield(resultats,'t0')
        Set.t0=resultats.t0;
    end
    ErrorFlag=MainFEM(Set);
    if(ErrorFlag==0)
        disp('Program Main of CNS_Direct succesfully finished.')
    else
        warning('Error encountered in DirectFEM.')
    end
    toc;
end

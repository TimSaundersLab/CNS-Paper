function Set=SetDefaults(Set)
%% Set defulats values for variables in Main.m
% Time integration
if ~isfield(Set,'AccumulateU')
    Set.AccumulateU=true; %Displacements and stresses are accumulated at each 
                          % time step (assumed true in viscous problems). Otherwise problem solved from previos 
                          % time-step with displaced corrdinates as 
                          % reference coordinates.
end
if ~isfield(Set,'Dimensions')
    Set.Dimensions=2; % time integration parameter for transient problems
end
if ~isfield(Set,'theta')
    Set.theta=0.5; % time integration parameter for transient problems
end
if ~isfield(Set,'dt')
    Set.dt=1; % time integration parameter for transient problems
end
if ~isfield(Set,'t0')
    Set.t0=0; % Reference time for PIV that translates time results. Analysis starts if at incrememtn t: t*Set.dt+Set.t0>Set.ts
end
% Mesh and PIV
if ~isfield(Set,'TrimPIV')
    Set.TrimPIV=nan*ones(1,4); %[xmin xmax ymin ymax] coordinates of trimmming box for PIV data. NaN=no trimming. Triming is applied after centering if CenterY=true
end
if ~isfield(Set,'TrimMesh')
    Set.TrimMesh=nan*ones(1,6); %[xmin xmax ymin ymax zmin zmax] coordinates of trimmming box for mesh. NaN=no trimming.
end
if ~isfield(Set,'PIVtol')
    Set.PIVtol=0.5; % PIV displacment applied on nodes closer to dx(i)*PIVtol from measured position
end
if ~isfield(Set,'MeshFile')
    Set.MeshFile='CNS_Straight.dat';
end
if ~isfield(Set,'Long')
    Set.Long=3; % Longest axis (A-P) is z (3rd coordinate).
end
if ~isfield(Set,'CenterY')
    Set.CenterY=true; % aling center on Y=0
end
if ~isfield(Set,'ConstrainY')
    Set.ConstrainY=true; % do not constraina long Y direction, just minimal direction
end
if ~isfield(Set,'DispFile')
    Set.MeshFile='results_velcotityt1.mat';
end
if ~isfield(Set,'ReductionX')
    Set.ReductionX=1; % No reduction
end
if ~isfield(Set,'TransX')
    Set.TransX=0; % No X translation of data in data Xlag
end
% Material properties 
if ~isfield(Set,'E')
    Set.E=1000;
end
if ~isfield(Set,'v')
    Set.v=0.3;
end
if ~isfield(Set,'eta')
    Set.eta=1; % Vscous coefficient
end
if ~isfield(Set,'Rheology')
    Set.Rheology='Maxwell'; % Elastic, Kelvin or Maxwell
end
if ~isfield(Set,'Rheology')
    Set.Lagrangian=true; % PIV read in Lagrangian coordinates
end
if ~isfield(Set,'WriteVTK')
    Set.WriteVTK=false; % =true: write VTK files
end
if ~isfield(Set,'WriteVTKdt')
    Set.WriteVTKdt=10; % Increment in minutes for writing VTK
end
% TESTS
if strcmp(Set.Rheology(1:6),'Elasti')==0 && strcmp(Set.Rheology(1:6),'Maxwel')==0 && strcmp(Set.Rheology(1:6),'Kelvin')==0
    warning('Rheology %s not recognised. Assumed elastic: Set.Rheology=Elastic',Set.Rheology);
    Set.Rheology='Elastic';
end
if strcmp(Set.Rheology(1:6),'Elasti')==1 && Set.eta>0
    warning('Elastic problem with eta>0. Value of viscous coefficient ignored and set to 0.');
    Set.eta=0;
end
if strcmp(Set.Rheology(1:6),'Elasti')==0 && ~Set.AccumulateU
    warning('In viscoelastic problems, displacements are always accumulated. Chaning parameter to Set.AccumulateU=true.');
    Set.AccumulateU=true;
end
if strcmp(Set.Rheology(1:6),'Elasti')==0 && Set.eta<eps
    warning('Rhology is set as viscoelastic but Set.eta =0. Rheology set to Elastic.');
    Set.Rheology='Elastic';
end
 
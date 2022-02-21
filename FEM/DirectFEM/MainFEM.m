function ErrorFlag=MainFEM(Set)
% Computes FEM analsyis, writes VTK files and builds kimographs database.
% For plotting kymographs, see function: 
%
% PlotFromKymo(dx,Ekimo,Skimo,Xkimo,Xkimo0,Velo)
%
ErrorFlag=0;
Set=SetDefaults(Set);
SetPaths();
fprintf('PIV File: %s\n',Set.DispFile);
fprintf('Reading Mesh\n');
[C,dx,X,Error]=ReadMesh(Set);
if ~isempty(Set.TrimMesh)
    if ~isnan(max(Set.TrimMesh))
        [C,X]=TrimMesh(C,Set.TrimMesh,X);
    end
end
if Error
    error('Error reading mesh from file %s\n',Set.MeshFile)
end
fprintf('Mesh Data (before translation):\n')
fprintf('[Xmin, Xmax]=[%f %f]\n',min(X(:,1)), max(X(:,1)))
fprintf('[Ymin, Ymax]=[%f %f]\n',min(X(:,2)), max(X(:,2)))
if abs(Set.TransX)>eps
        X(:,1)=X(:,1)+Set.TransX;
        fprintf('Mesh Data (after translation due to reference point):\n')
        fprintf('[Xmin, Xmax]=[%f %f]\n',min(X(:,1)), max(X(:,1)))
        fprintf('[Ymin, Ymax]=[%f %f]\n',min(X(:,2)), max(X(:,2)))
else
    fprintf('Mesh Data has not been translated.\n')
end
X0=X;
load(Set.DispFile);
% Preprocess. Convert Tim's data to Sham's struct
if ~isfield(resultats,'Xlag')
    if ~isfield(resultats,'XEul')
        error('Data Xlag and XEul not found.\n')
    end
    resultats.Xlag=resultats.XEul;
    resultats.Ylag=resultats.YEul;
    resultats.Ulag=resultats.UEul;
    resultats.Vlag=resultats.VEul;
    resultats.time=resultats.Time;
    Set.Lagrangian=false;
end
fprintf('PIV Data:\n')
fprintf('[Xmin, Xmax]=[%i %i]\n',min(min(resultats.Xlag)), max(max(resultats.Xlag)))
fprintf('[Ymin, Ymax]=[%i %i]\n',min(min(resultats.Ylag)), max(max(resultats.Ylag)))
fprintf('[Umin, Umax]=[%f %f]\n',min(min(resultats.Ulag)), max(max(resultats.Ulag)))
fprintf('[Vmin, Vmax]=[%f %f]\n',min(min(resultats.Vlag)), max(max(resultats.Vlag)))
% Contains data restultats struct with fields:
% Xlag(npoints, ntimes): x coordinates
% Ylag(npoints, ntimes): y coordinates
% Ulag(npoints, ntimes): x velocities
% Vlag(npoints, ntimes): y velocities
% time(ntimes): time values
Times=resultats.time;
nt=length(Times(Times>0));
Times=Times(Times>0);
Dt=zeros(nt,1); % Time increments
Dt(1:end-1)=(Times(2:end)-Times(1:end-1))*resultats.dt;
Dt(end)=Dt(end-1);
utot=zeros(size(X));
nnod=size(C,2);
nstrs=[1 1 3 3 3 3 3 6];
En=zeros(size(C,1),nnod,nstrs(nnod));
Sn=En;
Snv=En;
if Set.CenterY
    ymax=max(resultats.Ylag(:,1)*Set.ReductionX);
    ymin=min(resultats.Ylag(:,1)*Set.ReductionX);
    ymean=(ymax+ymin)/2; % Center to 0 y coordinate
end
if Set.WriteVTK
    tVTK=0;
    tVTKn=0;
    WriteVTK(C,En,Set,Sn,Snv,X0,utot,utot,utot,tVTK);
end
nx=floor((max(X(:,1))-min(X(:,1)))/dx(1));
Xmax0=max(X(:,1));
Xmin0=min(X(:,1));
Ekimo=nan(nt,nx);
Skimo=nan(nt,nx);
Xkimo=nan(nt,nx);
Xkimo0=nan(nt,nx);
ts=0;
if isfield(Set,'ts')
    ts=Set.ts;
end
for t=1:nt
    Set.dt=Dt(t);
    if t*Set.dt+Set.t0<=ts
        continue
    end
    Set.t=Times(t); % Time-step
    PIV.X=[resultats.Xlag(:,t), resultats.Ylag(:,t)];
    PIV.U=[resultats.Ulag(:,t), resultats.Vlag(:,t)]*resultats.pixels; % * No velocitites (no *Set.dt), Displacement in um
    PIV.U(isnan(PIV.X(:,1)),:)=[]; % Remove data that has nans in X
    PIV.X(isnan(PIV.X(:,1)),:)=[];
    PIV.X(isnan(PIV.U(:,1)),:)=[]; % Remove data that has nans in U
    PIV.U(isnan(PIV.U(:,1)),:)=[];
    if abs(Set.ReductionX)>eps
        PIV.X=PIV.X*Set.ReductionX;
        PIV.U=PIV.U*Set.ReductionX;
    end
    if abs(Set.TransX)>eps
        PIV.X(:,1)=PIV.X(:,1)+Set.TransX;
    end
    if Set.CenterY
        PIV.X(:,2)=PIV.X(:,2)-ymean; % Center to 0 y coordinate
    end
    % Apply trimming of PIV
    if ~isnan(max(Set.TrimPIV))
        if ~isnan(Set.TrimPIV(1))
            PIV.U(PIV.X(:,1)<Set.TrimPIV(1),:)=[];
            PIV.X(PIV.X(:,1)<Set.TrimPIV(1),:)=[];
        end
        if ~isnan(Set.TrimPIV(2))
            PIV.U(PIV.X(:,1)>Set.TrimPIV(2),:)=[];
            PIV.X(PIV.X(:,1)>Set.TrimPIV(2),:)=[];
        end
        if ~isnan(Set.TrimPIV(3))
            PIV.U(PIV.X(:,2)<Set.TrimPIV(3),:)=[];
            PIV.X(PIV.X(:,2)<Set.TrimPIV(3),:)=[];
        end
        if ~isnan(Set.TrimPIV(4))
            PIV.U(PIV.X(:,2)>Set.TrimPIV(4),:)=[];
            PIV.X(PIV.X(:,2)>Set.TrimPIV(4),:)=[];
        end
    end
    [BC,nApplied]=ApplyDispPIV(dx,PIV,Set,X,X0);
    % BC.u(BC.u(:,2)==2,:)=[];BC.u=[BC.u;1 2 0;2 2 0;4 3 0]; % Stress free stretching
    if nApplied==0
        warning('Time %i / %i. No nodes with applied displacement. No problem being solvedSolving.\n',t,nt);
        ui=zeros(size(X));
        uExp=ui;
    else
        fprintf('Solving time-step= %i / %i, time= %i. Applied %i exp. displ on %i/%i nodes\n',t,nt,Set.dt*t+Set.t0,size(PIV.X,1),nApplied,size(X,1))
        [aux,En,Sn,Snv,ui,uExp]=FEM(BC,C,En,Set,Sn,X,X0);
        fprintf('Constrained %i / %i dofs. \n',length(aux),numel(X))
    end
    if ~Set.AccumulateU % X_n equals current deformation, DX=0.
        utot=utot+ui;
    else % If accumulated, difference X=X_n, and DX=X-X0 is added when applying displacements
        utot=ui;
    end
    X=X0+utot;
    if Set.WriteVTK
        % WriteVTK(C,En,Set,Sn,Snv,X0,utot,utot,utot,0); % Deletes
        % all VTKs when called with t=0
        if t==tVTKn + floor((Set.ts-Set.t0)/Set.dt) + Set.WriteVTKdt
            tVTK=tVTK+1; % Filel number in VTK
            tVTKn=t-floor((Set.ts-Set.t0)/Set.dt);     % Previous time when VTK is written
            WriteVTK(C,En,Set,Sn,Snv,X,ui,uExp,utot,tVTK);
        end
    end
    [Ekimo(t,:),Skimo(t,:),Xkimo(t,:),Xkimo0(t,:)]=BuildKimo(C,dx(1),En,Ekimo(t,:),Sn,Snv,Skimo(t,:),X,X0,Xmax0,Xmin0,Xkimo(t,:),Xkimo0(t,:));
    save(Set.MatFile,'Ekimo','Skimo','Xkimo','Xkimo0','dx','Set')
end
end
function [Ekimo,Skimo,Xkimo,Xkimo0]=BuildKimo(C,dx,E,Ekimo,S,Sv,Skimo,X,X0,Xmax0,Xmin0,Xkimo,Xkimo0)
% Output of kimographs. 
% Ekimo = Strains for each current spatial position
% Skimo = Stresses for each current spatial position
% Xkimo = Positions for each current spatial position
% Positions:
% tmax=248;clim=[-1000 1500];cmap=[1 1 1;jet];imagesc(Xkimo0(1:tmax,:),clim);set(gca,'YDir','normal');colormap(cmap);colorbar;xlabel('A-P axis');ylabel('Time');title('x Pos');set(gca,'FontSize',24)
% Stresses:
% tmax=248;clim=[-1000 1500];cmap=[1 1 1;jet];imagesc(Skimo(1:tmax,:),clim);set(gca,'YDir','normal');colormap(cmap);colorbar;xlabel('A-P axis');ylabel('Time');title('\sigma_{xx}');set(gca,'FontSize',24)
% Strains:
% tmax=248;clim=[-1 1];cmap=[1 1 1;jet];imagesc(Ekimo(2:tmax,:)-Ekimo(1:tmax-1,:));set(gca,'YDir','normal');colormap(cmap);colorbar;xlabel('A-P axis');ylabel('Time');title('\epsilon_{xx}');set(gca,'FontSize',24)
nelem=size(C,1);
nstrs=size(E,3);
Ea=nan(nelem,nstrs);
Sa=nan(nelem,nstrs);
Sav=nan(nelem,nstrs);
nkimo=zeros(1,nelem);
nx0=max(0,floor(Xmin0/dx));
for e=1:nelem
    Ea(e,:)=sum(E(e,:,:))/size(E,2); % Average strain
    Sa(e,:)=sum(S(e,:,:))/size(S,2); % Average stress
    Sav(e,:)=sum(Sv(e,:,:))/size(S,2); % Average stress
    lnod=C(e,:);
    xe0=floor(sum(X0(lnod,1))/length(lnod)/dx)-nx0; % Mean position of element in units of dx, minus nx0. A-P axis for increasing x 
    xe=floor(sum(X(lnod,1))/length(lnod)/dx)-nx0; % Current mean position of element in units of dx, minus nx0. -P axis for increasing x 
    if xe<=length(Ekimo) && xe>0 && xe0<=length(Ekimo) && xe0>0
        Xkimo0(xe0)=(xe+nx0)*dx;
        Xkimo(xe)=(xe+nx0)*dx;
        if nkimo(xe)==0
            Skimo(xe)=Sa(e,1);
            Ekimo(xe)=Ea(e,1);
        else
            Skimo(xe)=Skimo(xe)+Sa(e,1);
            Ekimo(xe)=Ekimo(xe)+Ea(e,1);
        end
        nkimo(xe)=nkimo(xe)+1;
    end
end
Skimo(nkimo>0)=Skimo(nkimo>0)./nkimo(nkimo>0);
Ekimo(nkimo>0)=Ekimo(nkimo>0)./nkimo(nkimo>0);
end

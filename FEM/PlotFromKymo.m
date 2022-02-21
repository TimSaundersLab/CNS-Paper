function PlotFromKymo(dx,Ekimo,Skimo,Xkimo,Xkimo0,Velo)
% Plot curves from kymograph data
% 1. Stress profiles at differnet imes
% 2. Strains profiles at differnet imes
% 3. Condnesation profiles
% 4. Stress kymographs
% 5. Strains kymographs
% 6. Stress 3D data. Profiles at each time.
% USAGE Example:
%
% PlotFromKymo(dx,Ekimo,Skimo,Xkimo,Xkimo0)
% PlotFromKymo(dx,Ekimo,Skimo,Xkimo,Xkimo0,'V3')
%
% For directly plotting kymograph:
% Positions:
% Imax=248;caxis([-1000 1500]);cmap=[1 1 1;fire];pcolor(Xkimo0(1:Imax,:),clim);shading flat;colormap(cmap);colorbar;xlabel('A-P axis');ylabel('Time');title('x Pos');set(gca,'FontSize',24)
% PLOT Stress-Strain-Positions
% PIV and FEM analysis is done according to coordinates Post-Ant:
%  x=0: Posterior. Minimum value: x=80
%  x=400: Anterior. Maximum value: x=400.
% Kymo database has database according to Ant-Post:
%  x=80: Anterior (corresponding to x=114 from left)
%  x=313. Posterior (correspoding to x= 434 from left). Maximum value 313
% Therefore, Kymo coordinates must be incremented by 34 to get tru x
% coordinate from anterior (left)
%--------------------------------------------------------------------------
Plot=[0 0 0 1 1 0 0 1]; % Plot(i)=1, plot corresponding Figure, =0 Do not plot.
%                i=1: Stress profile
%                i=2: Strain profile
%                i=3: Positions
%                i=4: Strain kymograph
%                i=5: Stress kymograph
%                i=6: stress heights on time-space plane.
%                i=7: stress peaks on time-space plane.
%                i=8: Velocity kymograph from row data. Plot(8)=2, plot
%                also time evolution at two AP positiosn. Plot(8)=3=include also 3D kymograph
dataN={'WTZIP1','WTZIP2','ELAVZIP2','ELAVZIP4','REPOZIP1','REPOZIP4','TUBULIN','REPOGRIM'}; % '100KelvinE75eta75V3'; % Image 8
Folder='Results2D400ele'; % Results2D200 Results3Dclose all
Fire=true; % Colormap in kymographs. If false, 'jet' is used.
% Color limits
slim=[-4 6];       % [-4 6] Stress limits in Kymographs
elim=[-0.2 0.3]; % Orig: [-0.2 0.3]; Strain limits in Kymographs
vlim=[-0.13 0.03];  % Orig: [-1.3 0.3]; Velocitiy limits in velocity kymographs
if Fire
    slim=slim/2.5;
    elim=elim/2.5;
end
PositiveTimeUP=false; % if true, time in kyograph goes from down to up.
L=510; % Maximum AP length. Used when reversed and plotting stress peaks.
Imax0=690; % Maximum time interval in minutes in velocity kymogrpahs
Lmax=258; % MAximum length from T1 for velocity kymograph
Spline=1; % If >1, Spline msoothing of kymographs (not big change). Value is subdivision of kymograph for finer scales. Values may be 2,3,4,...Values 0 or 1 do not modify kymo resolution.
% Loop on datasets
k=0;
for ii=1:length(dataN)
    data=dataN{ii};
    if nargin==0
        Pwd=pwd;
        fileStress=strcat(Pwd,Esc,Folder,Esc,'KimoCNS_',data);
        fileVelocity=strcat(Pwd,Esc,'TimMatFiles/',data,'/resultats.mat');
        load(fileStress,'dx','Ekimo','Xkimo','Skimo','Xkimo0','Set');
        load(fileVelocity,'resultats');
        Velo=data;
    end
    I0=ceil((Set.ts-Set.t0)/Set.dt)+1; % Initial increment. ts=Global start time of analysis for all datasets [minutes], t0=initial time of this dataset [minutes]
    dt=Set.dt; % time-step size [minutes]
    Imax=I0+min(Set.t-I0,Imax0/Set.dt); % Last step being analysed. 
    Dt=floor(floor(Imax/10)*10/30)*10; % For creating intervals of stresses
    t=dt*[1 Dt 2*Dt 3*Dt]'; % Times where stresses and strains are plotted
    %t=dt*[90 140 190 240]'; % Times where stresses and strains are plotted
    t=[t;max(t)*100]; % Add 1 more position that will be removed
    nticks=5; % Number of labels in x and y axes of kymographs
    s=num2str(t,'  %5i\n');
    s(end,:)=[];
    t(end,:)=[];
    Set.LW=4;
    machine=computer;
    x0=0; % Initial position. Remove if should be retrieved from TrimMesh
    if isfield(Set,'TrimMesh') && ~exist('x0','var')
        if ~isnan(Set.TrimMesh(1))
            x0=Set.TrimMesh(1); % Value of initial x coordinate in output.
        end
        ReverseX=false; % True if x axis (AP) is reversed
    elseif exist('x0','var')
        ReverseX=false; % True if x axis (AP) is reversed
    else
        ReverseX=true; % True if x axis (AP) is reversed
    end
    xmin=x0; % Bounds for stress/strain profile
    xmax=max(Xkimo(I0,:)); % for stress/strain profile
    if strcmp(machine(1:3),'MAC')
        Set.FS=24;
    else
        Set.FS=18;
    end
    if nargin<6 && isempty(Velo)
        Velo='';
    end
    emax=size(Xkimo,2);
    emax=[1 ceil(emax/4+1) ceil(emax/2) ceil(emax*3/4) emax];
    nt=length(t);
    ne=length(emax);
    e=zeros(nt,ne);
    for i=1:nt
        if ReverseX
            x=L-Xkimo(t(i)/dt,:);
        else
            x=Xkimo(t(i)/dt,:);
        end
        if Plot(1)
            figure(1+k)
            plot(x,Skimo(t(i)/dt,:),'o-','LineWidth',Set.LW)
            set(gca,'fontSize',Set.FS)
            hold on;
        end
        if Plot(2)
            figure(2+k)
            plot(x,Ekimo(t(i)/dt,:),'o-','LineWidth',Set.LW)
            set(gca,'fontSize',Set.FS)
            hold on;
            e(i,:)=emax*dx(1);
            s(i,1:2)='t=';
        end
    end
    % 1. Stresses. Profile.
    if Plot(1)
        figure(1+k);
        hold off
        if ~isempty(Velo)
            title(Velo)
        end
        ylabel('\sigma_{xx}','FontSize',Set.FS);
        xlabel('x pos A-P [um]','FontSize',Set.FS);
        h=legend(s,'Location','NorthWest');
        h.Position(1)=h.Position(1)+0.4;
        set(h,'FontSize',Set.FS)
        if ReverseX
            Xk=L-Xkimo(t/dt,:);
        else
            Xk=Xkimo(t/dt,:);
        end
        Sk=Skimo(t/dt,:);
        ymin=min(Sk(Xk>xmin & Xk<xmax));
        ymax=max(Sk(Xk>xmin & Xk<xmax));
        axis([xmin xmax 1.1*ymin 1.1*ymax]); % ;-Inf Inf]); % E1E3, E=1E4=axis([80 325 -150 150]); % V3
    end
    % E=1E3: axis([80 325 -250 250]); % V3
    % 2. Strains. Profiles.
    if Plot(2)
        figure(2+k);
        hold off
        ylabel('\epsilon_{xx}','FontSize',Set.FS);
        xlabel('x pos A-P [m]','FontSize',Set.FS);
        if ~isempty(Velo)
            title(Velo)
        end
        h=legend(s,'Location','NorthWest');
        h.Position(1)=h.Position(1)+0.4;
        set(h,'FontSize',Set.FS)
        Ek=Ekimo(t/dt,:);
        ymin=min(Ek(Xk>xmin & Xk<xmax));
        ymax=max(Ek(Xk>xmin & Xk<xmax));
        axis([xmin xmax 1.1*ymin 1.1*ymax]);% -Inf Inf]); % axis([80 325 -0.5 0.5]); % V3
    end
    % 3. Displacements and condensation lines.
    if Plot(3)
        figure(3+k);
        t=1:size(Xkimo,1);
        nt=length(t);
        x=zeros(nt,ne);
        for i=1:nt
            if ReverseX
                x(i,:)=L-Xkimo0(t(i),emax);
            else
                x(i,:)=Xkimo0(t(i),emax);
            end
        end
        for i=1:size(e,2)
            plot(x(:,i),t*dt,'LineWidth',Set.LW);
            hold on;
        end
        hold off
        title(strcat('x-Position, ',Velo),'Interpreter','latex');
        xlabel('A-P axis','FontSize',Set.FS);
        ylabel('time [s]', 'FontSize',Set.FS);
        set(gca,'FontSize',Set.FS)
    end
    %4. Strain Kymograph
    if Plot(4)
        figure(4+k)
        clf;
        sKimo=Ekimo(I0+1:Imax,:)-Ekimo(I0:Imax-1,:);
        sKimo(end+1,:)=sKimo(end,:);
        if Spline>1
            x=size(sKimo);
            pp=csapi({1:x(1),1:x(2)},sKimo);
            x1=1:1/Spline:x(1); x2=1:1/Spline:x(2);
            sKimo=fnval(pp,{x1',x2'});
        end
        if Fire
            cmap=fire;
%            sKimo(isnan(sKimo))=elim(2);
        else
            colormap('jet');
            cmap=[colormap];
        end
        pcolor(sKimo);
        shading flat;
        if PositiveTimeUP
            set(gca,'YDir','normal');
        end
        colormap(cmap);
        colorbar;
        xlabel('A-P axis','FontSize',Set.FS);
        ylabel('Time','FontSize',Set.FS);
        title(strcat('Strain rates, $\dot{\epsilon}_{xx}$,',Velo),'Interpreter','latex');
        ndata=size(sKimo,2);
        xticklabels = x0:floor(ndata/(nticks-1)*dx(1)):x0+floor(ndata*dx(1));
        xticks = linspace(1,ndata,nticks); %linspace(1, ndata*dx(1), nticks*dx(1)); % Position
        nt=size(sKimo,1);
        yticklabels = 0:floor(nt/(nticks-1)*dt):floor(nt*dt); % values of y-labels
        yticklabels =round(yticklabels/20)*20;
        yticks = linspace(1, nt, nticks);               % index Position of y-tick labels 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
        caxis(elim)
        set(gca,'FontSize',Set.FS);
    end
    %Stress  Kymograph
    if Plot(5)
        figure(5+k)
        sKimo=Skimo(I0:Imax,:);
        if Fire
            cmap=[fire];
%c            sKimo(isnan(sKimo))=slim(2);        
        else
            colormap('jet');
            cmap=[colormap];
        end
        pcolor(sKimo);
        shading flat;
        %A caxis(0.25*[-1000 1500]); % eta=1E4, E=1E3, 0.25
        % eta=1E4, E=1E4, 0.5
        if PositiveTimeUP
            set(gca,'YDir','normal');
        end
        colormap(cmap);
        colorbar;
        xlabel('A-P axis');
        ylabel('Time');
        title(strcat('Stress, $\sigma_{xx}$, ',Velo),'Interpreter','latex');
        set(gca,'FontSize',Set.FS)
        if ~exist('xticks','var') || ~exist('yticks','var')
            ndata=size(sKimo,2);
            xticklabels = x0:floor(ndata/(nticks-1)*dx(1)):x0+floor(ndata*dx(1));
            xticks = linspace(1,ndata,nticks); %linspace(1, ndata*dx(1), nticks*dx(1)); % Position
            nt=size(sKimo,1);
            yticklabels = 0:floor(nt/(nticks-1)*dt):floor(nt*dt); % values of y-labels
            yticks = linspace(1, nt, nticks);               % index Position of y-tick labels
        end
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
        caxis(slim)
        set(gca,'FontSize',Set.FS);
    end
    % 3D cloud of stress profiles
    if Plot(6) && exist('pcshow','builtin')
        figure(6+k)
        nt=size(Skimo,1);
        nx=size(Skimo,2);
        xt=reshape(Xkimo,nx*nt,1);
        yt=kron(ones(nx,1),(1:nt)');
        zt=reshape(Skimo,nx*nt,1);
        x=xt;
        y=yt;
        z=zt;
        zmax=20;
        x(isnan(xt) | isnan(yt) | isnan(zt) | zt>zmax)=[];
        y(isnan(xt) | isnan(yt) | isnan(zt) | zt>zmax)=[];
        z(isnan(xt) | isnan(yt) | isnan(zt) | zt>zmax)=[];
        h=pcshow([x y z]);
        xlabel('AP axis [um]');
        ylabel('time [s]');
        zlabel('\sigma_{xx}');
        h.FontSize=Set.FS;
        h.LineWidth=Set.LW;
        %h.ZLim=[h.ZLim(1) 40];
    elseif Plot(6)
        warning('PCSHOW not available. No 3D plot of stress values')
    end
    if Plot(7) % Superimposed stress peaks and data
        %h.ZLim=[h.ZLim(1) 40];
        % Plot Least-Squares approximation of stress lines of stress minima
        % Settings:
        % Non-mutant. Image8 data.
        % x0=[135 155 210 240 275 310 340 375 400 420]; % Previous version. On top of lines.
        if sum(fileStress(end-1:end)=='V3')==2 
            x0v=[120 155 200 240 275 300 330 360 385 410];
        else
            x0v=[120 155 200 240 275 300 330 360 385 410]-119; % Initial coordinates of stress minima. Coordinates according to last plot.
        end
        Set.LeastSquares=true;
        Set.dt=dt; % Time-step length. Used in final plotting
        Set.Shift=0; % Shift for initial stress peaks
        Set.nLS=4; % For least-square fitting: degree of polynomial or exponential, exp(t^n)
        Set.Exponential=false; % If not exponential, use polynomail
        Set.tol=0; % Tolerance for computing stress minima. Minima in S(i) must satisfay S(i-1)>S(i)+tol<S(i+1)
        Set.PlotPeaks=false; % Plot dots of stress minima that are fitted
        Set.dxs=15; % Width of x average (horizontal), for retrivening mean lines in least squares
        Set.dts=5; % Width on t (vertical), for retrivening mean lines
        Set.x=[100 450]; % Anterior-Posterior interval [um] for triming images and output
        Set.y=[100 150]; % vetical (left to right) interval [um]
        Set.nt=[1 Set.t]; % Triming in time
        % Apply shoftingand offset
        x0v=x0v+Set.Shift;
        x0v=x0v-Set.x(1);
        Set.file=fileStress;
        Set.Figure=7+k;
        Set.SaveKymoData=false;
        [Kymograph,StressLines,StressMin]=StressPeaksOnKymo(Set,x0v);
        PlotStressPeaks(Kymograph,StressLines,StressMin)
    end
    if Plot(8) % Velocity kymograph
        % Compute dx
        x=resultats.XEul(:,1); % cooridnates in pixels
        xmin=min(x); % Minimum x in pixels
        x(x==xmin)=[];
        dx=(min(x)-xmin)*resultats.pixels; % Size of kymograph pixel x-size in um
        xref=resultats.xref; % Position of T1 in um
        xmin=xmin*resultats.pixels;  % Minimum coordinate in um
        xmax=max(x)*resultats.pixels;% MAximum coordinate in um
        % Kymograph on length form T1 (xref) to xref + Lmax
        nx=floor(Lmax/dx)+1; % Initial number of kymograph pixels
        nref1=max([1,(xmin-xref)/dx]); % Reference kymogrph x-position for left end point
        nref2=min((xmax-xref)/dx,floor(Lmax/dx)); % Reference kymograph x-position for right end point
        V=NaN(Imax-I0+1,nx);
        Ends=[nref1+1 nref2-1]; % Current position of endpoints in kyomograph xpoints (but real, not integer)
        for n=I0:Imax
            x=resultats.XEul(:,n)*resultats.pixels;
            v=resultats.UEul(:,n);
            v(isnan(x))=[];
            x(isnan(x))=[];
            x(isnan(v))=[];
            v(isnan(v))=[];
            x(x<resultats.xref)=[];
            v(x<resultats.xref)=[];
            x(x>Lmax+xref)=[];
            v(x>Lmax+xref)=[];
            nV=zeros(1,nx);
            for i=1:length(x)
                j=floor((x(i)-xref)/dx+1);
                if isnan(V(n-I0+1,j))
                    V(n-I0+1,j)=v(i);
                else
                    V(n-I0+1,j)=V(n-I0+1,j)+v(i);
                end
                nV(j)=nV(j)+1;
            end
            V(n-I0+1,:)=V(n-I0+1,:)./nV;
            % Trim wiht NaN accordin to end pixel velocities
            EndsNew(1)=max(2,Ends(1) + V(n-I0+1,floor(Ends(1)))*resultats.pixels/dx); % Velocity is in pixel/increment. Convert to kymograph-pixel/increment
            EndsNew(2)=min(nx-1,Ends(2) + V(n-I0+1,ceil(Ends(2)))*resultats.pixels/dx);
            V(n-I0+1,1:max(2,floor(Ends(1))))=NaN;
            V(n-I0+1,min(nx-1,ceil(Ends(2))):nx)=NaN;
            Ends=EndsNew;
        end
        figure(8+k)
        if Fire
            cmap=fire;
        else
            cmap=colormap('jet');
        end
        pcolor(V*resultats.pixels/resultats.dt)
        shading flat;
        if PositiveTimeUP
            set(gca,'YDir','normal');
        end
        colormap(cmap);
        colorbar;
        xlabel('A-P axis');
        ylabel('Time');
        title(strcat('Velocity [um/min], ',Velo),'Interpreter','latex');
        set(gca,'FontSize',Set.FS)
        xticklabels = x0:floor(nx/(nticks-1)*dx):x0+floor(nx*dx); % x-Labels at position xticks
        xticks = linspace(1,nx,nticks); % Position of xticklabels
        nt=size(V,1);
        yticklabels = Set.ts + (0:floor(nt/(nticks-1)*dt):floor(nt*dt)); % Position and values of y-labels
        yticks = linspace(1, nt, nticks);               % index Position of x-tick labels 
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
        set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
        set(gca,'FontSize',Set.FS);
        caxis(vlim)
        if Plot(8)>=2 && k==0
            figure(9)
            xlabel('Time [min]')
            ylabel('Velocity [um/min]')
            plot(1:nt,V(:,floor(xticks(2)))*resultats.pixels/resultats.dt,'LineWidth',Set.LW)
            hold on;
            plot(1:nt,V(:,floor(xticks(4)))*resultats.pixels/resultats.dt,'LineWidth',Set.LW)
            set(gca, 'XTick', yticks, 'XTickLabel', yticklabels)
            title(Velo,'Interpreter','latex');
            set(gca,'FontSize',Set.FS-2);
            legend(strcat('x=',num2str(xticklabels(2)),'um'),strcat('x=',num2str(xticklabels(3)),'um'),'Location','SouthEast');
        end
        if Plot(8)>=3 && k==0
            figure(10)
            [X,Y] = meshgrid(1:nt,1:nx);
            surf(X,Y,resultats.pixels/resultats.dt*V(:,end:-1:1)');
            xlabel('Time')
            ylabel('A-P')
            set(gca, 'XTick', yticks, 'XTickLabel', yticklabels)
            set(gca, 'YTick', xticks, 'YTickLabel', xticklabels(end:-1:1))
            title('Velocity surface')
            set(gca,'FontSize',Set.FS-2);
            colormap(cmap);
        end
    end
    k=k+sum(Plot);
end
end
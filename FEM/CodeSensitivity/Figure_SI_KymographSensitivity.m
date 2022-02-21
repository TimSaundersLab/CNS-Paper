function Figure_SI_KymographSensitivity()
% Analysis of visoelastic 1D domain with nelem elements
% Chooses apporach b) from following two possible extensions for 1 elements:
% a) Elemental equilibrium (internal viscosity): for each element [x_1^e x_2^e]
%    f=k1*(l^e-l0^e)+k2*(l^e-L^e)+eta*\dot l^e
%    x_1^e=x_2^{e-1}
%    x_2^e=x_1^e + l^e
% b) Nodal equilibrium (stres continuity, external viscosity):
%    l^e=x_2^e-x_1^e
%    sigma_1^e=-k1*(l^e-l0)-k2*(x2-x1-L)+eta*(x_1^{n+1}-x_1^n)/dt
%    sigma_2^e=k1*(l^e-l0)+k2*(x2-x1-L)+eta*(x_1^{n+1}-x_1^n)/dt
%    0=k1*(l^e-l0^e)+k2*(l^e-L^e)+eta*\dot l^e
%    x_1^e=x_2^{e-1}
%    x_2^e=x_1^e + l^e
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Geo.dl=-0.1; %-0.1; % If dl exist, imposed increment in current length. If =0, constatn value
%
% Material-Geometry Parameters
% Reference values
%   gamma k1   k2 eta delay   L0   l0  force
p0=[0.21 0.01 1.9 15  20     0.95 1.0 0.005];
par={'Remodelling rate \gamma','Stiffness k_1','Stiffness k_2','Viscosity \eta','delay \delta t'};
%Gamma: 0.6; 0.21 stable; 0.22 unstable
%eta: 15 stable , eta=16 unstable (free end),
%L0= Reference length per element
%l0= Initial apparent length per element
%f= Initial and final force at right end
%
Run=[0 1 1];
% Run(1)=1 Run single cases and show kymographsingle cases.
% Run(2)=2 Run sensisitivity for perturbation of each parameter alone. =1 just plot from Sensitivity.mat, not run
% Run(3)=2 Run condensation for combinations of 2 parameters: (k1,eta),(k1,tau). =2. no run, only plot for saved file 'Condensation'
% (tau,k1) sensisitivity for perturbation of each parameter alone.
% =1 or 2: sensisitivity analysis
% Time-integration and other settings
C=computer;
if C(1)=='M'
    FS=24; % Fontsize
else
    FS=18; % Fontsize
end
Set.LW=2; % Linewidth in sensitivity analysis
Set.gray=[0.5 0.5 0.5];
Geo.nelem=100; % number of elements
Geo.nx0=100; % Resolution in x: number of points per unit in x.
Geo.ny=1;% nt/tend; % Resolution on y. 1=maximum resolution (same as data).
Geo.ns=10; % Number of segments
Set.dt=1;
Set.tend=1250; % 200
% Single values kymographs Stable Unstable
if Run(1)>=1
    Set.Plot=1;
    % Stable values
    RunCase(Geo,p0,Set);
    % Single stable plots
    % Unstable due to gamma
    gm=p0(1);
    p0(1)=1.05*p0(1);
    RunCase(Geo,p0,Set);
    p0(1)=gm;
    % Unstable due to k2
    k2=p0(3);
    p0(3)=0.95*p0(3);
    RunCase(Geo,p0,Set);
    p0(3)=k2;
    % Unstable due to eta
    eta=p0(4);
    p0(4)=1.07*p0(4);
    RunCase(Geo,p0,Set);
    p0(4)=eta;
end
if Run(2)>=1
    % Sensitivity Plots
    Set.Plot=0;
    if Run(2)>1
        p=p0;
        np=5;
        % Factors that mulitply each nominal parameter value
        s=ones(np,1)*[0.9 0.94 0.98 1 1.03 1.07 1.08]; % gamma
        s(2,:)=[0.5 0.6 0.8  1 4 8 10]; % k1
        s(3,:)=[0.9 0.94 0.97 1 1.1 1.2 1.3]; % k2
        s(4,:)=[0.7 0.8 0.9 1 1.05 1.1 1.12]; % eta
        s(5,:)=[0.7 0.8 0.9 1 1.1 1.2 1.3];% delay
        nj=size(s,2);
        da=zeros(np,nj);
        da0=zeros(np,nj);
        dl=zeros(np,nj);
        fr=zeros(np,nj);
        for i=1:np
            for j=1:nj
                fprintf('Parameter %i/%i. Loop %i/%i\n',i,np,j,size(s,2))
                p(i)=s(i,j)*p0(i);
                [dl(i,j),da(i,j),da0(i,j),fr(i,j)]=RunCase(Geo,p,Set);
            end
            p=p0;
        end
        save('Sensitivity','da','da0','dl','fr','s','p0','Set');
    end
    if Run(2)>0
        load('Sensitivity','da','da0','dl','fr','s','p0','Set');
        np=size(da,1);
        for i=1:np
            % Plot sensitivity
            figure(i);clf
            x=s(i,:)*p0(i);
            y1=dl(i,:)+da(i,:)/2;
            y2=dl(i,:)-da(i,:)/2;
            patch([x x(end:-1:1)],[y1 y2(end:-1:1)],Set.gray)
            hold on
            plot(x,dl(i,:)+da0(i,:)/2,'k--') % Relative difference in amplitude of oscillation
            plot(x,dl(i,:),'LineWidth',Set.LW) % Relative difference in VNC length
            % plot(x,dl(i,:)-da0(i,:)/2,'k--') % Relative difference in amplitude of oscillation
            % yyaxis right
            % [AX,H1,H2]=plotyy([],[],x,10*fr(i,:),'r','LineWidth',Set.LW); % Frequency
            [AX,H1,H2]=plotyy(x,dl(i,:)-da0(i,:)/2,x,10*fr(i,:)); % Frequency
            H1.Annotation.LegendInformation.IconDisplayStyle = 'off';
            set(gca,'FontSize',FS)
            % yyaxis right
            ylabel(AX(2),'Frequency X 10','FontSize',FS)
            ylabel(AX(1),'Relative shortening ','FontSize',FS)
            xlabel(par{i},'FontSize',FS)
            legend('Final Osc. Amplitude ','Initial Osc. Amplitude','Relative VNC length','Frequency','Location','SouthWest')
            yminmax=[0 0.95 0.10 0.3]; % ylimits in length and freq
            for a=1:2
                AX(1).Position=[0.15 0.19 0.67 0.75];
                AX(a).FontSize=FS;
                ymin=yminmax((a-1)*2+1);
                ymax=yminmax((a-1)*2+2);
                AX(a).XLim=[min(x) max(x)];
                AX(a).YLim=[ymin ymax];
                AX(a).YTick=ymin:(ymax-ymin)/5:ymax;
                AX(a).YTickLabel={num2str(AX(a).YTick(1)),num2str(AX(a).YTick(2)),num2str(AX(a).YTick(3)),num2str(AX(a).YTick(4)),num2str(AX(a).YTick(5)),num2str(AX(a).YTick(6))};
                AX(a).YColor='k';
            end
            H1.LineStyle='--';
            H1.Color='k';
            H2.Color='g';
            H2.LineWidth=Set.LW;
            %f.Position=[360 357 739 565];
            hold off
            % axis([min(x) max(x) 0 0.95]);
        end
    end
end
if Run(3)>0 % Condensation rates
    Set.Plot=0;
    if Run(3)>1
        s=zeros(5,10);
        s(1,:)=[0.8  0.84 0.88 0.92 0.96 1    1.02 1.04 1.06 1.08]; % gamma
        s(2,:)=[0.1  0.3  0.6  1    2    4    6    8    10   12]; % k1
        s(3,:)=[0.85 0.88 0.91 0.94 0.97 1    1.05 1.1  1.15 1.2]; % k2 
        s(4,:)=[0.6  0.65 0.7  0.75 0.8  0.85 0.9  0.95 1    1.1]; % eta
        s(5,:)=[0.4  0.5  0.6  0.7  0.8  0.9  1    1.1  1.2  1.3];% delay
        ip=[2 4  % Parameters being compared in 1st (1st row) and 2nd loop (2nd row).
            2 5
            2 1
            1 5
            2 3];
        np=size(s,1);
        nj=size(s,2);
        ep=size(ip,1);
        dl=zeros(np,nj,ep);
        p=p0;
        % k1-eta
        for e=1:ep
            fprintf('Combination of parameters: %i/%i.(%i,%i). Total runs: %i\n',e,ep,ip(e,1),ip(e,2),nj*nj)
            k=0;
            for m=1:nj
                p1=ip(e,1);
                p(p1)=s(p1,m)*p0(p1);
                for n=1:nj
                    k=k+1;
                    fprintf('Combination %i/%i. Loop %i/%i. \n',e,ep,k,nj^2)
                    p2=ip(e,2);
                    p(p2)=s(p2,n)*p0(p2);
                    dl(m,n,e)=RunCase(Geo,p,Set);
                end
                p=p0;
            end
        end
        save('Condensation','dl','s','p0','Set','ip');
    end
    if Run(3)>0
        f0=0;
        if Run(2)>0
            f0=5;
        end
        load('Condensation','dl','s','p0','Set','ip');
        ep=size(ip,1);
        for e=1:ep
            % Plot condensation level
            figure(f0+e);clf
            p1=ip(e,1);
            p2=ip(e,2);
            x=s(p1,:)*p0(p1);
            y=s(p2,:)*p0(p2);
            surfc(y,x,dl(:,:,e));
            ylabel(par{p1});
            xlabel(par{p2});
            set(gca,'FontSize',FS)
            view(2)
            title('Relative Shortening')
            colorbar
            caxis([0.65 1.0]);
        end
    end
end
fprintf('Program finished\n')
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dl,da,da0,freq,l,L]=RunCase(Geo,p0,Set)
% INPUT
% p0: Set of material parameters: gamma, k1, k2, eta, dealy, L0, l0, f
% OUTPUT:
% dl = difference in length
% da = difference in amplitude
%
Mat.gm=p0(1)*ones(Geo.nelem,1); % 0.6; 0.21 stable; 0.22 unstable
Mat.k1=p0(2)*ones(Geo.nelem,1); % 2;
Mat.k2=p0(3)*ones(Geo.nelem,1); % 3;
Mat.eta=p0(4)*ones(Geo.nelem,1); %eta=15 stable , eta=16 unstable (free end),
Mat.T=p0(5); % 4. Delay
Mat.L0=p0(6)*ones(Geo.nelem,1); % Reference length per element
Mat.l0=p0(7); % Initial apparent length per element
Set.f=-p0(8)*[1 1]; % Initial and final force at right end
% Geo
nx0=Geo.nx0;
ny=Geo.ny;
nelem=Geo.nelem;
ns=Geo.ns;
% Mat
eta=Mat.eta;
k2=Mat.k2;
k1=Mat.k1;
gm=Mat.gm;
T=Mat.T;
L0=Mat.L0;
l0=Mat.l0;
% Set
tend=Set.tend;
dt=Set.dt;
f=Set.f;
% No homogenoeus case
if ns>1
    de=nelem/10:nelem/10:nelem;
    %k1(de)=2*k1(de);
    k2(de)=2*k2(de);
    %eta(de)=1.2*eta(de);
    % gm(de)=1.1*gm(de);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal Automatic settings
Ld=ones(round(T/dt),nelem);
ld=ones(round(T/dt),nelem);

t=0:dt:tend;
nt=length(t);

X=zeros(nt,nelem+1);
X(1,:)=0:l0:nelem*l0;

L=zeros(nt,nelem);
L(1,:)=L0;
Lt=zeros(nt,1); % Total length
Lt(1)=sum(L0);
l=zeros(nt,nelem);

df=f(2)-f(1);
if abs(df)<eps % Constant force
    f=f(1)*ones(nt,1);
else
    f=f(1):df/(nt-1):f(2); % Incremental end force
end
l(1,:)=l0;
if exist('dl','var') % Imposed displacements
    if abs(dl)<eps
        X(:,end)=l(1)*nelem*ones(nt,1);
    else
        X(:,end)=(l(1)*nelem):dl/(nt-1):(l(1)*nelem+dl);
    end
end
ld(1)=l(1);
lmax=nan;
lmin=nan;
nmax=nan; % Frequency
nmin=nan; % Frequency
freq=[];
lmin0=nan;
lmax0=nan;
% Loop in time-increments
for n=1:length(t)-1
    K=zeros(nelem+1);
    b=zeros(nelem+1,1);
    % Loop on elements
    for e=1:nelem
        % Build factors for equilibrium equations
        etae=eta(e)/nelem^2; % in order to avoid discretisation dependent results
        L0e=(k1(e)/etae)*L0(e);
        k1e=(k1(e)+k2(e))/etae;
        k2e=k2(e)/etae;
        fe=f(n+1)/etae;
        % Update L
        L(n+1,e)=dt*gm(e)*(ld(end,e)-Ld(end,e))  +L(n,e);
        % Define stresses and rhs
        s11=k1e+1/dt;
        s12=-k1e;
        s21=-k1e;
        s22=k1e+1/dt;
        b1=L0e+k2e*L(n+1,e)-X(n,e)/dt;
        b2=-L0e-k2e*L(n+1,e)-X(n,e+1)/dt;
        %Assemble
        K(e:e+1,e:e+1)=K(e:e+1,e:e+1) + [s11 s12 ; s21 s22];
        b(e:e+1)=b(e:e+1) + [b1 b2]';
    end
    Lt(n+1)=sum(L(n+1,:));
    b(end)=b(end)-fe;
    if ~exist('dl','var') % free right end
        X(n+1,2:end)=-K(2:end,2:end)\b(2:end);
        % l(n+1)=(dt/(dt*k1+1))*(f(n)+k2*L(n+1)+L0+l(n)/dt);
    else % contrained right end X(:,end)
        X(n+1,2:end-1)=-K(2:end-1,2:end-1)\(b(2:end-1)+X(n+1,end)*K(2:end-1,end));
    end
    l(n+1,:)=X(n+1,2:end)-X(n+1,1:end-1);
    % Update Ld, ld
    %   l(n+1)=dt*k2*L(n+1)+l(n);
    for i=size(Ld,1):-1:2
        Ld(i,:)=Ld(i-1,:);
        ld(i,:)=ld(i-1,:);
    end
    Ld(1,:)=L(n+1,:);
    ld(1,:)=l(n+1,:);
    % Compute length differences
    if n>1
        if X(n,end)>X(n-1,end) && X(n,end)>X(n+1,end)
            if isnan(lmax)
                lmax0=X(n,end);
            end
            if ~isnan(nmax)
                freq=[freq n-nmax];
            end
            lmax=X(n,end);
            nmax=n;
        end
        if X(n,end)<X(n-1,end) && X(n,end)<X(n+1,end)
            if isnan(lmin)
                lmin0=X(n,end);
            end
            lmin=X(n,end);
            if ~isnan(nmax)
                freq=[freq n-nmin];
            end
            nmin=n;
        end
    end
end
lm0=(lmin0+lmax0);
da0=(lmax0-lmin0)/lm0;
dl=(lmin+lmax)/lm0;  % Relative difference in length
da=(lmax-lmin)/lm0;  % Relative difference in amplitude
freq=1/mean(freq);
% Relative units
fprintf('Finished %i increments using %i elements.da0=%e, da=%e, dl=%e\n',length(t),nelem,da0, da,dl);
lmax=max(X(:,end));
% Time evolution of apparent length
if Set.Plot==1
    fprintf('Plotting time evolution. Max length=%f.\n',lmax);
    figure
    clf
    plot(t,X(:,end),'LineWidth',2)
    hold on
    plot(t,Lt,'LineWidth',2)
    grid on
    set(gca,'FontSize',FS)
    xlabel('Time (t)')
    ylabel('l(t), L(t)')
    legend('Current Length l(t)','Rest-Length L(t)')
    set(gcf,'color','w');
    
    % Kymograph
    figure
    clf
    lmax=max(X(:,end));
    nx=nx0*ceil(lmax);
    u=NaN(nx,tend/dt); % Kymograph
    k=1;
    for i=1:nt
        if i>k*ny
            % u(1:nx-floor((lmax-X(i,end))*nx0),k)=X(i,end);
            u(:,k)=MapResults(nx,nx0,X(i,:),l0);
            k=k+1;
        end
    end
    imagesc(u');
    c=colorbar;
    c.Label.String = 'Strain [%]';
    c.Label.VerticalAlignment='bottom';
    colormap( [1 1 1; parula(255)] )
    set(gca,'YDir','normal')
    nYTicks=5;
    nXTicks=4;
    LmaxKymo=500;
    LminKymo=100;
    set(gca, 'XTick', 1:(nx-1)/nXTicks:nx+1, 'XTickLabel', LminKymo:(LmaxKymo-LminKymo)/nXTicks:(LmaxKymo+1))
    set(gca, 'YTick', 1:(tend-1)/dt/nYTicks:tend/dt, 'YTickLabel', 0:tend/nYTicks:tend+1)
    ylabel('Time [s]','FontSize',FS)
    xlabel('AP-axis [um]','FontSize',FS)
    set(gca,'FontSize',FS)
    c.Label.Position=[-0.3 0.88];
    c.Label.FontSize=FS;
end
end
% Local functions
function u=MapResults(nx,nx0,X,l0)
% nx  = total dimension of results
% nx0 = number of components per X unit.
%
%
u=zeros(nx,1);
e=1;
for i=1:nx
    if i>X(end)*nx0
        u(i)=NaN;
    else
        if i>X(e+1)*nx0
            e=e+1;
        end
        u(i)=(X(e+1)-X(e))/l0;
    end
end
end

function createfigure(YData1, XData1, X1, YMatrix1, Y1)
%CREATEFIGURE(YData1, XData1, X1, YMatrix1, Y1)
%  YDATA1:  patch ydata
%  XDATA1:  patch xdata
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 17-Jun-2021 10:37:01

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.15251572327044 0.146596858638743 0.69182389937107 0.788887012328999]);
hold(axes1,'on');

% Create patch
patch('Parent',axes1,'DisplayName','Final Osc. Amplitude','YData',YData1,...
    'XData',XData1,...
    'FaceColor',[0.5 0.5 0.5]);

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineStyle','--','Color',[0 0 0]);
set(plot1(1),'DisplayName','Initial Osc. Amplitude');
set(plot1(2),'DisplayName','Relative VNC length','LineWidth',2,...
    'LineStyle','-',...
    'Color',[0 0.447 0.741]);

% Create ylabel
ylabel('Relative shortening ');

% Create xlabel
xlabel('k_1');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.001 0.1]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 0.95]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',18,'YColor',[0 0 0],'YTick',...
    [0 0.19 0.38 0.57 0.76 0.95],'YTickLabel',...
    {'0','0.19','0.38','0.57','0.76','0.95'});
% Create axes
axes2 = axes('Parent',figure1,...
    'ColorOrder',[0.85 0.325 0.098;0.929 0.694 0.125;0.494 0.184 0.556;0.466 0.674 0.188;0.301 0.745 0.933;0.635 0.078 0.184;0 0.447 0.741],...
    'Position',[0.15251572327044 0.146596858638743 0.69182389937107 0.788887012328999]);
hold(axes2,'on');

% Create plot
plot(X1,Y1,'Parent',axes2,'LineWidth',2,'Color',[0 1 0]);

% Create ylabel
ylabel('Frequency X 10');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes2,[0.001 0.1]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes2,[0.15 0.3]);
% Set the remaining axes properties
set(axes2,'Color','none','FontSize',18,'HitTest','off','YAxisLocation',...
    'right','YColor',[0 0 0],'YTick',[0.15 0.18 0.21 0.24 0.27 0.3],...
    'YTickLabel',{'0.15','0.18','0.21','0.24','0.27','0.3'});
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','southwest','AutoUpdate','off');

end
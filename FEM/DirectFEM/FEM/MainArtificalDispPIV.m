function resultats=MainArtificalDispPIV(Case)
% INPUT:
% Case :
%   =1: Initial linear (10 stesp) and the constant displacements (10
%   steps). For creep test.
%	=2: Constant displacements (20 steps). For stress relaxation test.
% OUTPUT:
% resultats.Xlag: x coordinates of applied x displacement
% resultats.Ylag: y coordinates of applied y displacement
% resultats.Ulag: x velocities
% resultats.Vlag: y velocities
% resultats.time: time points
% SINTAX:
% clearvars;resultats=MainArtificalDispPIV(1);save('results_creep.mat')
% clearvars;resultats=MainArtificalDispPIV(2);save('results_relaxation.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
nt=20; % number time points
Utot=10; % final applied displacement
Xmin=90;
Xmax=400;
Ymin=-30;
Ymax=30;
dY=2;
resultats.time=1:nt;
Y=Ymin:dY:Ymax;
nY=length(Y);
Xlag=zeros(2*nY,nt);
Ylag=zeros(2*nY,nt);
Ulag=zeros(2*nY,nt);
Vlag=zeros(2*nY,nt);
X1=Xmin*ones(1,nY);
X2=Xmax*ones(1,nY);
Zeros=zeros(1,nY);
Ones=ones(size(Zeros));
for t=1:nt
    if Case==1 % linear then constant, for creep testing.
        Ylag(:,t) = [Y Y];
        Vlag(:,t) = [Zeros Zeros];
        if t<=nt/2
            Xlag(:,t) = [X1 X2+Utot/(nt/2)*(t-1)];
            Ulag(:,t) = [Zeros Utot/(nt/2)*Ones];
        else
            Xlag(:,t) = [X1 X2+Utot];
            Ulag(:,t) = [Zeros Zeros];
        end
    elseif Case==2 % all constant, for stress relaxation.
        Ylag(:,t) = [Y Y];
        Vlag(:,t) = [Zeros Zeros];
        if t==1
            Xlag(:,t) = [X1 X2];
            Ulag(:,t) = [Zeros Utot*Ones];
        else
            Xlag(:,t) = [X1 X2+Utot];
            Ulag(:,t) = [Zeros Zeros];
        end
    end
end
resultats.Xlag=Xlag;
resultats.Ylag=Ylag;
resultats.Ulag=Ulag;
resultats.Vlag=Vlag;

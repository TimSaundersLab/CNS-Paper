function [aux,En,Sn,Snv,u,uExp]=FEM(BC,C,En,Set,Sn,X,X0)
%d,dim,E,h,Set,Tx,Ty,Tz,v)
% Solves Forward problem using FEM.
% Assumed (fully or partially) constrained problem.
% INPUT:
% E v
% En(e,gp,:) = Previous strains of element e at gauss  point gp. Exx,Eyy, Eyy, Exy,Exy,Exz,Eyz
% Sn(e,gp,:) = Previous stresses of element e at gauss point. Sxx,Syy, Syy,Sxy,Sxy,Sxz,Syz
% OUTPUT:
% aux = dofs of constrained dofs (applied displacements of PIV)
% uExp = experimental displacements
% En(e,gp,:) = Strains of element e at gauss  point gp. Exx,Eyy, Eyy, Exy,Exy,Exz,Eyz
% Sn(e,gp,:) = Stresses of element e at gauss point. Sxx,Syy, Szz,Sxy,Sxy,Sxz,Syz
% Snv(e,gp,:) = Viscous Stresses of element e at gauss point. Sxx,Syy, Szz,Sxy,Sxy,Sxz,Syz
if Set.AccumulateU % In viscoelasticity or accumulated elasticity
    x=X0;
else
    x=X;
end
G=[Set.E Set.v];
nodes=size(x,1);
dim=size(X,2);
nelem=size(C,1);
nnod=size(C,2);
%    A=sparse(nodes*dim,nelemb*dimT);
%    Ae=Aeq4e(X(C(1,1:nnod),:),Set.Tz); % Matrix such that fe=Ae*[Tx ; Ty; Tz]
dof=dim*nodes;
f=zeros(dof,1);
ig=zeros(nelem*nnod*dim*nnod*dim,1);
jg=ig;
vg=ig;
k=0;
for e=1:nelem
    lnod=C(e,1:nnod);
    if nnod==8
        [Ke,qe]=keh8e(En(e,:,:),G,Sn(e,:,:),x(lnod,:),Set);
    elseif nnod==4
        [Ke,qe]=keq4e(En(e,:,:),G,Sn(e,:,:),x(lnod,:),Set);
    else
        error('elements with %i nodes not implemented',nnod)
    end
    dofe=kron(lnod-1,ones(1,dim))*dim+kron(ones(1,nnod),1:dim);
    for i=1:nnod*dim
        for j=1:nnod*dim
            if abs(Ke(i,j))>eps
                k=k+1;
                ig(k)=dofe(i);
                jg(k)=dofe(j);
                vg(k)=Ke(i,j);
            end
        end
    end
    f(dofe)=f(dofe)-qe;
end
Kc=sparse(ig(1:k),jg(1:k),vg(1:k),dof,dof);
% Apply constraints
[aux,naux,uaux]=ApplyBC(BC,dim,dof);
% FORWARD PROBLEM. Compute nodal displacements
u=zeros(size(Kc,1),1);
u(aux)=uaux;
uExp=ones(size(u))*(min(uaux)+max(uaux))/2; % All set to average imposed values
uExp(aux)=uaux;
uExp=reshape(uExp,size(x'))';
f=f-Kc*u;
f(aux)=[];
Kc(aux,:)=[];
Kc(:,aux)=[];
us=Solve(Kc,f);
u(naux)=us;
u=reshape(u,size(x'))';
Ep=En; % Previous strains
Sp=Sn; % Previous stresses
Snv=0*Sn;
Een=zeros(size(En,2),size(En,3));
% Compute Stresses and Strains
for e=1:nelem
    lnod=C(e,1:nnod);
    Ue=u(lnod,:); % Total applied displacements at current time-step
    Xe=x(lnod,:);
    Een(:,:)=Ep(e,:,:); % elemental Previous strain
    Sen(:,:)=Sp(e,:,:); % elemental Previous stresses
    if nnod==8
        [~,Ee,Se,Sev] = qeh8e(Een,G,Ue,Sen,Set,Xe);
    elseif nnod==4
        [~,Ee,Se,Sev] = qeq4e(Een,G,Ue,Sen,Set,Xe);
    end
    En(e,:,:)=Ee;
    Sn(e,:,:)=Se;
    Snv(e,:,:)=Sev;
end
end
%%
function [aux,naux,uaux]=ApplyBC(BC,dim,dof)
% OUTPUT:
% aux = list of constratined dof
% naux = listof free dof
% uaxu = valus of displacement in constrainted dof
ndof=size(BC.u,1);
uaux=1:dof;
aux=-uaux; % Constrained dof
for i=1:ndof
    node=BC.u(i,1);
    dof=BC.u(i,2);
    u=BC.u(i,3);
    idof=(node-1)*dim+dof;
    aux(idof)=idof;
    uaux(idof)=u;
end
naux=-aux(aux<0);
uaux(naux)=[];
aux(naux)=[];
end


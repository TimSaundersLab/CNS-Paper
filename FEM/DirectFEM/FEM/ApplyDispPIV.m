function [BC,nApplied]=ApplyDispPIV(dx,PIV,Set,X,X0)
% Applies displacement of PIV data. PIV struct contains:
% PIV.X(i,:): x annd y positions of PIV node i
% PIV.U(i,:): x annd y displacement of PIV node i
% X(i,:) : x coordinates of mesh node i. 3rd dimension is computed.
nPIV=size(PIV.X,1);
BCu=zeros(0,3);
nnodes=1:size(X,1);
nApplied=0;
for i=1:nPIV
    XP=PIV.X(i,:);
    U=PIV.U(i,:);
    nodes=nnodes(abs(XP(1)-X(:,1))<Set.PIVtol*dx(1) & abs(XP(2)-X(:,2))<Set.PIVtol*dx(2));
    BCu=[BCu
        nodes'   ones(length(nodes),1) ones(length(nodes),1)*U(1)];
    if Set.ConstrainY
        BCu=[BCu
            nodes' 2*ones(length(nodes),1) ones(length(nodes),1)*U(2)];
    end
    nApplied=nApplied + length(nodes);
end
% Make sure that at least one dof is constrained along Y
xmin=min(X0(:,1));
xmax=max(X0(:,1));
nodes=1:size(X0,1);
nodeM=nodes(X0(:,1)==xmax);
nodem=nodes(X0(:,1)==xmin);
if ~Set.ConstrainY
    BCu=[BCu
        nodeM(1) 2 0 
        nodem(1) 2 0];
end
% Make sure that at x=xmin, boundary remains flat (Sensitive to some
% missing displacements). 
% Apply mean displacements of measurements for points that x0 < xmin +dx(1)
xmin=min(X0(:,1));
nodes=nnodes(X0(:,1)-xmin<dx(1));
u=0;
k=0;
iBC=zeros(size(BCu,1),1);
for i=1:size(BCu,1)
    if BCu(i,2)==1
        if min(abs(BCu(i,1)-nodes))==0
            u=u+BCu(i,3);
            k=k+1;
            iBC(k)=i;
        end
    end
end
iBC(k+1:end)=[];
u=u/k;
BCu(iBC,:)=[BCu(iBC,1) BCu(iBC,2) u*ones(size(iBC))];
% Add previous displacements if accumulated (u accounts for all accumulated displacements)
if Set.AccumulateU
    uPrev=zeros(size(BCu,1),1);
    for i=1:size(BCu,1)
        uPrev(i)=X(BCu(i,1),BCu(i,2))-X0(BCu(i,1),BCu(i,2));
    end
    BC.u=[BCu(:,1) BCu(:,2) BCu(:,3)+uPrev];
else
    BC.u=BCu;
end
if nApplied ==0
    warning('Could not apply displacements from PIV');
    BC.u=[1 3 0];
else
    ymin=min(X0(:,2));
    ymax=max(X0(:,2));
    nodes=1:size(X0,1);
    nodeM=nodes(X0(:,2)==ymax);
    nodem=nodes(X0(:,2)==ymin);
    BC.u=[BC.u
        nodem(1) 3 0 
        nodeM(1) 3 0];
end
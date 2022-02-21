function [C,X]=TrimMesh(C,Trim,X)
% Trims mesh according to bounding box in Trim
% INPUT:
% trim=[xmin xmax ymin ymax zmin zmax] coordinates of trimmming box for mesh. NaN=no trimming.
% OUTPUT:
% C(e,:)=connectivity of element e
% X(i,:)=(x,y,x) coordinates of node i
% Select nodes on x
nodesc=cell(6,1);
for id=0:2 % loop on dimension
    for s=1:2 % Loop on max min
        i=s+id*2;
        if ~isnan(Trim(i))
            if s==1
            nodesc{i}=find(X(:,id+1)<Trim(i));
            else
            nodesc{i}=find(X(:,id+1)>Trim(i));
            end
        end
    end
end
nodes=unique([nodesc{1}; nodesc{2}; nodesc{3}; nodesc{4}; nodesc{5}; nodesc{6}]);
% Remove nodes 
nnodes=size(X,1);
allnodes=(1:nnodes)';
OldNode=zeros(nnodes,1);
NewNode=zeros(nnodes,1);
allnodes(nodes)=0;
k=0;
r=0;
for i=1:nnodes
    if allnodes(i)==0
        r=r+1;
    else
        k=k+1;
        NewNode(k)=i;
    end
    OldNode(i)=k;
end    
NewNode(k+1:end)=[];
% Remove elements and renumber connectivity
nelem=size(C,1);
nnod=size(C,2);
OldElem=zeros(nelem,1);
NewElem=zeros(nelem,1);
k=0;
r=0;
for e=1:nelem
    Remove=0;
    for i=1:nnod
        if min(abs(nodes-C(e,i)))==0
            Remove=1;
            break;
        end
    end
    if Remove
        r=r+1;
    else
        k=k+1;
        NewElem(k)=e;
    end
    OldElem(e)=k;
    C(e,:)=OldNode(C(e,:)); % Renumber
end
NewElem(k+1:end)=[];
X=X(NewNode,:);
C=C(NewElem,:);

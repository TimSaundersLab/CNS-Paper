function [C,dX,X,Error]=ReadMesh(Set)
% Reads Mesh data from GID datafile using 'embryo' problemtype
% OUTPUT:
% C(e,:)=connectivity of element e
% dX=mean size of one element
% X(i,:)=(x,y,x) coordinates of node i
Error=0;
if Set.Dimensions==2
    dX=Set.dx*[1 1];
    nx=floor(400/dX(1));
    ny=floor(72/dX(2));
    [C,X]=GenerateMesh(dX(1),dX(2),[],nx,ny);
    C=C(:,1:end-1); % Remove elements number
    X(:,2)=X(:,2)-(max(X(:,2))-min(X(:,2)))/2; % Center on Y
else
    fid=fopen(Set.MeshFile);
    if fid<0
        fid=fopen(strcat(Set.MeshFile,'.dat'));
        if fid<0
            error('Unable to open file %s or %s',Set.MeshFile,strcat(Set.MeshFile,'.dat'));
        end
    end
    txt=fgets(fid);
    nnodes=0;
    nelem=0;
    reading=true;
    while reading && ~feof(fid)
        if length(txt)>15
            if strcmp(txt(1:16),'NUMBER OF POINTS')==1
                break;
            end
        end
        txt=fgets(fid);
    end
    if feof(fid)
        error(sprintf('%s not found in file %s\n','NUMBER OF POINTS',Set.MeshFile));
    end
    txt=fgets(fid);
    l=str2num(txt);
    if ~isnan(l)
        nnodes=l(1);
    end
    reading=true;
    while reading && ~feof(fid)
        if length(txt)>22
            if strcmp(txt(1:23),'NUMBER OF TYPE ELEMENTS')==1
                break;
            end
        end
        txt=fgets(fid);
    end
    if feof(fid)
        error(sprintf('%s not found in file %s\n','NUMBER OF TYPE ELEMENTS',Set.MeshFile));
    end
    txt=fgets(fid);
    l=str2num(txt);
    if ~isnan(l)
        if length(l)>3
            nelem=l(4);
        end
    end
    if nnodes==0 || nelem==0
        Error=1;
        return;
    end
    X=zeros(nnodes,3);
    C=zeros(nelem,8);
    reading=true;
    while reading
        if length(txt)>9
            if strcmp(txt(1:10),'NODE COORD')==1
                break;
            end
        end
        txt=fgets(fid);
    end
    % Read Nodes
    for i=1:nnodes
        txt=fgets(fid);
        l=str2num(txt);
        X(i,:)=l(2:4);
    end
    reading=true;
    while reading
        if length(txt)>8
            if strcmp(txt(1:9),'Hexahedra')==1
                break;
            end
        end
        txt=fgets(fid);
    end
    dX=zeros(3,1);
    for i=1:nelem
        txt=fgets(fid);
        l=str2num(txt);
        Ce=l(2:9);
        C(i,:)=Ce;
        dX=dX+[max(X(Ce,1))-min(X(Ce,1))
            max(X(Ce,2))-min(X(Ce,2))
            max(X(Ce,3))-min(X(Ce,3))];
    end
    dX=dX/nelem;
    fclose(fid);
    if Set.Long==3
        X=[X(:,3) X(:,1) X(:,2)];
    end
end
end
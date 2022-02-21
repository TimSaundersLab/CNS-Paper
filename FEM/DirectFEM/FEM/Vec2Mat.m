function t=Vec2Mat(dim,t,auxT)
%% Transforms vector of tractions in vector format (values of tractions 
%  for each dof) into Mat format (values of tractions in each direction per row):
%   t=[t1 t1 t1
%      t2 t2 t3
%       ...
%      tN tN tN]
% INPUT:
%  dim   = space dimensions
%  t     = array of traction values for each traction dof.
%  auxT  = degrees of freedom in t (in case not all dof are included in t)
% OUTPUT:
% t=Traction values in Mat format
dofT=size(t,1);
t_copy=t;
if ~exist('auxT','var') || min(size(auxT))==0
    auxT=1:dofT;
end
nodes=floor((max(auxT)-1)/dim)+1;
t=zeros(nodes,dim);
for i=1:dofT
    dof=auxT(i);
    node=floor((dof-1)/dim)+1;
    dof=dof-(node-1)*dim;
    t(node,dof)=t_copy(i);
end
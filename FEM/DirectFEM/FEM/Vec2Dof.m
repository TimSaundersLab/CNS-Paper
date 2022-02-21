function t=Vec2Dof(dim,t,auxT)
%% Transforms vector of tractions in array format (values of tractions 
%  for each dof) into dof format:
%      t=[node direction[1,2,3] value]
% INPUT:
% auxT  = array with displacement dof for each traction value
%         auxT(i)=displacement dof of traction in t(i)
% t     = array of traction values for each traction dof.
% OUTPUT:
% t=Traction values in dof format. Each row:
%   t=[node direction[1,2,3] value]
%%
dofT=size(t,1);
t_copy=t;
t=zeros(dofT,3);
if ~exist('auxT','var') || min(size(auxT))==0
    auxT=1:dofT;
end
for i=1:dofT
    dof=auxT(i);
    node=floor((dof-1)/dim)+1;
    dof=dof-(node-1)*dim;
    t(i,:)=[node dof t_copy(i)];
end
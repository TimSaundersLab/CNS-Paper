function [x,error]=SolveMA57m(A,b)
%% Calls solver function that uses MA57 subroutines
% from Central Laboratory
%               of the Research Councils
% ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
% Licence, see http://hsl.rl.ac.uk/acuk/cou.html
% Uses matrices in sparse component format: row(i),col(i),val(i)
%
% ONLY WORKS FOR SYMMETRIC MATRICES
% USAGE:
%
%n=50;A=sparse(rand(n));A=A+A';b=ones(n,1); [x,err]=SolveMA57m(A,b)
%
% INPUT:
% A
% acolt : array with column indeces
% arowt : array with row indeces
% avalt : array with values
%
% OUTPUT:
% x    : solution vector
% error ;
%       =1: PROBLEM IS INFEASIBLE
%       =2 CONVERGENCE FAILED TO PROCEED
%       =3 LOSS OF ACCURACY IN SOLVER
%

if ~issparse(A)
    warning('SolveMA57 only works for symmetric sparse matrices')
    x=zeros(size(A,1),1);
    return
end
[acol,arow,aval]=find(A);
[m,n]=size(A);
clearvars A;
li=zeros(length(acol),1);
k=0;
% Make sure that only lower part is sent
for i=1:length(acol)
    if acol(i)<=arow(i)
        k=k+1;
        li(k)=i;
    end
end
[x,error]=SolveMA57(acol(li(1:k))',arow(li(1:k))',aval(li(1:k))',b,m,n);
end
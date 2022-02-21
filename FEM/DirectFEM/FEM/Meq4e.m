function Me = Meq4e(Xe,Tz)
%***************************************************
% AeQ4E:
%   Creates the element mass matrix of elastic
%   quadrilateral 4-node element such that
% Syntax:
%   Me = Meq4e(Xe,Ge)
% Input:
%   Xe   : coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4]
%   Tz   : =0, z tractions assumed zero, Tz=0; if >0, Tz included.
% Output:
%   Ae   : element traction matrix.
% Date:
%   Version 1.0    04.05.13
%***************************************************

% Gauss abscissae and weights.
r = [-1 1]/sqrt(3);
w = [ 1 1];


% determine number of nodes per element
dimT=2;
if Tz
    dimT=3;
end
% Initialize stiffness matrix.
M=zeros(4,4);
I=eye(dimT);
% Gauss integration of stiffness matrix.
for i = 1:2
    for j = 1:2
        
        % Function values
        N= [ (1-r(i))*(1-r(j))
            (1+r(i))*(1-r(j))
            (1+r(i))*(1+r(j))
            (1-r(i))*(1+r(j)) ]/4;
        % Parametric derivatives:
        dN = [ -(1-r(j))  (1-r(j))  (1+r(j)) -(1+r(j))
            -(1-r(i)) -(1+r(i))  (1+r(i))  (1-r(i)) ]/4;
        
        % transform to global coordinates
        Jt = dN*Xe(:,1:2);
        
        M = M + w(i)*w(j)*(N*N')*det(Jt);
        
    end
end
Me=zeros(dimT*4,dimT*4);
for i=1:4
    for j=1:4
        Me((i-1)*dimT+1:i*dimT,(j-1)*dimT+1:j*dimT)=M(i,j)*I;
    end
end


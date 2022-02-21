function [Tx,Ty,Tz]=Elem2Node(dx,dy,Tx,Ty,Tz)
%***************************************************
% Transforms elemental tractions to nodal values
% Syntax:
%   Ae = Aeq4e(Xe,Ge)
% Input:
%   Tx,Ty   : elemental values of X and Y tractions
% Output:
%   Tx,Ty   : nodal forces in X and Y equivalent to Tx and Ty
% Date:
%   Version 1.0    04.05.14
%***************************************************
A=dx*dy;
Tx=Tx*A;
Ty=Ty*A;
Tz=Tz*A;
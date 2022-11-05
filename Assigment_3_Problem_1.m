%% Assignment 3 Problem 1
clear variables; clc; close all; format short
syms xi eta
N = [(1/4)*(1-xi)*(1-eta) (1/4)*(1+xi)*(1-eta) (1/4)*(1+xi)*(1+eta) (1/4)*(1-xi)*(1+eta)];

for i = 1:length(N)
    N_d_xi(i) = diff(N(1,i),xi);
    N_d_eta(i) = diff(N(1,i),eta);

    N_d(1,i) = N_d_xi(i);
    N_d(2,i) = N_d_eta(i);
end

%Global nodal coordinates
xy=[-15 -10;
    15 -10;
    15 10;
    -15 10]/2;

% Jacobian
J = N_d*xy;
Det_J = double(det(J));

% Strain Displacement Section %%%%%%%%%%%%%%%%%%%%%%%%%%
syms x y
% Gamma matrix (Inverse of Jacobian)
Gamma = inv(J);

% Conversion to x-y shape functions
for i = 1:length(N)
    N_xy(:,i) = Gamma*N_d(:,i);
end

% B strain-displacement matrix
B = cell2sym(cell(3,2*length(N)));
% Inputting the shape functions in the correct order
for i = 1:length(N)
    B(1,i*2-1) = N_xy(1,i);
    B(2,i*2-1) = 0;
    B(3,i*2-1) = N_xy(2,i);
    B(1,i*2) = 0;
    B(2,i*2) = N_xy(2,i);
    B(3,i*2) = N_xy(1,i);
end

% Stiffness matrix k integration section %%%%%%%%%%%%%%%%%%%
% Values
t = 1; %Thickness [mm]
E = 0.7*10^3; %Young's moduli [MPa]
v = 0.45; %Poissons ratio

% [B] has been found

% For plane stress (Sigma_z = 0) the Constitutive matrix is given by
E_m = E/(1-v^2).*[1 v 0; v 1 0; 0 0 (1-v)/2];

nGP = 2; %Number of Gauss Points
k = zeros(length(N)*2); %Size of k = Nodes*D.O.F.
[GP,W]  = GaussTable(nGP); % Gauss sampling points and weights
%Double nested for loop
for i = 1:nGP %Integration in the xi-direction
    for j = 1:nGP %Integration in the eta-direction
        B_temp = subs(B,xi,GP(i)); % Evaluating the strain-displacement matrix
        B_temp = double(subs(B_temp,eta,GP(j))); %And then for eta
        k = k + transpose(B_temp)*E_m*B_temp*t*Det_J*W(i)*W(j); %Calculating the product and summing
        % k [N/mm]
    end
end

%% GaussTable
function [xi,w] = GaussTable(n)
% Returns sampling points and weights for for Gauss quadrature of order n
%
% INPUTS:
%   nGP:      order of quadrature
%
% OUTPUTS:
%    xi:      sampling points
%     W:      weight functions


switch n
    case 1
        xi = 0;
        w  = 2;
    case 2
        xi = [-sqrt(1/3) sqrt(1/3)];
        w  = [1 1];
    case 3
        xi = [-sqrt(3/5) 0 sqrt(3/5)];
        w  = [5/9 8/9 5/9];
    case 4
        xi(1) = -sqrt(3/7 +2/7*sqrt(6/5));
        xi(2) = -sqrt(3/7 -2/7*sqrt(6/5));
        xi(3) =  sqrt(3/7 -2/7*sqrt(6/5));
        xi(4) =  sqrt(3/7 +2/7*sqrt(6/5));
        w(1) = (18-sqrt(30))/36;
        w(2) = (18+sqrt(30))/36;
        w(3) = (18+sqrt(30))/36;
        w(4) = (18-sqrt(30))/36;
    case 5
        xi(1) = -1/3*sqrt(5+2*sqrt(10/7));
        xi(2) = -1/3*sqrt(5-2*sqrt(10/7));
        xi(3) = 0;
        xi(4) =  1/3*sqrt(5-2*sqrt(10/7));
        xi(5) =  1/3*sqrt(5+2*sqrt(10/7));
        w(1) = (322-13*sqrt(70))/900;
        w(2) = (322+13*sqrt(70))/900;
        w(3) = 128/225;
        w(4) = (322+13*sqrt(70))/900;
        w(5) = (322-13*sqrt(70))/900;
end

end



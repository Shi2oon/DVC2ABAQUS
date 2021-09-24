function [disp_notchRef] = dispNotchRef(displacements_col, theta_deg)
% (c) MSL Jordan, University of Oxford 2015

%{
% TEST VAR
displacements_col = big_new_data_natural(:,4:6);
theta_deg = 9; % Angle <x-axis and x'-axis (|| notch)>

%% Process
% 0. Convert theta to radians
% 1. Calculate column vector 'a' (angular difference <vector and
% x'-axis> = <vector and x-axis> - theta
% 2. Calculate column vector 'M' (magnitude of vector)
% 3. Calculate column-wise displacements in notch reference frame (z-column
% constant)
%}
%% 0.
if size(theta_deg,2) ==3
    t = degtorad(-theta_deg(3)); 
    disp '##ISSUE 3-rot/1-rot. Identified 6 Sept 2014';
else
    t = degtorad(-theta_deg); 
end
%% 1.
a = atan2(displacements_col(:,2),displacements_col(:,1))-t;
ca = cos(a);
sa = sin(a);
%% 2.
M = sqrt(displacements_col(:,2).^2+displacements_col(:,1).^2);

%% 3.
disp_notchRef = displacements_col;
disp_notchRef(:,1) = M.*ca;
disp_notchRef(:,2) = M.*sa;

end
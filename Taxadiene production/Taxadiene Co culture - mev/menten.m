function [vg,vz,vi] = menten(G,Z,I)
kg = 0.0027; % % g/L 
vg_max = 10.5; % mmol/g/h 10.57.3
% Kig = 0.5; %g/l
% Kis = 10; %g/l
vg = vg_max.*(G./(kg+G));%%*(1/(1+(S/Kis)));
%
kz = 0.01; %mmol/l
vz_max = 9; % mmol/g/h  8
vz = vz_max.*Z./(kz+Z);
%
ki = 0.0001; % % g/L 
vi_max = 10.5; % mmol/g/h
vi = vi_max.*(I./(ki+I));%%*(1/(1+(S/Kis)));
end

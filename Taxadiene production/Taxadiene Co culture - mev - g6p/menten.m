function[vg,vi_m,vi_g,vz] = menten(G,Im,Ig,Z)
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
ki = 0.000000001; % % g/L
vi_max = 10.5; % mmol/g/h
vi_m = vi_max.*(Im./(ki+Im));%%*(1/(1+(S/Kis)));


ki = 0.00000027; % % g/L
vi_max = 10.5; % mmol/g/h
vi_g = vi_max.*(Ig./(ki+Ig));%%*(1/(1+(S/Kis)));

end
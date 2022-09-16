function [vg,vz,vi,vphe,vtyr] = menten(G,Z,I,phe,tyr)
kg = 0.0027; % % g/L 
vg_max = 10.5; % mmol/g/h 10.57.3
vg = vg_max.*(G./(kg+G));%%*(1/(1+(S/Kis)));
%
kz = 0.01; %mmol/l
vz_max = 9; % mmol/g/h  8
vz = vz_max.*Z./(kz+Z);
%
ki = 0.001; % % g/L 
vi_max = 10.5; % mmol/g/h
vi = vi_max.*(I./(ki+I));%%*(1/(1+(S/Kis)));
%
kph = 0.001; % % g/L 
vph_max = 10.5; % mmol/g/h
vphe = vph_max.*(phe./(kph+phe));%%*(1/(1+(S/Kis)));
%
kty = 0.001; % % g/L 
vtyr_max = 10.5; % mmol/g/h
vtyr = vtyr_max.*(tyr./(kty+tyr));%%*(1/(1+(S/Kis)));
end

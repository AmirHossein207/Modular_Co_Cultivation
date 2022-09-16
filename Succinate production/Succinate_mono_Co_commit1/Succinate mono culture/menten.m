function [vg,vo] = menten(G)
kg = 0.0027; % % g/L 
vg_max = 10.5; % mmol/g/h 10.57.3
% Kig = 0.5; %g/l
% Kis = 10; %g/l
vg = vg_max.*(G./(kg+G));%%*(1/(1+(S/Kis)));
%
ko = 0.024; %mmol/l
vo_max = 15; % mmol/g/h  8
vo = vo_max;%.*O./(ko+O);

end

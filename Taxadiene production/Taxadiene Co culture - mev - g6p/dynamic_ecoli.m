function [eqns] = dynamic_ecoli(t,y)
global model_up model_down Mw

if y(3) < 0
    y(3) = 0;
end
if y(5) < 0
    y(5) = 0;
end
if y(4) < 0
    y(4) = 0;
end
if y(6) < 0
    y(6) = 0;
end
if y(7) < 0
    y(7) = 0;
end

X_up = y(1); G = y(3); Z = y(4);
X_down = y(2); T = y(5);Im = y(6);Ig = y(7);

[vg,vi_m,vi_g,vz] = menten(G,Im,Ig,Z);
%% int production

model_up = changeRxnBounds(model_up,{'EX_glc__D_e','EX_o2_e',...
    'EX_xyl__D_e','EX_mva_e','EX_g6p_e'},...
    [-vg -0.8*15 0 0 0],'l');
model_up = changeObjective(model_up,'BIOMASS_Ec_iJO1366_core_53p95M');
[F_up,v_g_up,v_int_up_m,v_int_up_g,v_o,vz_up] = fba_mev(model_up);
% v_int_up

%% txdn production
model_down = changeObjective(model_down,'BIOMASS_Ec_iJO1366_core_53p95M');

model_down = changeRxnBounds(model_down,{'EX_glc__D_e','EX_o2_e',...
    'EX_mva_e','EX_g6p_e','EX_xyl__D_e'},...
    [0 -0.8*15 -vi_m -vi_g -vz],{'l','l','b','l','l'});

[F_down,v_g_down,v_int_down_m,v_int_down_g,v_txdn,v_o,vz_down] = fba_txdn(model_down);

%% equations

eqns = zeros(7,1);
eqns(1,1) = F_up*X_up;
eqns(2,1) = F_down*X_down;
eqns(3,1) = v_g_up*X_up*Mw(1) + v_g_down*X_down*Mw(1);
eqns(4,1) = vz_up*X_up*Mw(6) + vz_down*X_down*Mw(6);
eqns(5,1) = v_txdn*X_down*Mw(3);
eqns(6,1) = v_int_up_m*X_up*Mw(4) + v_int_down_m*X_down*Mw(4);
eqns(7,1) = v_int_up_g*X_up*Mw(5) + v_int_down_g*X_down*Mw(5);



end

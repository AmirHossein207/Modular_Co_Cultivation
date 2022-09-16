function [eqns] = dynamic_ecoli(t,y)
global model_up model_down Mw

if y(3) < 0
    y(3) = 0;
end
if y(4) < 0
    y(4) = 0;
end
if y(5) < 0
    y(5) = 0;
end
if y(6) < 0
    y(6) = 0;
end
if y(7) < 0
    y(7) = 0;
end
if y(8) < 0
    y(8) = 0;
end
X_up = y(1); G = y(3); Z = y(4); 
X_down = y(2);S = y(5); I = y(6);
phe = y(7);tyr = y(8);

[vg,vz,vi,vphe,vtyr] = menten(G,Z,I,phe,tyr);
%% int production

model_up = changeRxnBounds(model_up,{'EX_glc__D_e','EX_o2_e','EX_xyl__D_e',...
    'EX_phe__L_e'}...
    ,[0 -0.8*15 -vz -vphe],{'l','l','l','l'});
model_up = changeObjective(model_up,'BIOMASS_Ec_iJO1366_core_53p95M');
[F_up,v_g_up,vz_up,v_int_up,v_o_up,v_phe_up,v_tyr_up] = fba_tyr(model_up);
% v_g_up

%% txdn production

model_down = changeRxnBounds(model_down,{'EX_glc__D_e','EX_o2_e','EX_4_tyrosol_e'...
    ,'EX_xyl__D_e','EX_tyr__L_e'...
    },[-vg -0.8*15 -vi 0 -vtyr],{'l','l','l','l','l'});

[F_down,v_g_down,vz_down,v_int_down,v_salid,v_o_down,v_phe_down,...
    v_tyr_down] = fba_salid(model_down);

%% equations

eqns = zeros(8,1);
eqns(1,1) = F_up*X_up;
eqns(2,1) = F_down*X_down;
eqns(3,1) = v_g_up*X_up*Mw(1) + v_g_down*X_down*Mw(1);
eqns(4,1) = vz_up*X_up*Mw(5) + vz_down*X_down*Mw(5);
eqns(5,1) = v_salid*X_down*Mw(2);
eqns(6,1) = v_int_up*X_up*Mw(3) + v_int_down*X_down*Mw(3);
eqns(7,1) = v_phe_up*X_up*Mw(4) + v_phe_down*X_down*Mw(4);
eqns(8,1) = v_tyr_up*X_up*Mw(6) + v_tyr_down*X_down*Mw(6);

end

function [eqns] = dynamic_ecoli(t,y)
global model_up model_down Mw
if y(1) < 0
    y(1) = 0;
end
if y(2) < 0
    y(2) = 0;
end
if y(3) < 0
    y(3) = 0;
end
if y(5) < 0
    y(5) = 0;
end
if y(4) < 0
    y(4) = 0;
end
if y(5) < 0
    y(5) = 0;
end

X_up = y(1); G = y(3); Z = y(4); 
X_down = y(2); I = y(6);S = y(5);

[vg,vz,vi] = menten(G,Z,I);
%% pep production

model_up = changeRxnBounds(model_up,{'EX_glc__D_e','EX_o2_e','EX_xyl__D_e','EX_g6p_e'}...
    ,[-vg -0.8*15 0 0],{'l','l','l','l'});
model_up = changeObjective(model_up,'BIOMASS_Ec_iJO1366_core_53p95M');
[F_up,v_g_up,v_int_up,vz_up,~,v_o_up] = fba(model_up);

%% succ production

model_down = changeRxnBounds(model_down,{'EX_glc__D_e','EX_o2_e','EX_g6p_e'...
    ,'EX_xyl__D_e'},[0 -0.8*15 -vi -vz],{'l','l','l','l'});
[F_down,v_g_down,v_int_down,vz_down,v_succ,v_o_down] = fba(model_down);

% disp(v_int_down)
%% equations

eqns = zeros(6,1);
eqns(1,1) = F_up*X_up;
eqns(2,1) = F_down*X_down;
eqns(3,1) = v_g_up*X_up*Mw(1) + v_g_down*X_down*Mw(1);
eqns(4,1) = vz_up*X_up*Mw(4) + vz_down*X_down*Mw(4);
eqns(5,1) = v_succ*X_down*Mw(2);
eqns(6,1) = v_int_up*X_up*Mw(3) + v_int_down*X_down*Mw(3);


end

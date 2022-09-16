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
if y(5) < 0
    y(5) = 0;
end
if y(6) < 0
    y(6) = 0;
end

X_up = y(1); G = y(3); X = y(4); 
X_down = y(2); I = y(6); T = y(5);
%% int production
 [vg,vz,vi] = menten(G,X,I);
model_up = changeRxnBounds(model_up,{'EX_glc__D_e','EX_o2_e','EX_xyl__D_e'},...
    [-vg -0.8*15 0],'l');
[F_up,v_g_up,v_int_up,v_o_up,vz_up] = fba_mev(model_up);
% v_int_up

%% txdn production


model_down = changeRxnBounds(model_down,{'EX_glc__D_e','EX_o2_e',...
    'EX_xyl__D_e','EX_mevR_e'},...
    [0 -0.8*15 -vz -vi],{'l','l','l','b'});

[F_down,v_g_down,v_int_down,v_txdn,v_o_down,vz_down] = fba_txdn(model_down);

% v_int_down
%% equations

eqns = zeros(6,1);
eqns(1,1) = F_up*X_up;
eqns(2,1) = F_down*X_down;
eqns(3,1) = v_g_up*X_up*Mw(1) + v_g_down*X_down*Mw(1);
eqns(4,1) = vz_up*X_up*Mw(5) + vz_down*X_down*Mw(5);
eqns(5,1) = v_txdn*X_down*Mw(3);
eqns(6,1) = v_int_up*X_up*Mw(4) + v_int_down*X_down*Mw(4);



end

function [eqns] = dynamic_ecoli(t,y)
global model Mw
% DO = 0.30;
if y(2) < 0
    y(2) = 0;
end
if y(3) < 0
    y(3) = 0;
end

X = y(1); G = y(2); S = y(3);

[vg,vo] = menten(G,S);
model = changeRxnBounds(model, {'EX_glc__D_e','EX_xyl__D_e','EX_o2_e'...
    },[-vg,0,-15*0.8]...
    ,{'l','l','l'});
[F,v_g,v_txdn,v_o] = fba(model);
eqns = zeros(3,1);
eqns(1,1) = F*X;
eqns(2,1) = v_g*X*Mw(1);
eqns(3,1) = v_txdn*X*Mw(3);


end

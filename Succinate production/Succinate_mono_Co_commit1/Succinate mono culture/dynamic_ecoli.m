function [eqns] = dynamic_ecoli(t,y)
global model_1 Mw

if y(2) < 0
    y(2) = 0;
end
if y(3) < 0
    y(3) = 0;
end

X = y(1); G = y(2); S = y(3);

[vg,vo] = menten(G);
model_1 = changeRxnBounds(model_1,{'EX_glc__D_e','EX_o2_e'}...
    ,[-vg -15*0.8],{'l','l'});
[F,v_g,v_succ,v_o] = fba(model_1);

eqns = zeros(3,1);
eqns(1,1) = F*X;
eqns(2,1) = v_g*X*Mw(1);
eqns(3,1) = v_succ*X*Mw(2);



end

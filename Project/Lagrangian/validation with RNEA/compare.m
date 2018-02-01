function compare(NB)
    x = (pi/2)*rand(NB,1);
    y = 100*rand(NB,1);
    
    model_LE = getModel(NB);
    model_LE.q = x;
    model_LE.qd = y;
    
    model_RNEA = getModel_RNEA(NB);
    model_RNEA.q = x;
    model_RNEA.qd = y;
    model_RNEA.n0 = model_LE.n0;
    model_RNEA.l0 = model_LE.l0;
    model_RNEA.g = [0 0 0]';
    H1 = RNEA(model_RNEA);
    H2 = RNEA_bd(model_RNEA);
    
    model_RNEA.qd = [0;0;0];
    model_RNEA.g = [0 9.80665 0]';
    C1 = RNEA(model_RNEA);
    C2 = RNEA_bd(model_RNEA);
    
    model_RNEA.g = [0 0 0]';
    I = eye(NB);
    for i=1:NB
        model_RNEA.qdd = I(:,i);
        D1(:,i) = RNEA(model_RNEA);
        D2(:,i) = RNEA_bd(model_RNEA);
    end
    
    [D,H,C] = LE_solver(model_LE);
    
    H,H1,H2
    C,C1,C2
    D,D1,D2
end
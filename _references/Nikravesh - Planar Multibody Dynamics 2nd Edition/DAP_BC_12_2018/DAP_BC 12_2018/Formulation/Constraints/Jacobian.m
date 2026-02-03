function D = Jacobian
    include_global

    D = zeros(nConst,nB3);
for Ji = 1:nJ 
    switch (Joints(Ji).type);
        case {'rev'}
            J_rev
        case {'tran'}
            J_tran
        case {'rev-rev'}
            J_rev_rev
        case {'rev-tran'}
            J_rev_tran
        case {'rigid'}
            J_rigid
        case {'disc'}
            J_disc
        case {'rel-rot'}
            J_rel_rot
        case {'rel-tran'}
            J_rel_tran
    end
        rs  = Joints(Ji).rows;
        re  = Joints(Ji).rowe;
    if Joints(Ji).iBindex ~= 0
        cis = Joints(Ji).colis;
        cie = Joints(Ji).colie;
        D(rs:re,cis:cie) = Di;
    end
    if Joints(Ji).jBindex ~= 0
        cjs = Joints(Ji).coljs;
        cje = Joints(Ji).colje;
        D(rs:re,cjs:cje) = Dj;
    end
end

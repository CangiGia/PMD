function rhs = RHSAcc(t)
    include_global

    rhs = zeros(nConst,1);
for Ji = 1:nJ   
    switch (Joints(Ji).type);
        case {'rev'}
            A_rev
        case {'tran'}
            A_tran
        case {'rev-rev'}
            A_rev_rev
        case {'rev-tran'}
            A_rev_tran
        case {'rigid'}
            A_rigid
        case {'disc'}
            A_disc
        case {'rel-rot'}
            A_rel_rot
        case {'rel-tran'}
            A_rel_tran
    end
        rs = Joints(Ji).rows;
        re = Joints(Ji).rowe;
        rhs(rs:re) = f;
end

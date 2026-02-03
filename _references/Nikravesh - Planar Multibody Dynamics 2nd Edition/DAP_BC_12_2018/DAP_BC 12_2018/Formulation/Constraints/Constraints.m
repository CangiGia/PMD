function phi = Constraints(t)
    include_global

    phi = zeros(nConst,1);
for Ji = 1:nJ
    switch (Joints(Ji).type);
        case {'rev'}
            C_rev
        case {'tran'}
            C_tran
        case {'rev-rev'}
            C_rev_rev
        case {'rev-tran'}
            C_rev_tran
        case {'rigid'}
            C_rigid
        case {'disc'}
            C_disc
        case {'rel-rot'}
            C_rel_rot
        case {'rel-tran'}
            C_rel_tran
    end
        rs = Joints(Ji).rows;
        re = Joints(Ji).rowe;
        phi(rs:re) = f;
end

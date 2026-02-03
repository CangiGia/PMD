function rhs = RHSVel(t)
    include_global

    rhs = zeros(nConst,1);
for Ji = 1:nJ   
    switch (Joints(Ji).type);
        case {'rel-rot'}
            V_rel_rot
            rhs(Joints(Ji).rows:Joints(Ji).rowe) = f;
        case {'rel-tran'}
            V_rel_tran
            rhs(Joints(Ji).rows:Joints(Ji).rowe) = f;
    end
end

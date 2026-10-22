function Contact(Ci, Pi, Bi, k, e, Mi)
    include_global

    pen = -Points(Pi).rP(2);
if pen > 0
    pen_d = -Points(Pi).rP_d(2);
    if flags(Ci) == 0
        pen_d0(Ci) = pen_d;
        flags(Ci) = 1;
    end
    if Mi == 1
        fy = Contact_LN(pen, pen_d, pen_d0(Ci), k, e); % penetration force
    else
        fy = Contact_FM(pen, pen_d, pen_d0(Ci), k, e); % penetration force
    end
    fsd = [0; fy];
    Bodies(Bi).f = Bodies(Bi).f + fsd;
    Bodies(Bi).n = Bodies(Bi).n + Points(Pi).sP_r'*fsd;
else
    flags(Ci) = 0;
end
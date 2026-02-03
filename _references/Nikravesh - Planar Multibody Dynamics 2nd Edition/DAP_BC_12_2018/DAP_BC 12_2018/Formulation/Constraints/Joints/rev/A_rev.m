% revolute_A
% r-h-s of acc. constraints for a revolute joint

Pi = Joints(Ji).iPindex;     Pj = Joints(Ji).jPindex;
Bi = Points(Pi).Bindex;      Bj = Points(Pj).Bindex;

if Bi == 0
    f = s_rot(Points(Pj).sP_d)*Bodies(Bj).p_d;
elseif Bj == 0
    f = -s_rot(Points(Pi).sP_d)*Bodies(Bi).p_d;
else
    f = -s_rot(Points(Pi).sP_d)*Bodies(Bi).p_d + ...
         s_rot(Points(Pj).sP_d)*Bodies(Bj).p_d;
end

    if Joints(Ji).fix == 1
       f = [f
            0];
    end
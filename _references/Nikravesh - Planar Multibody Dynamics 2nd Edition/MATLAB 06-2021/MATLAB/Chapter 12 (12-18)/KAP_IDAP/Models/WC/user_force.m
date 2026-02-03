% ===== user supplied force ======
% ---- Web_cutter ----
% ***** If there are no forces, this file should be provided but left blank ******

    fcut = [0; 0];
    del = RP2(i,2) - RP3(i,2);
    if del > 0 && del < 0.005
        del_d = VP2(i,2) - VP3(i,2);
        if del_d < 0
            fcut = [0; 1000];
        end
    end
    
    h_a = h_a + ...
         [0
          0
          0
          fcut
          sP2r'*fcut
         -fcut
         -sP3r'*fcut];

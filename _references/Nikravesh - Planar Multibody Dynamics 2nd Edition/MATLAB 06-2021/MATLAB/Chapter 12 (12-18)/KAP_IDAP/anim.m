% Animate response
    
    nt = size(T,1);
    pause on
%%    
for i = 1:nt
    R_P1 = RP1(i,:); R_P2 = RP2(i,:); R_P3 = RP3(i,:);
    R_A = RA(i,:); R_B = RB(i,:);
    plot_system,     drawnow;         
    if T(i) == 0 && nt > 1
        disp('  Press a key to continue!')
        pause 
    end
    pause off
end
%%  Plot paths of points C and P
    hold on


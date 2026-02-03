% Animate response
        
    circle = zeros(2,40);
    [cx cy cz] = cylinder(R, 40); circle = [cx(1,:); cy(1,:)];
    O1 = [R2; 0]; O2 = [0; 0];
% Interpolate output for constant timesteps
    uTTi = interp1q(TT, uTT, (0:0.02:Tfinal)');
    nt = size(uTTi,1);
    
for i = 1:nt
    r1 = uTTi(i,1:2)'; r2 = uTTi(i,4:5)';
% Plot points 
    plot(O1(1),O1(2),'ro')
    hold on
    plot(O2(1),O2(2),'go')
% Draw lines   
    line([r1(1),O1(1)], [r1(2),O1(2)],'color','r')
    line([r2(1),O2(1)], [r2(2),O2(2)],'color','g')  
% Draw circles
    xx = r1(1) + circle(1,:); 
    yy = r1(2) + circle(2,:);
    for j = 1:40
        line([xx(j),xx(j+1)], [yy(j),yy(j+1)],'color','r')
    end
    xx = r2(1) + circle(1,:); 
    yy = r2(2) + circle(2,:);
    for j = 1:40
        line([xx(j),xx(j+1)], [yy(j),yy(j+1)],'color','g')
    end
% Set axes
    grid on
    axis manual
    axis equal
    axis([-0.12 0.14 -0.14 0.12]);
    hold off
    drawnow         
end

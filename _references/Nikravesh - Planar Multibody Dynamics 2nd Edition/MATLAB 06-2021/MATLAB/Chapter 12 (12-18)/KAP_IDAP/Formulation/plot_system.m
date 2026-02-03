function plot_system
%  Animation
global RQ R_A R_B R_P1 R_P2 R_P3 P1_yes P2_yes P3_yes axis_limits

% Plot points 
    plot(0,0,'ko','MarkerFaceColor','k','MarkerSize',8)
    hold on
    plot(RQ(1),RQ(2),'ko','MarkerFaceColor','k','MarkerSize',8)
    plot(R_A(1),R_A(2),'ko','color','k','MarkerSize',8)
    plot(R_B(1),R_B(2),'ko','color','k','MarkerSize',8)
% Include text
    text(0,0,'  O')
    text(RQ(1),RQ(2),'  Q')
    text(R_A(1),R_A(2),'  A')
    text(R_B(1),R_B(2),'  B')
% Draw lines between points
    line([0,R_A(1)],[0,R_A(2)],'color','r','LineWidth',2)
    line([R_A(1),R_B(1)],[R_A(2),R_B(2)],'color','b','LineWidth',2)
    line([RQ(1),R_B(1)],[RQ(2),R_B(2)],'color','g','LineWidth',2)
% Point P1 (crank point)
    if P1_yes == 1
        plot(R_P1(1),R_P1(2),'ko','MarkerFaceColor','b','MarkerSize',5)
        text(R_P1(1),R_P1(2),'  P1')
        line([R_P1(1),R_A(1)],[R_P1(2),R_A(2)],'color','b','LineWidth',1)
        line([R_P1(1),0],[R_P1(2),0],'color','b','LineWidth',1)    
    end
% Point P2 (coupler point)
    if P2_yes == 1
        plot(R_P2(1),R_P2(2),'ko','MarkerFaceColor','b','MarkerSize',5)
        text(R_P2(1),R_P2(2),'  P2')
        line([R_P2(1),R_A(1)],[R_P2(2),R_A(2)],'color','b','LineWidth',1)
        line([R_P2(1),R_B(1)],[R_P2(2),R_B(2)],'color','b','LineWidth',1)    
    end
% Point P3 (follower point)
    if P3_yes == 1
        plot(R_P3(1),R_P3(2),'ko','MarkerFaceColor','g','MarkerSize',5)
        text(R_P3(1),R_P3(2),'  P3')
        line([R_P3(1),RQ(1)],[R_P3(2),RQ(2)],'color','g','LineWidth',1)
        line([R_P3(1),R_B(1)],[R_P3(2),R_B(2)],'color','g','LineWidth',1)
    end
% Set axes
    axis manual; axis equal; axis(axis_limits);
    hold off
function plot_system
% 2D Animation

include_global

figure(1)  

% Plot body center points
for Bi = 1:nB
    plot(Bodies(Bi).r(1),Bodies(Bi).r(2),'ko', ...
        'MarkerFaceColor',Bodies(Bi).color,'MarkerSize',5)
    hold on
    text(Bodies(Bi).r(1),Bodies(Bi).r(2),sprintf('   (%.0f',Bi))
    text(Bodies(Bi).r(1),Bodies(Bi).r(2),sprintf('       )'))
end   

% Draw lines between body centers and points on those bodies
for Bi = 1:nB   
    npts = length(Bodies(Bi).pts);
    linecolor = Bodies(Bi).color;
    for j = 1:npts
        line([Bodies(Bi).r(1),Points(Bodies(Bi).pts(j)).rP(1)], ...
            [Bodies(Bi).r(2),Points(Bodies(Bi).pts(j)).rP(2)], ...
            'color',linecolor,'LineWidth',1)
    end
end

% Plot points that are defined by 's' vectors
for i = 1:nP
    plot(Points(i).rP(1),Points(i).rP(2),'ko', ...
        'MarkerFaceColor','k','MarkerSize',2)
end

% Draw lines between points that are connected by springs
for i = 1:nF
    switch (Forces(i).type);
        case {'ptp'}
            pt1 = Forces(i).iPindex;  pt2 = Forces(i).jPindex;
            line([Points(pt1).rP(1),Points(pt2).rP(1)], ...
             [Points(pt1).rP(2),Points(pt2).rP(2)], ...
             'color','m','LineStyle','--')
    end
end

% for Ji = 1:nJ 
%     switch (Joints(Ji).type);
% %         case {'rev'}
% %             Pi = Joints(Ji).iPindex;
% %             plot(Points(Pi).rP(1),Points(Pi).rP(2),'ko', ...
% %                 'MarkerFaceColor','k','MarkerSize',4)
% %         case {'tran'}
% %             J_tran
%         case {'rev-rev'}
%             Pi = Joints(Ji).iPindex; Pj = Joints(Ji).jPindex;
% %             plot(Points(Pi).rP(1),Points(Pi).rP(2),'ko', ...
% %                 'MarkerFaceColor','k','MarkerSize',4)
% %             plot(Points(Pj).rP(1),Points(Pj).rP(2),'ko', ...
% %                 'MarkerFaceColor','k','MarkerSize',4)
%             line([Points(Pi).rP(1),Points(Pj).rP(1)], ...
%                 [Points(Pi).rP(2),Points(Pj).rP(2)], 'color','c')
%        case {'rev-tran'}
%             Pi = Joints(Ji).iPindex;
%             plot(Points(Pi).rP(1),Points(Pi).rP(2),'ko', ...
%                 'MarkerFaceColor','k','MarkerSize',4)
%     end
% end

for Bi = 1:nB
    linecolor = Bodies(Bi).color;
    switch (Bodies(Bi).shape);
        case {'circle'}
            xx = Bodies(Bi).r(1) + Bodies(Bi).circ(1,:); 
            yy = Bodies(Bi).r(2) + Bodies(Bi).circ(2,:);
            for i = 1:40
                line([xx(i),xx(i+1)], [yy(i),yy(i+1)],'color',linecolor)
            end
        case {'rect'}
            P5 = zeros(2,5);
            for i = 1:4
                P5(:,i) = Bodies(Bi).r + Bodies(Bi).A*Bodies(Bi).P4(:,i);
            end
            P5(:,5) = P5(:,1);
            for i = 1:4
                line([P5(1,i),P5(1,i+1)], [P5(2,i),P5(2,i+1)],'color',linecolor)
            end
        case {'line'}
            P5 = zeros(2,2);
            for i = 1:2
                P5(:,i) = Bodies(Bi).r + Bodies(Bi).A*Bodies(Bi).P4(:,i);
            end
            line([P5(1,1),P5(1,2)], [P5(2,1),P5(2,2)],'color',linecolor)
    end
end

    grid on
% Set axes
    axis manual
    axis equal
    axis([xmin xmax ymin ymax]);
% axis square
%     xlabel('X');    ylabel('Y');
 
    hold off

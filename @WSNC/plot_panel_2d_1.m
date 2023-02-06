% »æÖÆ¶þÎ¬¸²¸ÇÍ¼
function plot_panel_2d_1(obj, agent, alg_name)
    point = obj.trans_agent(agent);
    N = size(point,1);
    cmap = hsv(N); %// define colors. You could change `hsv` to `jet`, `cool`, ...
    alpha = .5; %// define level of transparency
    t = 0 : .1 : 2 * pi;
    N = size(point,1);
    figure;
    hold on
    for i=1:obj.L
        for j=1:obj.L
            plot(i-obj.data/2,j-obj.data/2, '.', 'Color', 'k', 'MarkerSize', 2)
        end
    end
    for i=1:N
        x = obj.R * cos(t) + point(i,1);
        y = obj.R * sin(t) + point(i,2);
        patch(x, y, cmap(i,:), 'facealpha', alpha, 'edgecolor', cmap(i,:));
    end
    set(gca, 'XLim', [0,obj.L]);
    % set(gca, 'YLim', [0,obj.L]);
    if exist('alg_name','var')
        title(alg_name,'FontSize',20);
    end
    set(gca,'FontSize',16);
    axis equal;
end
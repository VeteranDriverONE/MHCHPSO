% ������ά����ͼ
function plot_panel_3d(obj, points, alg_name)
    N = size(points,1);
    figure;
    hold on;
    for i=1:N
        [x,y,z] = sphere(20);
        % �����뾶
        x = obj.R * x;
        y = obj.R * y;
        z = obj.R * z;
        % ��������
        x = x + points(i,1);
        y = y + points(i,2);
        z = z + points(i,3);
        % ʹ��mesh����
        % mesh(x,y,z);
        % ʹ��surf����
        surf(x,y,z);
    end
    set(gca, 'XLim', [0,obj.L]);
    if exist('alg_name','var')
        title(alg_name,'FontSize',20);
    end
    set(gca,'FontSize',16);
    axis equal;
end
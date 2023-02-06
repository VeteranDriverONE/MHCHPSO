% 计算3维环境下适应度
function [p] = wsn3d(obj, point)
    [N, ~] = size(point);                      % 节点总个数
    [X, Y, Z] = meshgrid(0:obj.data:obj.L,0:obj.data:obj.L,0:obj.data:obj.L);        % 离散化区域内的点
    [x_l,y_l, z_l] = size(X);
    for i = 1:N
        D = sqrt((point(i,1)-X).^2 + (point(i,2)-Y).^2 + (point(i,3)-Z).^2);
        [x0, y0, z0] = find(D <= obj.R);             % 检测出圆覆盖点的坐标
        Ind = (z0-1)*(x_l*y_l)+(x0-1)*y_l+y0;                % 坐标与索引转化
        M(Ind) = 1;                          % 改变覆盖状态
    end
    p = sum(M(1:end))/(x_l*y_l*z_l);         % 计算覆盖比例
    p = 1 - p;
end
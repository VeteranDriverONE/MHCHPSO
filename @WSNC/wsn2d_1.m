function [p] = wsn2d_1(obj, point)
    point = obj.trans_agent(point);
    N = size(point,1);  
    [m, n] = meshgrid(1:obj.data:obj.L);        % 离散化区域内的点
    m = m - obj.data/2;
    n = n - obj.data/2;
    [row, col] = size(m);
    M = zeros(row,col);
    for i = 1:N
        D = sqrt((m-point(i,1)).^2+(n-point(i,2)).^2);   % 计算坐标点到圆心的距离
        [m0, n0] = find(D <= obj.R);             % 检测出圆覆盖点的坐标
        M = M + (D<=obj.R)*1;
    end
    scale = sum(sum(M>0))/(row*col);         % 计算覆盖比例
    p = 1 - scale;
end

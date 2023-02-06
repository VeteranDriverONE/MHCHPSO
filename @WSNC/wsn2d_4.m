% 2维，无规则形状, 反Z字形
function [p] = wsn2d_4(obj, point)
    point = obj.trans_agent(point);
    N = size(point,1);
    [col, row] = meshgrid(1:obj.data:obj.L);
    row = row - obj.data/2;
    col = col - obj.data/2;
    M = zeros(size(col));
    for i = 1:N
        D = sqrt((row-point(i,1)).^2+(col-point(i,2)).^2);   % 计算坐标点到圆心的距离
        D = obj.ava_2d_3.*D + ~obj.ava_2d_3.*(obj.R+1);
        M = M + (D <= obj.R);             % 检测出圆覆盖点的坐标
    end
    scale = sum(sum(M>0))/sum(sum(obj.ava_2d_1));         % 计算覆盖比例
    p = 1- scale;
end
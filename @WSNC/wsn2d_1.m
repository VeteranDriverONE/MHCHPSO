function [p] = wsn2d_1(obj, point)
    point = obj.trans_agent(point);
    N = size(point,1);  
    [m, n] = meshgrid(1:obj.data:obj.L);        % ��ɢ�������ڵĵ�
    m = m - obj.data/2;
    n = n - obj.data/2;
    [row, col] = size(m);
    M = zeros(row,col);
    for i = 1:N
        D = sqrt((m-point(i,1)).^2+(n-point(i,2)).^2);   % ��������㵽Բ�ĵľ���
        [m0, n0] = find(D <= obj.R);             % ����Բ���ǵ������
        M = M + (D<=obj.R)*1;
    end
    scale = sum(sum(M>0))/(row*col);         % ���㸲�Ǳ���
    p = 1 - scale;
end

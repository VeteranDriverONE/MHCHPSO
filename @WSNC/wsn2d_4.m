% 2ά���޹�����״, ��Z����
function [p] = wsn2d_4(obj, point)
    point = obj.trans_agent(point);
    N = size(point,1);
    [col, row] = meshgrid(1:obj.data:obj.L);
    row = row - obj.data/2;
    col = col - obj.data/2;
    M = zeros(size(col));
    for i = 1:N
        D = sqrt((row-point(i,1)).^2+(col-point(i,2)).^2);   % ��������㵽Բ�ĵľ���
        D = obj.ava_2d_3.*D + ~obj.ava_2d_3.*(obj.R+1);
        M = M + (D <= obj.R);             % ����Բ���ǵ������
    end
    scale = sum(sum(M>0))/sum(sum(obj.ava_2d_1));         % ���㸲�Ǳ���
    p = 1- scale;
end
% ����3ά��������Ӧ��
function [p] = wsn3d(obj, point)
    [N, ~] = size(point);                      % �ڵ��ܸ���
    [X, Y, Z] = meshgrid(0:obj.data:obj.L,0:obj.data:obj.L,0:obj.data:obj.L);        % ��ɢ�������ڵĵ�
    [x_l,y_l, z_l] = size(X);
    for i = 1:N
        D = sqrt((point(i,1)-X).^2 + (point(i,2)-Y).^2 + (point(i,3)-Z).^2);
        [x0, y0, z0] = find(D <= obj.R);             % ����Բ���ǵ������
        Ind = (z0-1)*(x_l*y_l)+(x0-1)*y_l+y0;                % ����������ת��
        M(Ind) = 1;                          % �ı串��״̬
    end
    p = sum(M(1:end))/(x_l*y_l*z_l);         % ���㸲�Ǳ���
    p = 1 - p;
end
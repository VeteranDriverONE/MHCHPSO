% 构造函数
classdef WSNC
    properties
        L
        R
        data
        dim
        point_num
        ava_2d_1
        ava_2d_2
        ava_2d_3
    end
    methods
        % 构造函数
        function obj = WSNC(L,R,data,dim,point_num)
            obj.L=L;
            obj.R=R;
            obj.data=data;
            obj.dim=dim;
            obj.point_num=point_num;

            load('ava_2d_1'); % 加载Z字形
            obj.ava_2d_1 = ava_2d_1;
% 生成Z字形
%             obj.ava_2d_1 = zeros(40,40);
%             for i=1:40
%                 for j=1:40
%                     if i<=10
%                         if j<=20 && j>0
%                             obj.ava_2d_1(i,j) = 1;
%                         end
%                     elseif i<=30
%                         if j<=i+10 && j>i-10
%                             obj.ava_2d_1(i,j) = 1;
%                         end
%                     elseif i<=40
%                         if j<=40 && j>20
%                             obj.ava_2d_1(i,j) = 1;
%                         end
%                     else
%                         disp('Other')
%                     end
%                 end
%             end
            load('ava_2d_2');  % 加载菱形字形
            obj.ava_2d_2 = ava_2d_2;
% 生成菱形
%             obj.ava_2d_2 = zeros(79,79);
%             for i = 1:79
%                 index = 40-abs(40-i);
%                 len = 2 * index-1;
%                 s = 40 - floor(len / 2);
%                 e = 40 + floor(len / 2);
%                 for j = s:e
%                     obj.ava_2d_2(i,j) = 1;
%                 end
%             end
%             for i = 1:39
%                 index = 20-abs(20-i);
%                 len = 2 * index-1;
%                 s = 20 - floor(len / 2)+20;
%                 e = 20 + floor(len / 2)+20;
%                 for j = s:e
%                     obj.ava_2d_2(i+20,j) = 0;
%                 end
%             end
            load('ava_2d_3');  % 加载反Z字形
            obj.ava_2d_3 = ava_2d_3;
%             ava_2d_3 = zeros(size(ava_2d_1));
%             for i=1:size(ava_2d_1,1)
%                 for j=1:size(ava_2d_1,2)
%                     ava_2d_3(i,j) = ava_2d_1(i,size(ava_2d_1,2)+1-j);
%                 end
%             end
            disp('WSNC对象加载完成')
        end
    end
end
        
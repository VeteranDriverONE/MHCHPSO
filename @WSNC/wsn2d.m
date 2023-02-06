function [p] = wsn2d(obj, point, env_id)
    point = point'; % 兼容CEC格式
    if env_id == 1
        p = obj.wsn2d_1(point);
    elseif env_id == 2
        p = obj.wsn2d_2(point);
    elseif env_id == 3
        p = obj.wsn2d_3(point);
    elseif env_id == 4
        p = obj.wsn2d_4(point);
    else
        disp('other')
    end
end

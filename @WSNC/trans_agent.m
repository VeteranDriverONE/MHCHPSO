function new_agent = trans_agent(obj,agent)
    if obj.dim==2
        x = agent(1:obj.point_num);
        y = agent(obj.point_num+1:end);
        new_agent = [x',y'];
    elseif obj.dim==3
        x = agent(1:obj.point_num);
        y = agent(obj.point_num+1:obj.point_num*2);
        z = agent(obj.point_num*2+1:end);
        new_agent = [x',y',z'];
    end
end
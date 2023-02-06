% º∆À„  ”¶∂»
function [fit] = call_fit(obj, agent, func_id)
    agent = obj.trans_agent(agent);
    if func_id == 1
        fit = obj.wsn2d(agent);
    elseif func_id == 0
        fit = obj.wsn3d(agent);
    elseif func_id == 2
        fit = obj.wsn_2d_1(agent);
    end
end
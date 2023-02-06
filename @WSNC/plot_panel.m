% ªÊ÷∆∏≤∏«Õº
function plot_panel(obj, agent, env_id, alg_name)
    points = obj.trans_agent(agent);
    if env_id == 0
        obj.plot_panel_3d(points, alg_name)
    elseif env_id == 1
        obj.plot_panel_2d_1(points, alg_name)
    elseif env_id == 2
        obj.plot_panel_2d_2(points, alg_name)
    elseif env_id == 3
        obj.plot_panel_2d_3(points, alg_name)
    elseif env_id == 4
        obj.plot_panel_2d_4(points, alg_name)
    end
end
Base.@propagate_inbounds function interp_o2_at_half(us)
    (1 / 2) * (us[1] + us[2])
end

Base.@propagate_inbounds function interp_o4_at_half(us)
    (1 / 16) * (-us[1] + 9 * us[2] + 9 * us[3] - us[4])
end

Base.@propagate_inbounds function diff1_o2_at_half(us)
    (-us[1] + us[2])
end

Base.@propagate_inbounds function diff1_o4_at_half(us)
    (1 / 24) * (us[1] - 27us[2] + 27us[3] - us[4])
end

Base.@propagate_inbounds function diff2_o4_at_node(us)
    (1 / 12) * (-us[1] + 16 * us[2] - 30 * us[3] + 16 * us[4] - us[5])
end

Base.@propagate_inbounds function diff1_o4_at_node(us)
    (1 / 12) * (us[1] - 8us[2] + 8us[4] - us[5])
end

Base.@propagate_inbounds function diff1_o2_at_node(us)
    (1 / 2) * (us[3] - us[1])
end

Base.@propagate_inbounds function sum_flux_o4(j1s, j3s)
    (9 / 8) * (-j1s[1] + j1s[2]) - (1 / 8) * (1 / 3) * (-j3s[1] + j3s[2])
end

Base.@propagate_inbounds function diff2_o2_at_node_cross(us)
    (1 / 4) * (us[3, 3] - us[3, 1] - us[1, 3] + us[1, 1])
end

Base.@propagate_inbounds function diff2_o2_at_node_diag(us)
    us[3] - 2 * us[2] + us[1]
end

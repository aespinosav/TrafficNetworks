#This file will have functions for plotting TA flows of the networks
#and hopefully graph visualisations as well (nothing too ambitious... ha)
"""
Converts the flows of rn for the chosen regime (UE/SO) to a 
data frame, mostly for plotting, but could be useful for storing 
the data. (maybe this function should be defined in another file?)
"""
function make_data_frame_from_flows(rn::RoadNetwork, regime="UE")
    if regime=="UE"
        flows = rn.flows_ue
    else
        flows = rn.flows_so
    end
    
    names = append!([:q],[symbol("x$(i)") for i in 1:size(flows)[1]])
    data_frame = DataFrame(cat(2, rn.demand_range, flows'))
    names!(data_frame, names)
end

function plot_flows(rn::RoadNetwork, regime="UE")
    data_flows = make_data_frame_from_flows(rn, regime)
    #layers = [ layer( x=data_flows[1], y=data_flows[i], Geom.line, Theme(defaul_color=distinguishable_colors(size(data_flows)[2])[i])) for i in 2:length(flows) ]
    plot(melt(data_flows, :q), 
         x=:q, 
         y=:value, 
         color=:variable,
         Geom.line,
         Guide.XLabel("Demand"),
         Guide.YLabel("Flow"),
         Guide.colorkey("Flows"),
         Coord.Cartesian(ymin=0, ymax=rn.demand_range[end]),
         Guide.Title("Edge flows") )
end

# function plot_normalised_flows(rn::RoadNetwork, regime="UE")
#     data_flows = make_data_frame_from_flows(rn, regime)
#     #layers = [ layer( x=data_flows[1], y=data_flows[i], Geom.line, Theme(defaul_color=distinguishable_colors(size(data_flows)[2])[i])) for i in 2:length(flows) ]
#     plot(melt(data_flows, :q), 
#          x=:q, 
#          y=:value, 
#          color=:variable,
#          Geom.line,
#          Guide.XLabel("Demand"),
#          Guide.YLabel("Flow"),
#          Guide.colorkey("Flows"),
#          Coord.Cartesian(ymin=0, ymax=rn.demand_range[end]),
#          Guide.Title("Normalised edge flows") )
# end
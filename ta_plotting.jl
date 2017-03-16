#This file will have functions for plotting TA flows of the networks
#and hopefully graph visualisations as well (nothing too ambitious... ha)
"""
Converts the flows of rn for the chosen regime (UE/SO) to a 
data frame, mostly for plotting, but could be useful for storing 
the data. (maybe this function should be defined in another file?)
"""
function flows_data_frame(rn::RoadNetwork, regime="UE")
    if regime=="UE"
        flows = rn.flows_ue
    else
        flows = rn.flows_so
    end
    
    names = append!([:q],[symbol("x$(i)") for i in 1:size(flows)[1]])
    data_frame = DataFrame(cat(2, rn.demand_range, flows'))
    names!(data_frame, names)
end

"""
Converts a solution array from 'ta_solve' into a data frame. Array has to be  given as 'flows' to a data frame.
Needs the demand range as well...
"""
function flows_data_frame(flows::Array{Float64,2}, demand_range::Array{Float64,1})
    names = append!([:q],[symbol("x$(i)") for i in 1:size(flows)[1]])
    data_frame = DataFrame(cat(2, demand_range, flows'))
    names!(data_frame, names)
    data_frame
end


"""
normalises the flows of rn with respect to total demand for the chosen regime (UE/SO) 
and makes a data frame.

Be careful about divisions by zero (it is a waste of processor to actually optimise for 0 demand...)
"""
function normflows_data_frame(rn::RoadNetwork, regime="UE")
    flows = flows_data_frame(rn, regime)[1:end, :]
    for i in 2:size(flows)[2]
        flows[i] = flows[i] ./ flows[:q]
    end
    flows
end

"""
Plots the edge flows of a RoadNetwork that has been solved (with appropriate regime UE/SO).
Uses Gadfly.
"""
function plot_flows(rn::RoadNetwork, regime="UE")
    data_flows = flows_data_frame(rn, regime)
    #layers = [ layer( x=data_flows[1], y=data_flows[i], Geom.line, Theme(defaul_color=distinguishable_colors(size(data_flows)[2])[i])) for i in 2:length(flows) ]

    melted_data = melt(data_flows, :q)

    plot(melted_data, 
         x=:q, 
         y=:value, 
         color=:variable,
         Geom.line,
         Guide.XLabel("Demand"),
         Guide.YLabel("Flow"),
         Guide.colorkey("Flows"),
         Coord.Cartesian(ymin=0, ymax=maximum(melted_data[:value])),
         Guide.Title("Edge flows") )
end

"""
Plots the normalised edge flows.
"""
function plot_normalised_flows(rn::RoadNetwork, regime="UE")
    data_flows = normflows_data_frame(rn, regime)
    #layers = [ layer( x=data_flows[1], y=data_flows[i], Geom.line, Theme(defaul_color=distinguishable_colors(size(data_flows)[2])[i])) for i in 2:length(flows) ]
    plot(melt(data_flows, :q), 
         x=:q, 
         y=:value, 
         color=:variable,
         Geom.line,
         Guide.XLabel("Demand"),
         Guide.YLabel("Flow (normalised)"),
         Guide.colorkey("Flows"),
         Coord.Cartesian(ymin=0, ymax=1),
         Guide.Title("Normalised edge flows") )
end

#this file contains the definition for the TrafficAssignment type
# For now this is only on the assignment-branch of the repo, once everything
# is working properly, will merge into the master branch.
#
using DataFrames

"""
Traffic assignment type, that contains the network (I dont really want rn to be a
deep copy though) and will have the flows for a given traffic assignment on all links
for a given demand range.

It has functions to construct data frames. This way it doesnt store these things in 
memory and only provides them when needed. Most likely when saving simulation data
or doing explratory data analysis.

Having this type is justified, to make the whole process make more sense and it makes
all other functions more modular instead of having a monster type that stores
everything but is restricted because of the way it was coded originally for only 
UE or SO. This makes the code more flexible and reusable... If i get it right that is.
"""
type TrafficAssignment
    rn::RoadNetwork
    demand_range::Array{Float64,1} # Array compattible with rn's OD matrix
    flows::DataFrame# Flow assignment over the demand_range (will be dataframe)
    regime::AbstractString #Some label (metadata) on the regime being solved
end

"""
Constructor for TrafficAssignment. Makes an empty dataframe with a column for 
demand and a column for each of the links to hold the flows.

Handles the size of the network automatically.
"""
function TrafficAssignment(rn::RoadNetwork, demand_range::Array{Float64,1}, flow_label="x", regime="custom TA regime")
    m = num_edges(rn)
    flow_col_labels = [flow_label*"$i" for i in 1:m]

    ex = "flows = DataFrame("
    ex *= "q = Float64[]"
    for i in 1:length(flow_col_labels)
        ex *= ", $(flow_col_labels[i]) = Float64[]"
    end
    ex *= ")"
    ex = parse(ex)
    eval(ex)
    
    TrafficAssignment(rn, demand_range, flows, regime)
end
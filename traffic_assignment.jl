#this file contains the definition for the TrafficAssignment type
# For now this is only on the assignment-branch of the repo, once everything
# is working properly, will merge into the master branch.
#

type TrafficAssignment
    rn::RoadNetwork
    demand_range # Array compattible with rn's OD matrix
    flow # Flow assignment over the demand_range
    regime::AbstractString
end
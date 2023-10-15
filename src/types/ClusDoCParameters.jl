struct ClusterParameters
    epsilon::Float64               # default 20
    minpoints::Int                 # default 3
    uselocalradius_threshold::Bool # default true
    smoothingradius::Int           # default 15
    minsigclusterpoints::Int          # default 10
end

struct DoCParameters
    localradius::Float64               # default 20
    radiusmax::Float64                 # default 500
    radiusstep::Float64                # default 10
    colocalized_threshold::Float64     # default 0.4
end

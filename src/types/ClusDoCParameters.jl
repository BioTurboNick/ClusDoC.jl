struct ClusterParameters
    epsilon::Float64               # default 20
    minpoints::Int                 # default 3
    uselocalradius_threshold::Bool # default true
    smoothingradius::Int           # default 15
    pointscutoff::Int                    # default 10 ADD TO UI
end

struct DoCParameters
    localradius::Float64               # default 20
    radiusmax::Float64                 # default 500
    radiusstep::Float64                # default 10
end

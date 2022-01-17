struct ClusDoCParameters
    doc_localradius::Float64               # default 20
    doc_radiusmax::Float64                 # default 500
    doc_radiusstep::Float64                # default 10
    cluster_epsilon::Float64               # default 20
    cluster_minpoints::Int                 # default 3
    cluster_uselocalradius_threshold::Bool # default true
    cluster_smoothingradius::Float64       # default 15
end

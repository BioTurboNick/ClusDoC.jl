using ClusDoC
using Test

@testset "Lr" begin
    channel1points = [(0.0, 0.0)]
    channel2pointsbase = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]
    lr_radius = 40

    for i ∈ 1:10:41
        channel2points = channel2pointsbase .* i
        lr(channel2points, channel1points, lr_radius)
    end
end


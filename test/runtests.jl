using ClusDoC
using Test

@testset "Lr" begin
    channel1pointsX = [28739.6, 29635, 28894.8, 29151.3, 28824.9, 29491.4, 29129.6, 29131.9, 28574.5, 29492.4, 29132.5,
                       28580.6, 28561.3, 28577.5, 28739.6, 29633.5, 29131.6]
    channel1pointsY = [37181.4, 37849.4, 37446.7, 37957, 37026.1, 37287.2, 37786.3, 37959, 37558.8, 37309.7, 37771.8,
                       37555, 37571.5, 37568.2, 37194.2, 37851.2, 37761]
    coordinates = permutedims([channel1pointsX channel1pointsY])
    lr_radius = 20
    roiarea = 1373584

    expected_density = [0.0015915494309189533, 0.0015915494309189533, 0.0007957747154594767, 0.0015915494309189533, 0.0007957747154594767,
                        0.0007957747154594767, 0.0015915494309189533, 0.0015915494309189533, 0.0031830988618379067, 0.0007957747154594767,
                        0.00238732414637843, 0.00238732414637843, 0.00238732414637843, 0.0031830988618379067, 0.0015915494309189533,
                        0.0015915494309189533, 0.0015915494309189533]

    expected_lr = [193740.44623196532, 193740.44623196532, 0.0, 193740.44623196532, 0.0, 0.0, 193740.44623196532,
                   193740.44623196532, 335568.2963548302, 0.0, 273990.36664146074, 273990.36664146074, 273990.36664146074, 
                   335568.2963548302, 193740.44623196532, 193740.44623196532, 193740.44623196532]

    Lr, density = lr(coordinates, lr_radius, roiarea)

    @test expected_lr ≈ Lr
    @test expected_density ≈ density
end

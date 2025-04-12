using Droplets
using Test

@testset begin
    # export esat,sat,drdtcondensation1,drdtcondensation2,drdtcondensation3
    # export FK,FD,drkohler,θcondenseupdate!,qvcondenseupdate!,dXkohler_function_of_radius
    # export dXkohler_function_of_radius_activated,drkohler_activated

    T = 273 # Temperature in K
    S = 1.01 #Saturation
    smallS = 0.95 #Saturation < Sactivation for these settings
    M = 1e-16 #Mass of solute
    a = 3.3*10^(-7)/T #(m)
    m = 58.44 #molecular weight of NaCL
    b = 4.3 *2/m/1e6 #m^3 for NaCL
    denom = FK(T)+FD(T) #denominator for the Köhler equation
    radius = 2e-7
    @test esat(T) ≈ 610.94
    @test esat(280) > esat(273)

    @test drkohler(radius,M,m,T,S,1) > 0.0

    qv = 0.22*1e-3
    P = 101300
    sat(qv,P)
    # a, b, S, M,denom = p #1
    # M,m,T,qv,P = p #2
    # M,m,T,Senv = p#

    # drkohler(R,M,m,T,qv,P,timestep)
    # drkohler(R,M,m,T,Senv,timestep)

    
end
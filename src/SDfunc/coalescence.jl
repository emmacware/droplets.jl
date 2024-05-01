#---------------------------------------------------------
#Collision Coalescence Functions
#---------------------------------------------------------
#currently, this collision-coalescence is only set up to handle radius, 
#multiplicity, volume, and solute mass as superdroplet attributes


#functions:

#INITIALIZATION
    #init_ξ_const(Ns,ΔV,n0,R0),         returns ξstart, Rstart, Xstart
    #init_ξ_const(Ns,ΔV,n0,R0,M0),      returns ξstart, Rstart, Xstart, Mstart

#KERNEL FUNCTIONS
    #terminal_v(r),                     returns tv
    #collision_efficiency(R1,R2),       returns E
    #hydrodynamic(droplet1,droplet2),   returns E*π*Rsum^2*vdif
    #hydrodynamic(R1,R2,X,Xs),          returns E .*π .* Rsum^2 .*vdif
    #golovin(droplet1,droplet2),        returns b*Xsum
    #golovin(R1,R2,X,Xs),               returns b*Xsum

#SDM FUNCTIONS
    #calc_Ps(droplet1,droplet2,Δt,ΔV;kernel=golovin),                       returns Ps
    #calc_Ps(R1,R2,X,Xs,Δt,ΔV,ξj,ξk;kernel=golovin),                        returns Ps
    #coalescence_timestep!(droplets,Ns,Δt,ΔV;kernel=golovin),               returns droplets
    #coalescence_timestep!(ξ,R,X,Ns,Δt,ΔV;kernel=golovin),                  returns ξ,R,X
    #coalescence_timestep!(ξ,R,M,X,Ns,Δt,ΔV;kernel=golovin),                returns ξ,R,M,X
    #all_or_nothing!(R1,R2,X1,X2,ξ1,ξ2,M1,M2,scale,Δt,ΔV;kernel=golovin),   returns R1,R2,X1,X2,ξ1,ξ2,M1,M2
    #all_or_nothing!(R1,R2,X1,X2,ξ1,ξ2,scale,Δt,ΔV;kernel=golovin),         returns R1,R2,X1,X2,ξ1,ξ2
    #all_or_nothing!(droplet1,droplet2,scale,Δt,ΔV;kernel=golovin),         returns droplet1,droplet2


#SMALL ALPHA STUDY
    #coalescence_timestep_small_alpha!(droplets,Ns,y,Δt,ΔV;kernel=golovin), returns droplets

#UNIT TEST RUN, RETURNS GRAPH 0:1200:3600 seconds
    #coalescence_unittest_graph!(droplets,Ns,Δt,ΔV;smooth=true,label=true,kernel=golovin,
        #radius_bins_edges=10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128))
    #coalescence_unittest_graph!(ξstart,Rstart,Xstart,Ns,Δt,ΔV,smooth=smooth,label=true,kernel=golovin,
        #radius_bins_edges=10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128))


#---------------------------------------------------------
# INITIALIZATION
#---------------------------------------------------------
# init_ξ_const initializes arrays based on the constant-multiplicity initialization method
# described in Shima et al. (2009) https://doi.org/10.5194/acp-9-4491-2009
# n0 is the number concentration of droplets, R0 an initial radius of the droplets
# ΔV is the volume of the domain, Ns is the number of superdroplets
# the second option adds mass of solute as an attribute, M0 is the initial mass of solute
# the function draws from an exponential distribution

function init_ξ_const(Ns,ΔV,n0,R0)
    ξstart = n0*ΔV/Ns*ones(Float64,Ns)
    # R0 = Float64(30.531e-6) # meters
    X0 = Float64(4*π/3*R0^3) # initial volume m3    
    Xstart = rand(Exponential(X0), Ns)
    Rstart = (3 .*Xstart./(4*π)).^(1/3)
    return ξstart, Rstart, Xstart
end

function init_ξ_const(Ns,ΔV,n0,R0,M0)
    ξstart = n0*ΔV/Ns*ones(Float64,Ns)
    # R0 = Float64(30.531e-6) # meters
    X0 = Float64(4*π/3*R0^3) # initial volume m3
    Xstart = rand(Exponential(X0), Ns)
    Rstart = (3 .*Xstart./(4*π)).^(1/3)
    Mstart = rand(Exponential(M0), Ns)
    return ξstart, Rstart, Xstart,Mstart
end






#---------------------------------------------------------
# HYDRODYNAMIC KERNEL
#---------------------------------------------------------

# terminal velocity of droplets
function terminal_v(r) # terminal velocity 
    # Tables from
    # THE TERMINAL VELOCITY OF FALL FOR WATER DROPLETS IN STAGNANT AIR
    # Gunn, R., Kinzer, G. D. (1949), https://doi.org/10.1175/1520-0469(1949)006<0243:TTVOFF>2.0.CO;2

    if 2*r*100<0.0078
        tv=1.2*10e6*(r*100)^2
        tv = tv/100
    else    
        d_table = [0.078,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,
            2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8]./10

        v_table = [18,27,72,117,162,206,247,287,327,367,403,464,517,565,609,649,690,
        727,757,782,806,826,844,860,872,883,892,898,903,907,909,912,914,916,917]

        interpo_extrapo = linear_interpolation(d_table,v_table,extrapolation_bc=Line())
        tv = interpo_extrapo(2*r*100)/100 # radius in meters to diameter in cm, then velocity cm->m
    end
    return tv
end

#collision efficiency function
function collision_efficiency(R1,R2)
    #Parameterization from Berry 1967
    #https://doi.org/10.1175/1520-0469(1967)024<0688:CDGBC>2.0.CO;2
    r = max(R1,R2)*1e6
    rs = min(R1,R2)*1e6

    p = rs/r
    D = (-27)/(r^1.65)
    E = (-58)/(r^1.9)
    F = (15/r)^4 +1.13 
    G = (16.7/r)^8 +1 +0.004*r

    Y = 1+p+D/p^F+E/(1-p)^G
    if Y<0
        Y=0
    end

    return Y
end

###########################################################
# The hydrodynamic kernel function two methods: one for superdroplet 
# structures and one for two radius values -- it takes dummy volumes so 
# that the call could be generalized in the coalescence_timestep function

function hydrodynamic(droplet1,droplet2)
    E = collision_efficiency(droplet1.R,droplet2.R)
    Rsum = (droplet1.R+droplet2.R) # sum of radius is in meters
    vdif = abs(terminal_v(droplet1.R)-terminal_v(droplet2.R))
    return E*π*Rsum^2*vdif
end

function hydrodynamic(R1,R2,X1,X2)
    E = collision_efficiency(R1,R2)
    Rsum = (R1+R2) # sum of radius is in meters
    vdif = abs(terminal_v(R1).-terminal_v(R2))
    return E .*π .* Rsum^2 .*vdif
end

#---------------------------------------------------------
# GOLOVIN KERNEL (Golovin,1963)
#---------------------------------------------------------
###########################################################
# The golovin kernel function two methods: one for superdroplet 
# structures and one for two volume values -- it takes dummy radius values so 
# that the call could be generalized in the coalescence_timestep function

function golovin(droplet1,droplet2)
    b = 1.5e3 # seconds^-1
    Xsum = droplet1.X+droplet2.X # m3
    return b*Xsum
end

function golovin(R1,R2,X1,X2)
    b = 1.5e3 # seconds^-1
    Xsum = X1+X2 # m3
    return b*Xsum
end

#---------------------------------------------------------
# PROBABILITIES
#---------------------------------------------------------
#calc_PS calculates the probability of two superdroplets colliding in a 
#given timestep and volume,
#kernel is a kwarg that can be set to either golovin or hydrodynamic
#the function has two methods: one for superdroplet structures and one for
#free radius, volume, and multiplicity values


function calc_Ps(droplet1,droplet2,Δt,ΔV;kernel=golovin)
    Pjk = kernel(droplet1,droplet2)*Δt/ΔV
    Ps = max(droplet1.ξ,droplet2.ξ)*Pjk 
    return Ps
end

function calc_Ps(R1,R2,X,Xs,Δt,ΔV,ξj,ξk;kernel=golovin)
    # This function takes the following arguments:
    Pjk = kernel(R1,R2,X,Xs)*Δt/ΔV
    Ps = max(ξj,ξk)*Pjk
    return Ps
end

#----------------------------------------------------------
# COALESCENCE
#----------------------------------------------------------
#the coalescence_timestep function takes superdroplets and runs them
#through an All-or-nothing update timestep based on the SDM method 
#described in Shima et al. (2009)
#when the lowest multiplicity of superdroplets is less than 1, the 
#largest superdroplet is split into two equal parts, as proposed by
#Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017

#currently all are parallelized on CPUs with Threads.@threads

#there are three methods: one for superdroplet structures,
#the other two are for vectorized attribute inputs, one has radius,multiplicity, and volume
#as attributes and the other has (ξ,R,X,M) as attributes
#coalescence_timestep!(droplets,Ns,Δt,ΔV;kernel=golovin),   returns droplets
#coalescence_timestep!(ξ,R,X,Ns,Δt,ΔV;kernel=golovin),      returns ξ,R,X
#coalescence_timestep!(ξ,R,M,X,Ns,Δt,ΔV;kernel=golovin),    returns ξ,R,M,X

#coalescence_timestep_small_alpha! can handle sublinear sampling parameter y 
#coalescence_timestep_small_alpha!(droplets,Ns,y,Δt,ΔV;kernel=golovin), returns droplets


function coalescence_timestep!(droplets,Ns,Δt,ΔV;kernel=golovin)
    I =shuffle(1:Ns)
    L= [(I[l-1],I[l]) for l=2:2:length(I)]
    scale = Ns*(Ns-1)/2/(Ns/2)

    Threads.@threads for α=1:Int(floor(Ns/2))
    # for α in 1:Ns/2

        jj=L[α][1] #struggled with choosing j>k so its hard coded
        kk=L[α][2]

        if droplets[jj].ξ<droplets[kk].ξ
            k=jj
            j=kk
        else
            j=jj
            k=kk
        end			

        ϕ = rand()
        pα = scale*calc_Ps(droplets[j],droplets[k],Δt,ΔV,kernel=kernel)


        if ϕ < pα - floor(pα)
            γ = floor(pα)+1
        else ϕ >= pα - floor(pα)
            γ = floor(pα)
        end

        if γ == 0
            continue
        end

        γ_tilde = min(γ,floor(droplets[j].ξ/droplets[k].ξ))
        # if γ_tilde == floor(ξ[j]/ξ[k])
        #     println("limit exceeded")
        # end

        if droplets[j].ξ - γ_tilde*droplets[k].ξ > 0
            droplets[j].ξ = droplets[j].ξ-γ_tilde*droplets[k].ξ
            droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
            droplets[k].X = 4/3*π * droplets[k].R^3
            droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M

        elseif droplets[j].ξ - γ_tilde*droplets[k].ξ == 0
            droplets[j].ξ = floor(droplets[k].ξ/2)
            droplets[k].ξ = droplets[k].ξ - floor(droplets[k].ξ/2)
            droplets[j].R = droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
            droplets[j].X = droplets[k].X = 4/3*π * droplets[j].R^3
            droplets[j].M = droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M
        elseif droplets[j].ξ - droplets[k].ξ < 0
            print("nooooo")
        end 
    end

    #this takes a lot of searching.. could we put the failure condition somewhere else/
    # min_ξ = findmin(droplet -> droplet.ξ, droplets)
    # max_ξ = findmax(droplet -> droplet.ξ, droplets)
    # 
    # if (min_ξ[1] <= 0 && max_ξ[1] > 1)
    #     while (min_ξ[1] <= 0 && max_ξ[1] > 1)
    #         # println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)])
    #         # Split superdroplet with highest multiplicity, half goes to argmin(ξ) and half stays
    #         # Idea from Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017
    #         min_ξ = findmin(droplet -> droplet.ξ, droplets)
    #         max_ξ = findmax(droplet -> droplet.ξ, droplets)
            
    #         argmin_i = min_ξ[2]
    #         argmax_i = max_ξ[2]
    #         droplets[argmin_i].ξ = floor(max_ξ[1]/2)
    #         droplets[argmin_i].R = droplets[argmax_i].R
    #         droplets[argmin_i].X = droplets[argmax_i].X
    #         droplets[argmin_i].M = droplets[argmax_i].M
    #         droplets[argmax_i].ξ = floor(max_ξ[1]/2)
    #     end
    # elseif (min_ξ[1] <= 0 && max_ξ[1] <= 1)
    #     # Cannot split highest multiplicity superdroplet, have to remove superdroplet from system
    #     println("Superdroplet ", min_ξ[2], " has multiplicity of ", min_ξ[1], ", removing from system")
    #     deleteat!(droplets,droplets[min_ξ[2]])
    #     Ns=Ns-1 #do I need to return this??
    # end
return droplets
end



function coalescence_timestep!(ξ,R,M,X,Ns,Δt,ΔV;kernel=golovin)
    I =shuffle(1:Ns)
    L= [(I[l-1],I[l]) for l=2:2:length(I)]
    scale = Ns*(Ns-1)/2/(Ns/2)

    Threads.@threads for α=1:Int(floor(Ns/2))
    # for α in 1:Int(floor(Ns/2))

        jj=L[α][1] #struggled with choosing j>k so its hard coded
        kk=L[α][2]

        if ξ[jj]<ξ[kk]
            k=jj
            j=kk
        else
            j=jj
            k=kk
        end			

        ϕ = rand()
        pα = scale*calc_Ps(R[j],R[k],X[j],X[k],Δt,ΔV,ξ[j],ξ[k],kernel=kernel)

        if ϕ < pα - floor(pα)
            γ = floor(pα)+1
        else ϕ >= pα - floor(pα)
            γ = floor(pα)
        end

        if γ == 0
            continue
        end

        γ_tilde = min(γ,floor(ξ[j]/ξ[k]))
        # if γ_tilde == floor(ξ[j]/ξ[k])
        #     println("limit exceeded")
        # end

        if ξ[j] - γ_tilde*ξ[k] > 0
            ξ[j] = ξ[j]-γ_tilde*ξ[k]
            R[k] = (γ_tilde*R[j]^3+R[k]^3)^(1/3)
            M[k] = γ_tilde*M[j]+M[k]

        elseif ξ[j] - γ_tilde*ξ[k] == 0
            ξ[j] = floor(ξ[k]/2)
            ξ[k] = ξ[k] - floor(ξ[k]/2)
            R[j] = R[k] = (γ_tilde*R[j]^3+R[k]^3)^(1/3)
            M[j] = M[k] = γ_tilde*M[j]+M[k]
        elseif ξ[j] - ξ[k] < 0
            print("nooooo")
        end

    end
    X = 4/3 * π .* R.^3

    # #this takes a lot of searching.. could we put the failure condition somewhere else?
    # if (minimum(ξ) <= 0 && maximum(ξ) > 1)

    #     while (minimum(ξ) <= 0 && maximum(ξ) > 1)
    #         # println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)])
    #         # Split superdroplet with highest multiplicity, half goes to argmin(ξ) and half stays
    #         # Idea from Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017
    #         argmin_i = argmin(ξ)
    #         argmax_i = argmax(ξ)
    #         ξ[argmin_i] = floor(ξ[argmax_i]/2)
    #         R[argmin_i] = R[argmax_i]
    #         X[argmin_i] = X[argmax_i]
    #         M[argmin_i] = M[argmax_i]
    #         ξ[argmax_i] = floor(ξ[argmax_i]/2)
    #     end
    # elseif (minimum(ξ) <= 0 && maximum(ξ) <= 1)
    #     # Cannot split highest multiplicity superdroplet, have to remove superdroplet from system
    #     println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)], ", removing from system")
    #     deleteat!(R,argmin(ξ))
    #     deleteat!(M,argmin(ξ))
    #     deleteat!(X,argmin(ξ))
    #     deleteat!(ξ,argmin(ξ))
    #     Ns=Ns-1
    # end
return ξ,R,M,X
end


function coalescence_timestep!(ξ,R,X,Ns,Δt,ΔV;kernel=golovin)
    I =shuffle(1:Ns)
    L= [(I[l-1],I[l]) for l=2:2:length(I)]
    scale = Ns*(Ns-1)/2/(Ns/2)
    

    Threads.@threads for α=1:Int(floor(Ns/2))
    # for α in 1:Int(floor(Ns/2))


        jj=L[α][1] #struggled with choosing j>k so its hard coded
        kk=L[α][2]

        if ξ[jj]<ξ[kk]
            k=jj
            j=kk
        else
            j=jj
            k=kk
        end			

        ϕ = rand()
        pα = scale*calc_Ps(R[j],R[k],X[j],X[k],Δt,ΔV,ξ[j],ξ[k],kernel=kernel)

        if ϕ < pα - floor(pα)
            γ = floor(pα)+1
        else ϕ >= pα - floor(pα)
            γ = floor(pα)
        end

        if γ == 0
            continue
        end

        γ_tilde = min(γ,floor(ξ[j]/ξ[k]))
        # if γ_tilde == floor(ξ[j]/ξ[k])
        #     println("limit exceeded")
        # end

        if ξ[j] - γ_tilde*ξ[k] > 0
            ξ[j] = ξ[j]-γ_tilde*ξ[k]
            R[k] = (γ_tilde*R[j]^3+R[k]^3)^(1/3)

        elseif ξ[j] - γ_tilde*ξ[k] == 0
            ξ[j] = floor(ξ[k]/2)
            ξ[k] = ξ[k] - floor(ξ[k]/2)
            R[j] = R[k] = (γ_tilde*R[j]^3+R[k]^3)^(1/3)
        elseif ξ[j] - ξ[k] < 0
            print("nooooo")
        end

    end
    X = 4/3 * π .* R.^3

    # #this takes a lot of searching.. could we put the failure condition somewhere else/
    # if (minimum(ξ) <= 0 && maximum(ξ) > 1)

    #     while (minimum(ξ) <= 0 && maximum(ξ) > 1)
    #         # println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)])
    #         # Split superdroplet with highest multiplicity, half goes to argmin(ξ) and half stays
    #         # Idea from Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017
    #         argmin_i = argmin(ξ)
    #         argmax_i = argmax(ξ)
    #         ξ[argmin_i] = floor(ξ[argmax_i]/2)
    #         R[argmin_i] = R[argmax_i]
    #         X[argmin_i] = X[argmax_i]
    #         ξ[argmax_i] = floor(ξ[argmax_i]/2)
    #     end
    # elseif (minimum(ξ) <= 0 && maximum(ξ) <= 1)
    #     # Cannot split highest multiplicity superdroplet, have to remove superdroplet from system
    #     println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)], ", removing from system")
    #     deleteat!(R,argmin(ξ))
    #     deleteat!(X,argmin(ξ))
    #     deleteat!(ξ,argmin(ξ))
    #     Ns=Ns-1
    # end
return ξ,R,X
end


function coalescence_timestep_small_alpha!(droplets,Ns,y,Δt,ΔV;kernel=golovin)#,M)
    I =(sample(1:Ns, (y*2), replace = false))
    L= [(I[l-1],I[l]) for l=2:2:length(I)]
    scale = Ns*(Ns-1)/2/(y)


    Threads.@threads for α=1:y
    # for α in 1:y

        jj=L[α][1] #struggled with choosing j>k so its hard coded
        kk=L[α][2]

        if droplets[jj].ξ<droplets[kk].ξ
            k=jj
            j=kk
        else
            j=jj
            k=kk
        end			

        ϕ = rand()
        pα = scale*calc_Ps(droplets[j],droplets[k],Δt,ΔV,kernel=kernel)


        if ϕ < pα - floor(pα)
            γ = floor(pα)+1
        else ϕ >= pα - floor(pα)
            γ = floor(pα)
        end

        if γ == 0
            continue
        end

        γ_tilde = min(γ,floor(droplets[j].ξ/droplets[k].ξ))
        # if γ_tilde == floor(ξ[j]/ξ[k])
        #     println("limit exceeded")
        # end

        if droplets[j].ξ - γ_tilde*droplets[k].ξ > 0
            droplets[j].ξ = droplets[j].ξ-γ_tilde*droplets[k].ξ
            droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
            droplets[k].X = 4/3*π * droplets[k].R^3
            droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M

        elseif droplets[j].ξ - γ_tilde*droplets[k].ξ == 0
            droplets[j].ξ = floor(droplets[k].ξ/2)
            droplets[k].ξ = droplets[k].ξ - floor(droplets[k].ξ/2)
            droplets[j].R = droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
            droplets[j].X = droplets[k].X = 4/3*π * droplets[j].R^3
            droplets[j].M = droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M
        elseif droplets[j].ξ - droplets[k].ξ < 0
            print("nooooo")
        end 
    end
        #this takes a lot of searching.. could we put the failure condition somewhere else/
        # min_ξ = findmin(droplet -> droplet.ξ, droplets)
        # max_ξ = findmax(droplet -> droplet.ξ, droplets)
        
        # if (min_ξ[1] <= 0 && max_ξ[1] > 1)
        #     while (min_ξ[1] <= 0 && max_ξ[1] > 1)
        #         # println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)])
        #         # Split superdroplet with highest multiplicity, half goes to argmin(ξ) and half stays
        #         # Idea from Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017
        #         min_ξ = findmin(droplet -> droplet.ξ, droplets)
        #         max_ξ = findmax(droplet -> droplet.ξ, droplets)
                
        #         argmin_i = min_ξ[2]
        #         argmax_i = max_ξ[2]
        #         droplets[argmin_i].ξ = floor(max_ξ[1]/2)
        #         droplets[argmin_i].R = droplets[argmax_i].R
        #         droplets[argmin_i].X = droplets[argmax_i].X
        #         droplets[argmin_i].M = droplets[argmax_i].M
        #         droplets[argmax_i].ξ = floor(max_ξ[1]/2)
        #     end
        # elseif (min_ξ[1] <= 0 && max_ξ[1] <= 1)
        #     # Cannot split highest multiplicity superdroplet, have to remove superdroplet from system
        #     println("Superdroplet ", min_ξ[2], " has multiplicity of ", min_ξ[1], ", removing from system")
        #     deleteat!(droplets,droplets[min_ξ[2]])
        #     Ns=Ns-1 #do I need to return this??
        # end
return droplets
end

#---------------------------------------------------------
# all_or_nothing! is the update scheme (Shima et al, 2009) given two superdroplets
# can take two superdroplets or two sets of radius, volume, multiplicity, and ξ values
# only set to update multiplicity, radius, and volume and solute mass
# this is the parallelized part of the coalescence_timestep function above

function all_or_nothing!(R1,R2,X1,X2,ξ1,ξ2,M1,M2,scale,Δt,ΔV;kernel=golovin)

        if ξ1<ξ2
            large1 = false
            Rj = R2
            Rk = R1
            Xj = X2
            Xk = X1
            ξj = ξ2
            ξk = ξ1
            Mj = M2
            Mk = M1
        else
            large1 = true
            Rj = R1
            Rk = R2
            Xj = X1
            Xk = X2
            ξj = ξ1
            ξk = ξ2
            Mj = M1
            Mk = M2
        end			

        ϕ = rand()
        pα = scale*calc_Ps(Rj,Rk,Xj,Xk,Δt,ΔV,ξj,ξk,kernel=kernel)

        if ϕ < pα - floor(pα)
            γ = floor(pα)+1
        else ϕ >= pα - floor(pα)
            γ = floor(pα)
        end

        if γ == 0
            return R1,R2,X1,X2,ξ1,ξ2,M1,M2 #return in same order as input, no changes
        end

        γ_tilde = min(γ,floor(ξj/ξk))
        # if γ_tilde == floor(ξ[j]/ξ[k])
        #     println("limit exceeded")
        # end

        if ξj - γ_tilde*ξk > 0
            ξj = ξj-γ_tilde*ξk
            Rk = (γ_tilde*Rj^3+Rk^3)^(1/3)
            Mk = γ_tilde*Mj+Mk

        elseif ξj - γ_tilde*ξk == 0
            ξj = floor(ξk/2)
            ξk = ξk - floor(ξk/2)
            Rj = Rk = (γ_tilde*Rj^3+Rk^3)^(1/3)
            Mj = Mk = γ_tilde*Mj+Mk
        elseif ξj - ξk < 0
            print("nooooo")
        end

    Xj = 4/3 * π * Rj^3
    Xk = 4/3 * π * Rk^3

    if large1 == true
        return R1,R2,X1,X2,ξ1,ξ2,M1,M2
    else
        return R2,R1,X2,X1,ξ2,ξ1,M2,M1
    end
end

function all_or_nothing!(R1,R2,X1,X2,ξ1,ξ2,scale,Δt,ΔV;kernel=golovin)

    if ξ1<ξ2
        large1 = false
        Rj = R2
        Rk = R1
        Xj = X2
        Xk = X1
        ξj = ξ2
        ξk = ξ1

    else
        large1 = true
        Rj = R1
        Rk = R2
        Xj = X1
        Xk = X2
        ξj = ξ1
        ξk = ξ2

    end			

    ϕ = rand()
    pα = scale*calc_Ps(Rj,Rk,Xj,Xk,Δt,ΔV,ξj,ξk,kernel=kernel)

    if ϕ < pα - floor(pα)
        γ = floor(pα)+1
    else ϕ >= pα - floor(pα)
        γ = floor(pα)
    end

    if γ == 0
        return R1,R2,X1,X2,ξ1,ξ2 #return in same order as input, no changes
    end

    γ_tilde = min(γ,floor(ξj/ξk))
    # if γ_tilde == floor(ξ[j]/ξ[k])
    #     println("limit exceeded")
    # end

    if ξj - γ_tilde*ξk > 0
        ξj = ξj-γ_tilde*ξk
        Rk = (γ_tilde*Rj^3+Rk^3)^(1/3)
 
    elseif ξj - γ_tilde*ξk == 0
        ξj = floor(ξk/2)
        ξk = ξk - floor(ξk/2)
        Rj = Rk = (γ_tilde*Rj^3+Rk^3)^(1/3)
    elseif ξj - ξk < 0
        print("nooooo")
    end

    Xj = 4/3 * π * Rj^3
    Xk = 4/3 * π * Rk^3

if large1 == true
    return R1,R2,X1,X2,ξ1,ξ2
else
    return R2,R1,X2,X1,ξ2,ξ1
end
end

function all_or_nothing!(droplet1,droplet2,scale,Δt,ΔV;kernel=golovin)
    droplets = [droplet1,droplet2]
    if droplet1.ξ<droplet2.ξ    
        large1 = false
        j=2
        k=1
    else
        large1 = true
        j=1
        k=2
    end

        ϕ = rand()
        pα = scale*calc_Ps(droplets[j],droplets[k],Δt,ΔV,kernel=kernel)


        if ϕ < pα - floor(pα)
            γ = floor(pα)+1
        else ϕ >= pα - floor(pα)
            γ = floor(pα)
        end

        if γ == 0
            return droplet1,droplet2 #return in same order as input, no changes
        end

        γ_tilde = min(γ,floor(droplets[j].ξ/droplets[k].ξ))
        # if γ_tilde == floor(ξ[j]/ξ[k])
        #     println("limit exceeded")
        # end

        if droplets[j].ξ - γ_tilde*droplets[k].ξ > 0
            droplets[j].ξ = droplets[j].ξ-γ_tilde*droplets[k].ξ
            droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
            droplets[k].X = 4/3*π * droplets[k].R^3
            droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M

        elseif droplets[j].ξ - γ_tilde*droplets[k].ξ == 0
            droplets[j].ξ = floor(droplets[k].ξ/2)
            droplets[k].ξ = droplets[k].ξ - floor(droplets[k].ξ/2)
            droplets[j].R = droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
            droplets[j].X = droplets[k].X = 4/3*π * droplets[j].R^3
            droplets[j].M = droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M
        elseif droplets[j].ξ - droplets[k].ξ < 0
            print("nooooo")
        end 

    return droplet1,droplet2 #return in same order as input, changes made in place
end





###########################################################
# The coalescence_unittest_graph! function is the coalescence unit test function that runs the 
# coalescence_timestep function, plotting every 1200 seconds for 1 hour


function coalescence_unittest_graph!(droplets,Ns,Δt,ΔV;smooth=true,label=true,
    kernel=golovin,radius_bins_edges=10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128))
    
    println("Running simulation...")
    seconds = 0

    plot1 = plot()

    coal_func_time = 0
    bins = zeros(length(radius_bins_edges)-1)

    simtime = @elapsed begin
    while seconds <= 3600
   
        
        if seconds >= 0 && mod(seconds,1200) == 0
            println("Time: ", seconds, " seconds")
            # droplets = sort(droplets, by = droplet -> droplet.R)

            xx,yy = binning_dsd(droplets,radius_bins_edges,seconds,smooth=smooth,scope_init=2)
            bins = hcat(bins,yy)

                plot1 = plot!(xx,yy,lc=:black,
                    # linetype=:steppost,
                    # linewidth=2,
                    xaxis=:log,xlims=(10,10000), 
                    label = label == true ? "t = "*string(seconds)*" sec" : nothing,legend=:topright,
                    legendtitle = string(kernel)*", "*s"$N_s$=2^"*string(Int(log2(Ns))),guidefontsize=10)

                title!("Time Evolution of the Mass Density Distribution")
                xlabel!("Droplet Radius R (μm)")
                ylabel!("Mass Density Distribution g(lnR)(g/m^3/unit ln R)")

                plot!(plot1)
        end

        time = @elapsed begin
        droplets = coalescence_timestep!(droplets,Ns,Δt,ΔV,kernel=kernel)
        # println("Time: ", seconds)
        end # coalescence_timestep function timing
        coal_func_time += time


        seconds = seconds + Δt
    end
    end # simulation time
    println("simtime =",simtime)
    println("coal_func_time =",coal_func_time)

    bins = bins[:,2:end]

    return plot1
end


function coalescence_unittest_graph!(ξ,R,X,Ns,Δt,ΔV;smooth=true,label=true,kernel=golovin,radius_bins_edges=10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128))
    println("Running simulation...")
    seconds = 0

    plot1 = plot()
    coal_func_time = 0
    bins = zeros(length(radius_bins_edges)-1)

    simtime = @elapsed begin
    while seconds <= 3600
    
        if seconds >= 0 && mod(seconds,1200) == 0
            println("Time: ", seconds, " seconds")

            xx,yy = binning_dsd(X,ξ,radius_bins_edges,seconds,smooth=smooth,scope_init=2)
            bins = hcat(bins,yy)

            plot1 = plot!(xx,yy,lc=:black,
                # linetype=:steppost,
                # linewidth=2,
                xaxis=:log,xlims=(10,10000), 
                label = label == true ? "t = "*string(seconds)*" sec" : nothing,legend=:topright,
                legendtitle = string(kernel)*", "*s"$N_s$=2^"*string(Int(log2(Ns))),guidefontsize=10)
            title!("Time Evolution of the Mass Density Distribution")
            xlabel!("Droplet Radius R (μm)")
            ylabel!("Mass Density Distribution g(lnR)(g/m^3/unit ln R)")
        
        end

        time = @elapsed begin
        ξ,R,X = coalescence_timestep!(ξ,R,X,Ns,Δt,ΔV,kernel=kernel)
        end # coalescence_timestep function timing
        coal_func_time += time
        X = 4/3*π * R.^3

        seconds = seconds + Δt
    end
    end # simulation time
    println("simtime =",simtime)
    println("coal_func_time =",coal_func_time)

    return plot1
end




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
    #calc_Ps(droplet1,droplet2,Δt,ΔV;kernel=golovin),           returns Ps
    #calc_Ps(R1,R2,X,Xs,Δt,ΔV,ξj,ξk;kernel=golovin),            returns Ps
    #coalescence_timestep!(droplets,Ns,Δt,ΔV;kernel=golovin),   returns droplets
    #coalescence_timestep!(ξ,R,X,Ns,Δt,ΔV;kernel=golovin),      returns ξ,R,X
    #coalescence_timestep!(ξ,R,M,X,Ns,Δt,ΔV;kernel=golovin),    returns ξ,R,M,X


#SMALL ALPHA STUDY
    #coalescence_timestep_small_alpha!(droplets,Ns,y,Δt,ΔV;kernel=golovin), returns droplets



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
        tv=1.2*10e8*(r*100)^2
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









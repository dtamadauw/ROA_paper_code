using Statistics
using LinearAlgebra

function build_motion(mode::AbstractString, vel_sig, freq, offset, t)
    if mode == "linear"
        v = vel_sig
        r = vel_sig * t
    elseif mode == "sine"
        omega = 2*pi*freq
        v = vel_sig * cos(omega*t + offset)
        r = vel_sig/omega * sin(omega*t + offset)
    elseif mode == "random"
        r = 0.25 * vel_sig * (rand() - 0.5)
        v = vel_sig * (rand() - 0.5)
    else
        error("Unknown motion mode: $mode")
    end
    return r, v
end

function zrot(phi)
    return [cos(phi) -sin(phi) 0.0;
            sin(phi)  cos(phi) 0.0;
            0.0       0.0      1.0]
end

function xrot(phi)
    return [1.0 0.0 0.0;
            0.0 cos(phi) -sin(phi);
            0.0 sin(phi)  cos(phi)]
end

function throt(phi, theta)
    Rz = zrot(-theta)
    Rx = xrot(phi)
    return inv(Rz) * Rx * Rz
end

function freeprecess(T, T1, T2, df)
    phi = 2*pi*df*T/1000.0
    E1 = exp(-T / T1)
    E2 = exp(-T / T2)
    Afp = Diagonal([E2, E2, E1]) * zrot(phi)
    Bfp = [0.0, 0.0, 1.0 - E1]
    return Afp, Bfp
end

function RF_spoiling_2D_fast(T1, T2, TR, TE, flip, inc, Nex, Nf_in, moment, motion_in, flip_on)
    # motion_in: [pixScale, peak_velocity, motion_frequency, Gradient_Moment, Delta_t]
    motion = Dict()
    motion[:pixScale] = motion_in[1]
    motion[:fun] = (t)->build_motion("sine", motion_in[2], motion_in[3], motion_in[4], t)
    motion[:G] = motion_in[5]
    motion[:delta_t] = motion_in[6]

    phix = moment[1]
    phiy = moment[2]

    Nf_sq = round(Int, sqrt(Nf_in))
    Nf = Nf_sq^2
    mesh_temp = (collect(1:Nf_sq) ./ Nf_sq) .- 0.5
    gx_vec = Float64[]
    gy_vec = Float64[]
    for y in mesh_temp
        for x in mesh_temp
            push!(gx_vec, x)
            push!(gy_vec, y)
        end
    end
    gx = [v * phix for v in gx_vec]
    gy = [v * phiy for v in gy_vec]

    M = zeros(3, Nf)
    Msig = Vector{ComplexF64}(undef, Nex)
    Msig2 = Vector{ComplexF64}(undef, Nex)

    df = 0.0

    Ate, Bte = freeprecess(TE, T1, T2, df)
    Ate2, Bte2 = freeprecess(TR - 2*TE, T1, T2, df)
    Atr, Btr = freeprecess(TE, T1, T2, df)

    M .= 0.0
    M[3, :] .= 1.0
    on = ones(Nf)
    phase_motion = zeros(Nex)
    Rfph = zeros(Nex+1)
    Rfph[1] = 0.0
    Rfinc = inc

    fl_sign = 1
    if flip_on
        fl_sign = -1
    end

    motion_fun = motion[:fun]
    delta_t = motion[:delta_t]
    pixScale = motion[:pixScale]

    # Calculate phase shift by motion for each excitation
    for n in 1:Nex
        if flip_on
            G = motion[:G] * ((-1)^n)
        else
            G = motion[:G]
        end
        r0, v = motion_fun((n-1)*TR*1e-3)
        phase_motion[n] = pixScale * (G * r0 * delta_t + 0.5 * v * G * delta_t * delta_t)
    end

    for n in 1:Nex
        A = Ate * throt(flip * pi/180.0, Rfph[n] * pi/180.0)
        B = Bte
        # A * M + B * on
        M = A * M .+ B .* ones(1, Nf)

        Msig[n] = mean((M[1, :] .+ im .* M[2, :])) * exp(-im * Rfph[n] * pi/180.0)

        M = Ate2 * M .+ Bte2 .* ones(1, Nf)

        M = zrot(phase_motion[n]) * M

        for k in 1:Nf
            angle = gx[k] * (fl_sign^n) + gy[k] * (1^n)
            M[:, k] = zrot(angle) * M[:, k]
        end

        dess_dem = Rfph[n] + n * inc
        Msig2[n] = mean((M[1, :] .+ im .* M[2, :])) * exp(-im * dess_dem * pi/180.0)

        M = Atr * M .+ Btr .* ones(1, Nf)

        Rfph[n+1] = Rfph[n] + Rfinc
        Rfinc = (n+1) * inc
    end

    return Msig, Msig2
end

# module wrapper removed for direct inclusion

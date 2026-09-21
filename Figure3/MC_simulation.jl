using Statistics

include("RF_spoiling_2D_fast.jl")

# Quadratic vs Linear RF phase increment for spoiling (transplanted from MATLAB)
df = 0.0
T1 = 809.0 # ms
T2 = 34.0  # ms
TE = 0.0   # ms
Nf = 256

TR_DESS = 6.8 # ms
TR_PB = 6.3   # ms
flip_DESS = 35
flip_PB = 20
dphi = 2

numRO_DESS = 128
numRO_PB = 160
FOV = 32.0
Gamma = 267.52e6
Target_MonGx = 2*pi

res_DESS = FOV / numRO_DESS
Gx_DESS = Target_MonGx / (Gamma * res_DESS * TR_DESS * 1e-3)
res_PB = FOV / numRO_PB
Gx_PB = Target_MonGx / (Gamma * res_PB * TR_PB * 1e-3)

NumDummy = 100
Nex = 1024

# Set the freq range (use fixed samples for testing or generate N_freq random values in [a,b])
#Chourpiliadis, Charilaos, and Abhishek Bhardwaj. "Physiology, respiratory rate." StatPearls [Internet]. StatPearls Publishing, 2022.
a = 0.2
b = 1.2
# Configure number of frequency samples and whether to use fixed test values
N_freq = 100 # change as needed
use_fixed_freq = false
fixed_freq = [1.2]
if use_fixed_freq
    freq = fixed_freq
else
    freq = a .+ (b - a) .* rand(N_freq)
end
c = 0.0
d = 2*pi
offset = c .+ (d - c) .* rand(N_freq)


mov_int = collect(0.0:0.025:1.0)

SR = zeros(length(mov_int), length(freq))
SR_ROA = zeros(length(mov_int), length(freq))
PH = zeros(length(mov_int), length(freq))
PH_ROA = zeros(length(mov_int), length(freq))

for (ii, mi) in enumerate(mov_int)
    for (jj, fj) in enumerate(freq)
        motion_params = [numRO_DESS / FOV, mi, fj, offset[jj], Gamma * Gx_DESS, 2.56e-3]
        F0, F1 = RF_spoiling_2D_fast(T1, T2, TR_DESS, TE, flip_DESS, 0, Nex, Nf, [2*pi, 2*pi], motion_params, false)
        F0_flip, F1_flip = RF_spoiling_2D_fast(T1, T2, TR_DESS, TE, flip_DESS, 0, Nex, Nf, [2*pi, 1*pi], motion_params, true)

        val_F0 = F0[NumDummy:end]
        val_F1 = F1[NumDummy:end]
        val_F0_flip = F0_flip[NumDummy:end]
        val_F1_flip = F1_flip[NumDummy:end]

        SR[ii, jj] = abs(mean(val_F1)) / abs(mean(val_F0))
        SR_ROA[ii, jj] = abs(mean(val_F1_flip)) / abs(mean(val_F0_flip))

        motion_params = [numRO_PB / FOV, mi, fj, offset[jj], Gamma * Gx_DESS, 2.56e-3]
        S0, _ = RF_spoiling_2D_fast(T1, T2, TR_PB, TE, flip_PB, dphi, Nex, Nf, [2*pi, 2*pi], motion_params, false)
        S0_flip, _ = RF_spoiling_2D_fast(T1, T2, TR_PB, TE, flip_PB, dphi, Nex, Nf, [2*pi, 1*pi], motion_params, true)

        val_S0 = S0[NumDummy:end]
        val_S0_flip = S0_flip[NumDummy:end]

        PH[ii, jj] = angle(mean(val_S0)) + pi/2
        PH_ROA[ii, jj] = angle(mean(val_S0_flip)) + pi/2
    end
end

# Save results to MATLAB .mat file (like the original MATLAB code)
try
    using MAT
    matwrite("Simulation_Results.mat", Dict("PH"=>PH, "PH_ROA"=>PH_ROA, "SR"=>SR, "SR_ROA"=>SR_ROA))
    println("Wrote Simulation_Results.mat")
catch err
    println("Could not save .mat file: ensure MAT.jl is installed. Error: ", err)
    println("To install: import Pkg; Pkg.add(\"MAT\")")
end



function compare_and_print(name, A, B)
    d = mean(abs.(vec(A .- B)))
    if d < 1e-4
        println("Pass: $name")
    else
        println("Fail: $name (mean abs diff = $(d))")
    end
end

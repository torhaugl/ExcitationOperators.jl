using Test
using ExcitationOperators
using ExcitationOperators.BasicStuff.StandardIndices
using ExcitationOperators.BasicStuff.StandardOperators

a1 = ind(vir, "a1") #₁
i1 = ind(occ, "i1") #₁
b1 = ind(vir, "b1") #₁
j1 = ind(occ, "j1") #₁

a2 = ind(vir, "a2") #₂
i2 = ind(occ, "i2") #₂
b2 = ind(vir, "b2") #₂
j2 = ind(occ, "j2") #₂

a3 = ind(vir, "a3") #₃
i3 = ind(occ, "i3") #₃
a4 = ind(vir, "a4") #₄
i4 = ind(occ, "i4") #₄
a5 = ind(vir, "a5") #₅
i5 = ind(occ, "i5") #₅
a6 = ind(vir, "a6") #₆
i6 = ind(occ, "i6") #₆
a7 = ind(vir, "a7") #₃
i7 = ind(occ, "i7") #₃
a8 = ind(vir, "a8") #₄
i8 = ind(occ, "i8") #₄
a9 = ind(vir, "a9") #₅
i9 = ind(occ, "i9") #₅

o = ind(occ, "o")
d = ind(vir, "d")

T2 = summation(1//2 * cc_amp_tensor("t", a1, i1, b1, j1) * E(a1, i1) * E(b1, j1), [a1, i1, b1, j1])
T22 = summation(1//2 * cc_amp_tensor("t", a2, i2, b2, j2) * E(a2, i2) * E(b2, j2), [a2, i2, b2, j2])
T2u = summation(1//2 * (2//3 * cc_amp_tensor("u", a1, i1, b1, j1) + 1//3 * cc_amp_tensor("u", a1, j1, b1, i1)) * E(a1, i1) * E(b1, j1), [a1, i1, b1, j1])
S1 = summation(cc_amp_tensor("s", a3, i3) * E(a3, i3), [a3, i3])
S11 = summation(cc_amp_tensor("s", a4, i4) * E(a4, i4), [a4, i4])
S2 = summation(cc_amp_tensor("s", a5, i5, a6, i6) * E(a5, i5) * E(a6, i6), [a5, i5, a6, i6])
S22 = summation(cc_amp_tensor("s", a7, i7, a8, i8) * E(a7, i7) * E(a8, i8), [a7, i7, a8, i8])
μ = summation(real_tensor("μ", p, q) * E(p, q), [p,q])
γ = real_tensor("γ")

hF = summation( (real_tensor("F", p, q) + summation(-2//1 * psym_tensor("g", p,q,o,o) + psym_tensor("g", p,o,o,q), [o])) * E(p, q), [p,q])
H0 = h + g + μ * (S1 + γ + S2)
HF = hF + g + μ * (S1 + γ + S2)

# PROJECTION MANIFOLD
P1 = 1//2 * E(i, a)
P2 = (1//3 * E(i, a) * E(j, b) + 1//6 * E(i, b) * E(j, a)) 

# ENERGY
ECCSD = exval( 1//2 * H0 + 1//2 * HF + comm(HF, T2u)) |> simplify
@show ECCSD;

# SINGLES
Ω1 = exval( P1 * (HF + comm(HF, T2u))) |> simplify
@show Ω1;

# DOUBLES
Ω2 = exval( P2 * (HF + comm(HF, T2) + 1//2 * comm(comm(HF, T2), T22))) |> simplify
@show Ω2;

# PHOTONS
HF1 = comm(HF, S11) + comm(HF, S22) + real_tensor("ω") * (S1 + γ + S2) + μ
Ω0P = exval( HF1 + comm(HF1, T2u)) |> simplify
@show Ω0P;

# PHOTONS SINGLES
Ω1P = exval( P1 * (HF1 + comm(HF1, T2u))) |> simplify
@show Ω1P;

# PHOTONS DOUBLES
Ω2P = exval( P2 * (HF1 + comm(HF1, T2))) |> simplify
@show Ω2P;

# Introduce Paibj back into the equations
#PΩ2 = collapse_perm(Ω2, [(a,i),(b,j)])
#@show PΩ2;

#PΩ2P = collapse_perm(Ω2P, [(a,i),(b,j)])
#@show PΩ2P;
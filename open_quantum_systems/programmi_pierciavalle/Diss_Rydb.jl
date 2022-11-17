using QuantumOptics
using PyPlot


pi=3.1415926535897

k = 8 # Decay rate (in MHz) 
Ω = 0.1 # Rabi freq. = Transverse field (MHz)
Δ = 0 # detuning = longitudinal field (MHz)
γ =2*pi*0.7  # Rydberg-ground dephasing (MHz)
Γecc = 10 #Excitation rate (in MHz)
omega0=0
omega1=749000000
omega2=1048000000


#for Ω=0:100
#Ω=Ω*0.1

#for k=1:20


T = [0:0.005:5;] # Array of times (from 0 to 10 with 0.1-length steps)

ba=NLevelBasis(3)

psi1=nlevelstate(ba,1)
psi2=nlevelstate(ba,2)
psi3=nlevelstate(ba,3)

proj1=projector(psi1)
proj2=projector(psi2)
proj3=projector(psi3)

sigma13=transition(ba,1,3)
sigma31=transition(ba,3,1)
sigma23=transition(ba,2,3)
sigma12=transition(ba,1,2)
sigma32=transition(ba,3,2)
sigma13=transition(ba,1,3)

sigmaz=proj1-proj3


H = omega0*proj1 + omega1*proj2 + (omega1-Δ)*proj3 + (Ω/2)*sigma23 + (Ω/2)*sigma32

J = [sqrt(k)*sigma12, sqrt(Γecc)*sigma31, sqrt(γ)*sigmaz]

psi0=psi1

# Master Equation evolution

tout, ρt_master = timeevolution.master(T, psi0, H, J) # Master eq. evolution


excitations = expect(proj1, ρt_master) 
exc = real(excitations)


m = open("N1_delta0_Om40_k8_psi0_1_EXC.dat", "w")


for t=0:1000
 out=0
 out=exc[t+1]
 t=t*0.005
 write(m, "$t   $out",  "\n")
 #println(t, "   ", exc[t+1])
end

close(m)

#tout, ρ_master = steadystate.master(H, J)
#println(k, " ",real(expect(proj3, ρ_master[end])))

#plot(T, exc, label=L"$Δ=-1, Ω = 0.5$")
#plot(T, exc_r-exc_l, label="n_{5}, Ω = 1")
#xlim(0, 2000)
#ylim(-1, +1)
#xlabel(L"\mathrm{Time}")
#ylabel(L"$N_{exc}$")
#title(L"C_{6}=1, N=8, k=α=β=0.5, Rydberg")
#legend()

#savefig("./Excis.pdf", dpi = 100, format = "pdf", transparent = false)
#end

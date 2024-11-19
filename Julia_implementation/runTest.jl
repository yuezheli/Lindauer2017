using DifferentialEquations
using Plots
using CSV
using DataFrames

# The ode file containing the system of equations in the model
# note: this file needs to be included only ONCE in each Julia session, 
include("Lindauer_mus.jl")

# define simulation time span
tspan = (0.0, 480); # hr

# define dosing events
# use dose = 10mg/kg, dose at day 0, day7, and day 24
MW = 149000; # antibody molecular weight
mouseweight = 20; # mouse weight = 20g
dose = 10 * mouseweight *1000/(MW); 
dosetimes = [7*24, 14*24]; 
condition(u,t,integrator) = t ∈ dosetimes
affect!(integrator) = integrator.u[1] += dose
cb = DiscreteCallback(condition,affect!)

# define initial condition
u0 = [dose,0,0,0,0, 49800, 0,0,0,7e-6, 170]; 

prob = ODEProblem(ode_Lindauer!, u0, tspan, callback=cb, tstops=dosetimes);
num_sol = solve(prob, Tsit5(), reltol=1e-8);

V1 = 1.26/1000;


plasmamAb = [(u[1]/V1)*MW/1000000 * (1-0.197)  - 0.0658  for u in num_sol.u]
tumorsize = [u[11] for u in num_sol.u]

# read in data from mrgsolve
df = CSV.read("../data/mrgresult.csv", DataFrame);


p1 = plot(num_sol.t/24, plasmamAb,
            linewidth=1, xlabel = "time (days)", legend=:bottomleft,
            ylabel="plasma mAb concentation (mg/L)", label = "Julia solution")
plot!(df[:, 1]/24, df[:, 2], label = "mrgsolve solution", color = 2, seriestype = :scatter)


p2 = plot(num_sol.t/24, tumorsize,
            linewidth=1, xlabel = "time (days)",  
            ylabel="tumor size (mm^3)",  label = "Julia solution")
plot!(df[:, 1]/24, df[:, 3], label = "mrgsolve solution", color = 2, seriestype = :scatter)


plotd = plot(p1, p2, layout = @layout [a b])
savefig(plotd,"../img/julia_CP_TV.png")


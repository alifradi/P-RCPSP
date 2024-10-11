include("libraries.jl")
include("data.jl")
# PICK SOLVER
m = Model(CPLEX.Optimizer)
# Import constraints and variables
include("Model.jl")

# SOLVE

optimize!(m)

# get values

solution_summary(m, verbose=true)


#objective function values

objective_value(m)



#EXPORT DATA 

include("plotGantt.jl")
CSV.write("D:\\00- Maitrise ALI\\Projet 2\\Implémentation\\schedule.csv", df_gantt)
CSV.write("D:\\00- Maitrise ALI\\Projet 2\\Implémentation\\ressources.csv", df_ressource)
CSV.write("D:\\00- Maitrise ALI\\Projet 2\\Implémentation\\ressources_job.csv", df_ressource_job)
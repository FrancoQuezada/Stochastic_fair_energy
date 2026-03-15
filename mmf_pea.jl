include("structuresMulti.jl")

function solve_lex_pea_step(inst::InstanceM, iter::Int64, optim::Vector{Float64})
    scenarioTree=inst.tree
    model = Model(CPLEX.Optimizer)

    @variable(model, p[1:inst.J,1:scenarioTree.V] >= 0)
    @variable(model, zeta[1:iter])
    @variable(model, daux[1:iter,1:inst.J] >= 0)

    @constraint(model, [n in 1:scenarioTree.V], sum(p[j,n] for j in 1:inst.J) == inst.c_pv[n])
    @constraint(model, [j in 1:inst.J, s in eachindex(scenarioTree.scenarios)],
        sum(p[j,n] for n in scenarioTree.scenarios[s]) <= sum(inst.d[j,n] for n in scenarioTree.scenarios[s]))

    scenario_prob=[scenarioTree.rho[scenario[inst.T]] for scenario in scenarioTree.scenarios]
    @expression(model, alloc[j in 1:inst.J],
        sum(scenario_prob[s]*sum(p[j,n] for n in scenarioTree.scenarios[s]) for s in eachindex(scenarioTree.scenarios)))

    @constraint(model, [n in 1:iter, j in 1:inst.J], zeta[n] - daux[n,j] <= alloc[j])
    lex_eps_abs = 1e-3
    @constraint(model, [n in 1:(iter-1)], n*zeta[n] - sum(daux[n,:]) >= optim[n]- lex_eps_abs)
    @objective(model, Max, iter*zeta[iter] - sum(daux[iter,:]))

    set_silent(model)
    optimize!(model)
    step_time = round(solve_time(model), digits=2)

    if has_values(model)
        return objective_value(model), value.(p), step_time
    end
    return 0.0, zeros(inst.J, scenarioTree.V), step_time
end

function lexico_mmf_pea_targets(inst::InstanceM)
    ω=zeros(inst.J)
    p_last=zeros(inst.J, inst.tree.V)
    total_solve_time=0.0
    for k in 1:inst.J
        obj, p_sol, step_time=solve_lex_pea_step(inst, k, ω)
        total_solve_time+=step_time
        ω[k]=obj
        p_last.=p_sol
    end

    scenarioTree=inst.tree
    scenario_prob=[scenarioTree.rho[scenario[inst.T]] for scenario in scenarioTree.scenarios]
    scenario_target=[sum(p_last[j,n] for n in scenario) for j in 1:inst.J, scenario in scenarioTree.scenarios]
    expected_target=[sum(scenario_prob[s]*scenario_target[j,s] for s in eachindex(scenarioTree.scenarios)) for j in 1:inst.J]
    return expected_target, scenario_target, round(total_solve_time, digits=2)
end

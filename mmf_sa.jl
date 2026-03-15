include("structuresMulti.jl")
include("parametersMS.jl")


function solve_lex_step(inst::InstanceM, iter::Int64, optim::Vector{Float64},sol::SolutionM)
    """
    function to solve model L_k in paper
    """
    scenarioTree=inst.tree
    model = Model(CPLEX.Optimizer)
    timePeriods=createTime(scenarioTree)
    ################Variables#######################
    @variable(model,s[1:scenarioTree.V]>=0) #total battery
    @variable(model, I[1:inst.J,1:scenarioTree.V]>=0) #grid
    @variable(model, G[1:inst.J,1:scenarioTree.V]>=0) #grid
    @variable(model, z[1:inst.J,1:scenarioTree.V]>=0) #charge of battery
    @variable(model, y[1:inst.J,1:scenarioTree.V]>=0) #discharge of battery
    @variable(model, p[1:inst.J,1:scenarioTree.V]>=0) #photovoltaic production
    @variable(model, lambda[1:inst.J,1:scenarioTree.V]>=0)
    ##################dual variables##################
    @variable(model,zeta[1:iter]) 
    @variable(model,d[1:iter,1:inst.J]>=0)
    ################Constraints###################
    @constraint(model, [j in 1:inst.J, n in 1:scenarioTree.V], p[j,n] ==lambda[j,n]*inst.c_pv[n])
    @constraint(model, [t in 1:scenarioTree.V], sum(lambda[:,t]) == 1)
    # @constraint(model, [n in 1:scenarioTree.V], sum(p[j,n] for j in 1:inst.J)==inst.c_pv[n])
    @constraint(model, [n in 1:scenarioTree.V], s[n]<=inst.s_max)
    @constraint(model, [n in 2:scenarioTree.V], s[n]==s[scenarioTree.parents[n]]+inst.delta*inst.e_c*sum(z[j,n] for j in 1:inst.J)-inst.delta*sum(y[j,n] for j in 1:inst.J)/inst.e_d)
    @constraint(model, s[1]==inst.s_I+inst.delta*inst.e_c*sum(z[j,1] for j in 1:inst.J)-inst.delta*sum(y[j,1] for j in 1:inst.J)/inst.e_d)
    @constraint(model, [n in 1:scenarioTree.V], sum(y[j,n] for j in 1:inst.J)<=inst.f_bar)
    @constraint(model, [n in 1:scenarioTree.V], sum(z[j,n] for j in 1:inst.J)<=inst.f_under)
    @constraint(model, [j in 1:inst.J, n in 1:scenarioTree.V], inst.d[j,n]==p[j,n]+y[j,n]+I[j,n]-z[j,n]-G[j,n])
    @constraint(model, [n in 1:scenarioTree.V ; timePeriods[n]==inst.T],s[n]==inst.s_I)
    @constraint(model,[n in 1:scenarioTree.V], s[n]>=inst.s_min)
    ########################costs###########################
    @expression(model, costs[j in 1:inst.J], (inst.delta)*sum(scenarioTree.rho[n]*(inst.mu*y[j,n] +inst.nu[j,timePeriods[n]]*I[j,n] - inst.beta*G[j,n])  for n in 1:scenarioTree.V)) 
    @constraint(model, [j in 1:inst.J], costs[j]>=0)
    ########################Lexico Constraints##############
    @constraint(model,[n in 1:iter, j in 1:inst.J],zeta[n]-d[n,j]<=(inst.delta)*sum(scenarioTree.rho[t]*inst.nu[j,timePeriods[t]]*inst.d[j,t] for t in 1:scenarioTree.V)-costs[j])
    # Numerical safeguard: keep previous lex levels with a tiny tolerance
    # to avoid infeasibilities caused by floating-point noise in optim[n].
    lex_eps_abs = 1e-3
    @constraint(model,[n in 1:(iter-1)],
        n*zeta[n]-sum(d[n,:]) >= optim[n] - (lex_eps_abs))
    ###########################Objective function###########
    @objective(model, Max, iter*zeta[iter]-sum(d[iter,:]))


    # set_attribute(model, "CPXPARAM_TimeLimit", 60*10)
    # set_attribute(model,"CPXPARAM_Simplex_Tolerances_Feasibility",1e-4)
    # set_attribute(model, "CPXPARAM_Preprocessing_Presolve", 0)
    # set_attribute(model, "CPXPARAM_Emphasis_Numerical", 1)

    
    
    set_silent(model)
    optimize!(model)
    term=termination_status(model)
    pstat=primal_status(model)

    if has_values(model)
        println("Lex step ", iter, " objective: ", round(objective_value(model),digits=3), " (", term, ")")
        sol.costs=value.(costs)
        sol.y=value.(y)
        sol.z=value.(z)
        sol.I=value.(I)
        sol.p=value.(p)
        sol.s=value.(s)
        sol.G=value.(G)
        sol.w=zeros(inst.J,inst.T)
        sol.x=zeros(inst.J,inst.T)
        step_time=round(solve_time(model), digits=2)
        sol.time=step_time
        sol.status=true
        sol.id=inst.id
        # zio=[sol.I[j,t]*sol.G[j,t]==0 for j in 1:inst.J, t in 1:inst.T]
        # println(all(zio))
        # sol=Solution(sTotAux,iAux,GAux,xAux,wAux,zAux,yAux,pAux,costsAux,true,solTime,inst.id)
        return objective_value(model), step_time
    else
        step_time=round(solve_time(model), digits=2)
        println("Lex step ", iter, " failed: term=", term, ", primal=", pstat, ", time=", step_time)
        sol.costs=fill(Inf,inst.J)
        sol.y=zeros(inst.J,inst.T)
        sol.z=zeros(inst.J,inst.T)
        sol.I=zeros(inst.J,inst.T)
        sol.p=zeros(inst.J,inst.T)
        sol.s=zeros(inst.T)
        sol.G=zeros(inst.J,inst.T)
        sol.w=zeros(inst.J,inst.T)
        sol.x=zeros(inst.J,inst.T)
        sol.time=step_time
        sol.status=false
        sol.id=inst.id
        # sol=Solution(sTotAux,iAux,GAux,xAux,wAux,zAux,yAux,pAux,costsAux,false,solTime,inst.id)   
        return 0, step_time
    end
end 

function lexico(inst::InstanceM)
    """
    Algorithm lexicographic in paper
    """
    println("##########################################################")
    println("Lexicographic Algorithm")
    ω=zeros(inst.J) #initialize omega
    x_sol=SolutionM() #initial empty solution
    total_solve_time=0.0
    for k in 1:inst.J
        obj, step_time=solve_lex_step(inst,k,ω,x_sol)
        total_solve_time+=step_time
        if !x_sol.status
            println("Lexicographic algorithm stopped at step ", k, " (no solution values).")
            break
        end
        ω[k]=deepcopy(obj)
    end
    x_sol.time=round(total_solve_time, digits=2)
    println("end")
    println("##########################################################")
    # println([(inst.delta)*sum(inst.nu[j,t]*inst.d[j,t] for t in 1:inst.T) for j in 1:inst.J]) #for debug
    return x_sol    
end

lexico_mmf_sa(inst::InstanceM)=lexico(inst)

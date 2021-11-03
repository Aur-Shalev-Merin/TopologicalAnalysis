using CSV,DataFrames
using LightGraphsFlows
import LightGraphs
const lg = LightGraphs
import StatsBase
using SparseArrays
#using Convex, Mosek, MosekTools
using Base.Threads
import Clp  # use your favorite LP solver here
#using MathProgBase,Gurobi


function ret_weights(dict_w,N,W_code_to_idx,vmap)
    # Find the distribution of networks in the connected comp of the flip graph
    weight_arr = zeros(Float64,N)
    key_arr = collect(keys(dict_w))
    for k in key_arr
        if haskey(W_code_to_idx,k)
            weight_arr[W_code_to_idx[k]] = dict_w[k]
        end
    end
    weight_arr = weight_arr[vmap]
    return weight_arr./sum(weight_arr)
end

function load_w_graph(network_save_file)
    # Wrapper to load the flip graph, and to reduce to the largest connected comp
    g = lg.loadgraph(network_save_file*".lgz")
    #gplot(g)
    dat_in = CSV.read(network_save_file*".txt",DataFrame)
    W_code_to_idx = Dict(dat_in.codes .=> dat_in.index)
    W_idx_to_code = Dict(dat_in.index .=> dat_in.codes)

    ccomp = lg.connected_components(g)
    idx = argmax([length(c) for c in ccomp])
    g_con, vmap = lg.induced_subgraph(g, ccomp[idx])

    return g_con,vmap,lg.nv(g),W_code_to_idx,W_idx_to_code
end

function W_dist(g_undirected,S,ret_flow = false;tol=5e-5)

    # First create a directed version of g_undirected with edges in both
    # directions and capacities/costs equal to 1
    Nv = lg.nv(g_undirected)
    g = lg.DiGraph(Nv+2) # Create a flow-graph
    w = spzeros(Nv+2,Nv+2)
    capacity = spzeros(Nv+2,Nv+2)

    for e in lg.edges(g_undirected)
        lg.add_edge!(g,lg.src(e),lg.dst(e))
        lg.add_edge!(g,lg.dst(e),lg.src(e))
        capacity[lg.src(e),lg.dst(e)] = 1.
        capacity[lg.dst(e),lg.src(e)] = 1.
        w[lg.src(e),lg.dst(e)] = 1.
        w[lg.dst(e),lg.src(e)] = 1.
    end

    # Take the Nv+1th vertex as a source and Nv+2th as a sink.
    demand = spzeros(Nv+2,Nv+2)
    for i in 1:Nv
        if S[i] > 0
            lg.add_edge!(g,Nv+1,i)
            capacity[Nv+1,i] = S[i]*(1+tol)
        elseif S[i] < - 0.1*tol
            lg.add_edge!(g,i,Nv+2)
            demand[i,Nv+2] = -S[i]*(1-tol)
            capacity[i,Nv+2] = 1.

        end
    end
    # call min cost flow
    flow = mincost_flow(g,spzeros(lg.nv(g)), capacity , w, Clp.Optimizer,
                edge_demand=demand, source_nodes=[Nv+1], sink_nodes=[Nv+2])
    if ret_flow
        return flow
    else
        return sum(abs.(flow[1:Nv,1:Nv]))
    end
end

function sample_dist(probs,N)
    w = StatsBase.Weights(probs)
    p_emp = zeros(length(probs))
    for i = 1:N
        p_emp[StatsBase.sample( w)] += 1
    end
    return p_emp/sum(p_emp)
end

function geo_core(w_in,flow,α,N)
    weight_geo = deepcopy(w_in)

    I,J,V = findnz(flow)
    for i = 1:length(I)
        if (I[i] > length(weight_geo)) || (J[i] > length(weight_geo))
            continue
        end
        weight_geo[I[i]] -= α*V[i]
        weight_geo[J[i]] += α*V[i]
    end
    return sample_dist(weight_geo,N)
end



function geodesic(w1,w2,α,network_save_file,N)
    w_vec_in = [w1;w2]
    g,vmap,N,W_code_to_idx,W_idx_to_code = load_w_graph(network_save_file)
    weight = [ret_weights(w_vec_in[i],N,W_code_to_idx,vmap) for i in 1:length(w_vec_in)]
    flow = W_dist(g,weight[1] .- weight[2],true)
    w_geo = [geo_core(weight[1],flow,α0,N) for α0 in α]
    return w_geo
end

function find_reg_geo(p0,p1,g,k)
        flow = W_dist(g,p0 .- p1,true;tol=0.0)
        flow[findall(flow .< 0)] .= 0
        flow = flow[1:length(p0),1:length(p0)]
        I,J,V = findnz(flow)
        W1_cost = sum(abs.(flow))
        verts = unique(vcat(I,J))
        verts_inv = Dict(verts .=> 1:length(verts))
        num_v = length(verts)
        num_e = length(I)
        p0_v = p0[verts]
        p1_v = p1[verts]

        q = Variable(num_v,k+1)
        F = Variable(num_e,k)

        obj_ = sum([quadoverlin(F[e,i], q[verts_inv[I[e]],i]) +
                quadoverlin(F[e,i], q[verts_inv[J[e]],i+1])  for i = 1:k for e = 1:num_e])

        con1 = q >= 0
        con2 = [-sum(F[findall(I .== j),i]) + sum(F[findall(J .== j),i]) == q[verts_inv[j],i+1]-q[verts_inv[j],i] for i = 1:k for j in verts]
        con3 = q[:,1] == p0_v
        con4 = q[:,end] == p1_v
        con5 = F >= 0.0
        con6 = sum(F) <= W1_cost*1.08
        problem = minimize(obj_, [con1;con2;con3;con4;con5;con6])
        solve!(problem, () -> Mosek.Optimizer())
        println(problem.status)

        q_full = zeros(length(p0),k+1)
        for i = 1:(k+1)
                q_full[verts,i] .= q.value[:,i]
        end
        return q_full, k*problem.optval
end

function geodesic_reg(network_save_file,w1,w2,k;α=1.0,β=0.0)
        g,vmap,N,W_code_to_idx,W_idx_to_code = load_w_graph(network_save_file)
        w_vec_in = [w1;w2]
        weight = [ret_weights(w_vec_in[i],N,W_code_to_idx,vmap) for i in 1:length(w_vec_in)]
        q,val = find_reg_geo(α*weight[1]+(1.0-α)*weight[2],β*weight[1]
                                +(1.0-β)*weight[2],g,k)
        return q,val
end

function distance_OT(g,W)

    try
    	d =  W_dist(g,W)
        return d
    catch
    	d = W_dist(g,-W)
        return d
    end
end

function rem_self_edges!(g)
    # Removes self loops
    for e in collect(lg.edges(g))
        if lg.src(e) == lg.dst(e)
            lg.rem_edge!(g,e)
        end
    end
end

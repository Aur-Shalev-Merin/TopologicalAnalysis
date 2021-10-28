using Statistics: mean
#using PyCall
using LinearAlgebra: norm
#scipy = pyimport("scipy.spatial.qhull")



function shape_find(index_points,vert,regions,ridge_pts,ridge_vert,ind)
    # This function takes in an index and returns the voronoi cell associated
    # with that index, and a graph of the vertices and edges of the polytope
    kk = index_points[ind]
    if any(regions[kk] .== -1)
        println("error in shape_find")
        return -1
    end

    poly_pts = vert[regions[kk].+1,:] # the actual points of the polytope
    faces = ridge_vert[ (ridge_pts[:,1] .== ind) .| (ridge_pts[:,2] .== ind)]
    point_lookup = Dict(regions[kk].+1 .=> 1:length(regions[kk]) )

    g = SimpleGraph(length(regions[kk]))

    for s in 1:length(faces)
        face = faces[s].+1
        for i = 1:(length(face)-1)
            add_edge!(g, point_lookup[face[i]], point_lookup[face[i+1]])
        end
        add_edge!(g,point_lookup[face[end]],point_lookup[face[1]])
    end
    face = faces[1].+1
    face_points = [point_lookup[face[i]] for i in 1:length(face)]

    return g, poly_pts, face_points
end



function right_label(Points,g,ind)
    # This function orders the neighbouring points around a vertex, ind, according to their
    # counter clockwise orientation, where the the out of plane vector is defined
    # to be from the origin to the vertex. Shift origin so it lies within the shape
    Points = Points - repeat(mean(Points,dims=1),size(Points,1),1)
    renorm = x-> x/norm(x)
    e1 = renorm(Points[ind,:])
    nbs = neighbors(g,ind)
    a = Points[nbs[1],:]
    e2 = renorm(a - dot(a,e1)*e1) #/norm(a - dot(a,e1)*e1)
    e3 = cross(e1,e2)
    plane_angle = [angle(dot(Points[n,:],e2) + dot(Points[n,:],e3)*im) for n in nbs[2:end]]
    plane_angle = [(p < 0) ? p + 2π : p for p in plane_angle]
    order =  nbs[ vcat([1], (sortperm(plane_angle).+1))]
    return order
end

function weinberg_vect(g_aug,e,order_mat; mirror=false)
    # This function is the core part of the Weinberg algorithm. It takes the
    # augmented graph (with forward and backward edges) and a starter edge and
    # finds the Weinberg code. The mirror option either picks the right most (default) or
    # left most available edge (or right most in the mirrored representation)
    vertex_to_code = Dict(src(e) => 1 ) # convert vertex numbers to their code numbering
    code = [1]
    current_vertex = dst(e)
    previous_vertex = src(e)
    while ne(g_aug) > 0

        if haskey(vertex_to_code, current_vertex)
            push!(code,vertex_to_code[current_vertex])
            if has_edge(g_aug,current_vertex,previous_vertex)
                rem_edge!(g_aug, previous_vertex,current_vertex)
                previous_vertex,current_vertex = current_vertex,previous_vertex
                continue
            end
        else
            vertex_to_code[current_vertex] = maximum(code)+1
            push!(code,vertex_to_code[current_vertex])
        end


        edge_choices = order_mat[current_vertex]

        ind = findfirst(x->x==previous_vertex,edge_choices)
        if mirror
            edge_choices = circshift(edge_choices, 1-ind)
            next_vertex_ind = findlast(ed -> has_edge(g_aug, current_vertex,ed),edge_choices)
        else
            edge_choices = circshift(edge_choices, -ind)
            next_vertex_ind = findfirst(ed -> has_edge(g_aug, current_vertex,ed),edge_choices)
        end
         if isnothing(next_vertex_ind)
             return code
         end

        rem_edge!(g_aug, previous_vertex,current_vertex)
        previous_vertex,current_vertex = current_vertex,edge_choices[next_vertex_ind]
    end
    println("Error occured")
end

function augmented_graph(g)
    g_aug = SimpleDiGraph(size(g)[1])
    for e in edges(g)
        add_edge!(g_aug, src(e), dst(e))
        add_edge!(g_aug, dst(e), src(e))
    end
    return g_aug
end

#=
function weinberg_find(g,order_mat)
    g_aug = augmented_graph(g)

    vecs = Array{Int64}[]
    for e in edges(g_aug)
        append!(vecs,[weinberg_vect(augmented_graph(g),e,order_mat;mirror=true)])
        append!(vecs,[weinberg_vect(augmented_graph(g),e,order_mat;mirror=false)])

    end
    sort!(vecs)
    S = Int64(length(vecs)/length(unique(vecs))) # size of the symmetry group
    return vecs[1], S
end
=#
function weinberg_find!(code_tot,S_tot,kk,g,order_mat,cent_node = -1)
    g_aug = augmented_graph(g)
    vecs = Array{Int64}[]
    for e in edges(g_aug)
        if (cent_node < 0) || (src(e) == cent_node)
            append!(vecs,[weinberg_vect(augmented_graph(g),e,order_mat;mirror=true)])
            append!(vecs,[weinberg_vect(augmented_graph(g),e,order_mat;mirror=false)])
        end
    end
    sort!(vecs)
    try
        S_tot[kk] = Int64(length(vecs)/length(unique(vecs))) # size of the symmetry group
        code_tot[kk] =  copy(vecs[1])

    catch
        savegraph("debug.lgz", g)
        S_tot[kk] =  -1

    end

end

#=
function weinberg_find_total(index_points,vert,regions,ridge_pts,ridge_vert)

    code_tot = []
    S_tot = []

    for k in 1:length(index_points)
        if any(regions[index_points[k]] .== -1)
            append!(code_tot,[0])
            append!(S_tot,0)
            continue
        end

        g, poly_pts,side = shape_find(index_points,vert,regions,ridge_pts,ridge_vert,k)
        order_mat = [right_label(poly_pts,g,n) for n in 1:size(poly_pts,1)]
        code_unique,S = weinberg_find(g,order_mat)
        append!(code_tot,[code_unique])
        append!(S_tot,S)
    end

    return code_tot, S_tot
end
=#

function edge_neighbors(g,edge_index)
    e_neighbors = Array{Int64}(undef,0)
    for e in edge_index
        append!(e_neighbors,neighbors(g,e))
    end
    return unique!(e_neighbors)
end


function weinberg2D_core(g,order_mat,N,r=2)
    code_tot = Vector{Array{Int64}}(undef,N)
    S_tot = Array{Int64}(undef,N)

    for k in 1:N
        nbh  = neighborhood(g,k,r)
        g_ego, vmap = induced_subgraph(g,nbh)
        vmap_inv = Dict(vmap[k] => k for k in 1:length(vmap))
        order_local = order_mat[vmap]

        for i = 1:length(vmap)
            total_order = map.(x -> get(vmap_inv, x, -1), order_mat[vmap[i]])
            order_local[i] = total_order[total_order .> 0]
        end
        try
            weinberg_find!(code_tot,S_tot,k,g_ego,order_local,vmap_inv[k])
        catch
            println(k)
            #println(order_mat[vmap])
            #println([e for e in edges(g_ego)])
            i = 1
            println(vmap[i])
            println(order_mat[vmap[i]])
            println(vmap_inv)
            #total_order = map.(x -> get(vmap_inv, x, -1), order_mat[vmap[i]])
            #println(total_order)

            error()
        end
    end
    return code_tot,S_tot
end

function motif_size_find(Pos,r=2;periodic=true)
    Positions = deepcopy(Pos)
    N = length(Positions)
    index_ref = [x for x in 1:N]
    if periodic; periodic_extend!(Positions,index_ref); end
    # All other operations are performed both periodic and non-periodic systems,
    # but for non-periodic systems should have no effect
    p, simplices, neighbrs, edge_index = Delaunay_find(Positions)
    g_full = graph_construct(simplices,length(Positions))
    motif_lens = [length(neighborhood(g_full,k,r)) for k in 1:N]
    return motif_lens
end

function weinberg2D_old(Positions,periodic=false,r=2; α = 0)
    N = length(Positions)
    index_ref = [x for x in 1:N]
    if periodic; periodic_extend!(Positions,index_ref); end
    # All other operations are performed both periodic and non-periodic systems,
    # but for non-periodic systems should have no effect
    p, simplices, neighbrs, edge_index = Delaunay_find(Positions, α = α)
    g_full = graph_construct(simplices,length(Positions))


    order_mat = Vector{Array{Int64}}(undef,nv(g_full))
    for kk = 1:nv(g_full)
        N_list = neighbors(g_full,kk)
        theta = [atan(p[s,2]-p[kk,2], p[s,1]-p[kk,1]) for s in N_list]
        order_mat[kk] = N_list[sortperm(theta)]
    end

    order_mat = order_mat[1:N]
    for i = 1:N
        order_mat[i] = map(x -> index_ref[x],order_mat[i])
    end

    rows_keep = [minimum(simplices[k,:]) <= N for k in 1:size(simplices,1)]
    simplices_period = map.(x -> index_ref[x],simplices[rows_keep,:])
    g_period = graph_construct(simplices_period,N)

    println(N)
    println(nv(g_period))

    if periodic
        return weinberg2D_core(g_period,order_mat,r)
    else
        e_2 = edge_neighbors(g_period,edge_index)
        idx = setdiff(1:nv(g_period), e_2)

        code_tot,S_tot = weinberg2D_core(g_period,order_mat,r)
        return code_tot[idx], S_tot[idx], idx
    end
end



function weinberg2D(path_to_dir_in,params_in,r)

    g = loadgraph(path_to_dir_in*"_graph.lgz")
    order_mat = map.(s -> parse(Int64, s), split.(eachline(path_to_dir_in*"_order_mat.txt"),'\t'))
    edge_index = readdlm(path_to_dir_in*"_edge_nodes.txt", '\t', Int, '\n')
    periodic = (params_in["Periodic"] == "True")

    if periodic
        N = Int(params_in["Original vertex number"])
    else
        N = nv(g)
    end

    code_tot,S_tot = weinberg2D_core(g,order_mat,N,r)
    idx = 1:N

    if !periodic
        for i = 1:r-1
            edge_index = edge_neighbors(g,edge_index)
        end
        idx = setdiff(1:nv(g), edge_index)
    end

    return idx, code_tot[idx], S_tot[idx]
end

using Combinatorics
using JuMP
using Gurobi


"""
solve(file)

# Arguments
- `file`: name of the file to read

# Example
solve("puzzles/sheet.txt")
"""
function solve(file::String)
    islands, intersections = readFile(file)
    n = length(islands)

    model = Model(Gurobi.Optimizer)
    @variable(model, x[1:n, 1:n, 1:2], Bin)  # x[i, j, l] iff at least l bridges connected island i and island j
    @variable(model, y[1:n, 1:n, 1:n], Bin)  # y[i, j, l] iff edge (i, j) can be reached from source edge in at most l steps

    # by defintion x[i, j, 2] implies x[i, j, 1]
    @constraint(model, [i=1:n, j=1:n], x[i,j,2] <= x[i,j,1])

    # by definition x must be syymetric
    @constraint(model, [i=1:n, j=1:n, l=1:2], x[i,j,l] == x[j,i,l])

    # 1) Correct bridge count
    @constraint(model, [i=1:n], sum(x[i,j,1] + x[i,j,2] for j in islands[i].neighbours) == islands[i].k)

    # 2) No intersections
    for (a,b,c,d) in intersections
        @constraint(model, x[a,b,1] + x[c,d,1] <= 1)
    end

    # 3) Connected bridges

    optimize!(model)

    # TODO: Pretty Print Solution
    # Return runtime
end


"""
Island

# Arguments
- `id`: island id (one-indexed)
- `r`: row containing island
- `c`: column containg island
- `k`: number of bridges to be connected to the island
- `neighbours`: vector of the ids of islands to which a bridge could be build
"""
struct Island
    id::Int64
    r::Int64
    c::Int64
    k::Int64
    neighbours::Vector{Int64}
end


"""
readFile(file)

# Arguments
- `file`: name of the file to read

# Example
islands = readFile("puzzles/sheet.txt")
"""
function readFile(file::String)
    lines = readlines(file)
    nrows, ncols, nislands = parse.(Int, split(lines[1]))
    
    grid = zeros(Int, nrows, ncols)
    islandMap = Dict{Tuple{Int,Int}, Int}()  # (r,c) => id
    islands = Island[]
    intersections = Vector{Tuple{Int64, Int64, Int64, Int64}}()
    
    # Create islands and build position map
    currId = 1
    for r in 1:nrows
        grid[r,:] = parse.(Int, split(lines[r+1]))

        for c in 1:ncols
            if grid[r,c] > 0
                push!(islands, Island(currId, r, c, grid[r,c], Vector{Int64}()))
                islandMap[(r,c)] = currId
                currId += 1
            end
        end
    end
    
    # Find neighbors
    for island in islands
        for (dr, dc) in [(0, -1), (0, 1), (1, 0), (-1, 0)]
            r, c = island.r + dr, island.c + dc
            while 1 <= r <= nrows && 1 <= c <= ncols
                if haskey(islandMap, (r,c))
                    push!(island.neighbours, islandMap[(r,c)])
                    break
                end
                r, c = r + dr, c + dc
            end
        end
    end

    # Find intersections as non-island grid cells with four island neighbours
    for r in 1:nrows
        for c in 1:ncols
            if !haskey(islandMap, (r,c))
                neighIds = Int64[]
                # Look in all four directions
                for (dr, dc) in [(0, -1), (0, 1), (1, 0), (-1, 0)]
                    rr, cc = r + dr, c + dc
                    while 1 <= rr <= nrows && 1 <= cc <= ncols
                        if haskey(islandMap, (rr, cc))
                            push!(neighIds, islandMap[(rr, cc)])
                            break
                        end
                        rr, cc = rr + dr, cc + dc
                    end
                end
                
                # If we found islands in all directions, it's an intersection
                if length(neighIds) == 4
                    push!(intersections, (neighIds[1], neighIds[2], neighIds[3], neighIds[4]))
                end
            end
        end
    end

    return islands, intersections
end
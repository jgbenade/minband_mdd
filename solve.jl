#!julia
using JuMP
using Gurobi

#f = open("myciel4.col")
f = open("testgraph1.5")

counter = 0
lines = readlines(f)

n_vertices = -1
n_edges = -1
adj_mat = zeros(Int64, 1,1)
for l in lines
	
	a = split(l)
	if a[1] == "p"
		n_vertices = int64(a[3])
		n_edges = int64(a[4])
		adj_mat = zeros(Int64, n_vertices, n_vertices)
	end
	if a[1] == "e"
		x = int32(a[2])
		y =int32(a[3])
		
		adj_mat[x,y]=1 
		adj_mat[y,x]=1 
	end

end

println("read info minband\n")
println(n_edges)
println(n_vertices)
println(adj_mat)


m = Model(solver = GurobiSolver(TimeLimit=1200))
@defVar(m,  x[i=1:n_vertices, j=1:n_vertices], Bin)
@defVar(m,  y[i=1:n_vertices], Int)
@defVar(m, z, Int)
	
@setObjective(m, Min, z)

@addConstraint(m, onevalpervar[i=1:n_vertices], sum{x[i,j], j=1:n_vertices} ==1)
@addConstraint(m, onevarperpos[i=1:n_vertices], sum{x[j,i], j=1:n_vertices} ==1)

@addConstraint(m, link[i=1:n_vertices, j=1:n_vertices], y[i] == j*x[i,j])

for from=1:n_vertices-1 
	for to=from+1:n_vertices
		if adj_mat[from, to] == 1 
			println(from, to)
			@addConstraint(m, edge[from,to], z >= y[from]-y[to] )
			@addConstraint(m, edge[to, from], z >= y[to]- y[from])
		end
	end

end


status = solve(m)


println("Objective value: ", getObjectiveValue(m))

xx = getValue(x)
println(xx)
yy = getValue(y)
println(yy)



module Graphs

import Base.+, Base.-, Base.==, Base.in, Base.union, Base.intersect, Base.setdiff

# Export Declarations
export AbstractEdge, AbstractGraph, BK, connectedvertices, DEdge, Δ, DFS, DFS!, Edge, Γ, Graph, incident, isincident, otherend, Vertex

# Abstract Type Declarations
abstract type AbstractEdge end
abstract type AbstractGraph end

# Type Declarations
struct Vertex
	id
end

struct Edge <: AbstractEdge
	p::Tuple{Vertex,Vertex}
	w::Number
end

struct DEdge <: AbstractEdge
	p::Tuple{Any,Any}
	w::Number
end

struct Graph{T<:AbstractEdge} <: AbstractGraph
	V::Set{Vertex}
	E::Set{T}

	# Modifies input variable V to contain any vertex incident to an
	# edge in E but not in the original V.
	function Graph{T}(V::Set{Vertex},E::Set{T}) where T <: AbstractEdge
		union!(V,Set([vertex for edge in E for vertex in edge.p]))
		return new{T}(V,E)
	end
end
Graph(V::Set{Vertex}, E::Set{T}) where {T <: AbstractEdge} = Graph{T}(V,E);

# Graph functions
function union(a::Graph, b::Graph)::Graph
	V = union(a.V,b.V);
	E = union(a.E,b.E);
	return Graph(V,E);
end

function intersect(a::Graph, b::Graph)::Graph
	V = intersect(a.V,b.V);
	E = intersect(a.E,b.E);
	return Graph(V,E);
end

function +(G::Graph{T}, E::T)::Graph{T} where T <: AbstractEdge
	V = union(G.V,E.p);
	E = union(G.E,Set([E]))
	return Graph(V,E)
end

function -(G::Graph{T}, E::T)::Graph{T} where T <: AbstractEdge
	V = G.V;
	E = setdiff(G.E,[E]);
	return Graph(V,E)
end

function DFS!(u::Vertex, G::Graph, Explored::Array=[])::Set{Vertex}
	# Returns all Vertices in Graph G connected to a Vertex u.
	# Only modifies the Explored Array, which should be initially
	# input as empty.
	append!(Explored,[u])
	for edge in incident(u, G.E)
		if (v=otherend(u,edge)) ∉ Explored
			println(u.id, " -> ", v.id)
			DFS!(v,G,Explored)
		end
	end
	return Set(Explored)
end
DFS(u::Vertex,G::Graph) = DFS!(u,G)

# Edge Functions
function ==(a::T, b::T)::Bool where T <: AbstractEdge
	if typeof(a) != typeof(b)
		return false
	elseif typeof(a) == DEdge
		return (a.p == b.p) & (a.w == b.w)
	elseif typeof(a) == Edge
		return (a.p[1] in b.p) & (a.p[2] in b.p) & (a.w == b.w)
	end
end

function in(E::T, iter::Set{T})::Bool where T <: AbstractEdge
	for edge in iter
		if edge == E
			return true
		end
	end
	return false
end

function union(a::Set{T}, b::Set{T})::Set{T} where T <: AbstractEdge
# I feel like there's a better way to do this using filter.
	edges = collect(a);
	for edge in b
		if edge ∉ edges
			append!(edges,[edge])
		end
	end
	return Set(edges)
end

function intersect(a::Set{T}, b::Set{T})::Set{T} where T <: AbstractEdge
	return filter(edge->edge ∈ b, a)
end

function setdiff(a::Set{T}, b::Set{T})::Set{T} where T <: AbstractEdge
	return filter(edge->edge ∉ b, a)
end

function otherend(v::Vertex, e::T)::Vertex where T <: AbstractEdge
	# Return the other vertex in the edge incident to a given vertex.
	if !isincident(v,e)
		throw("Edge is not incident to the vertex!")
	else
		if typeof(e) == Edge
			if e.p[1] == v
            	v2 = e.p[2]
        	else
            	v2 = e.p[1]
        	end
		elseif typeof(e) == DEdge
			v2 = e.p[2]
		end
	end
	return v2
end

# This could fix my set of edges issue
#function Set(edges::Array{Edge,1})
#	edgeset = [];
#	for edge in edges
#		if edge ∉ edgeset
#			append!(edgeset,[edge]);
#		end
#	end
#	return Set(edgeset)
#end


# Vertex Functions
function isincident(v::Vertex, e::T)::Bool where T <: AbstractEdge
	# Takes as arguments a vertex and an edge and returns
	# whether the edge is incident to the vertex.
	if typeof(e) == Edge
		return v in e.p
	elseif typeof(e) == DEdge
		return v == e.p[1]
	end
end

function incident(v::Vertex, edges::Set{T})::Set{T} where T <: AbstractEdge
	# Return a set of all edges in the given set of edges incident to the
	# given vertex.
	return filter(edge->isincident(v, edge), edges);
end

function connectedvertices(V::Set{Vertex}, G::Graph)::Set{Vertex}
	# Takes a set of vertices V and returns a set of all vertices v such that
	# an edge (u,v) exists for a vertex u in V.
	connected::Set{Vertex} = Set();
	for v in V
		edges = incident(v, G.E);
		vertices = [otherend(v,edge) for edge in edges];
		union!(connected,vertices);
	end
	return connected
end
Γ(V::Set{Vertex}, G::Graph) = connectedvertices(V,G);

function +(a::Set{Vertex}, b::Set{Vertex})::Set{Vertex}
	return union(a,b)
end

function -(a::Set{Vertex}, b::Set{Vertex})::Set{Vertex}
	return setdiff(a,b)
end

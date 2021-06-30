#pragma once

#include "DMesh3.h"
#include "g3types.h"
#include <limits>
#include <string>

namespace g3 {

// find the vtx that is the same in both ev0 and ev1
static int find_shared_edge_v(const Vector2i &ev0, const Vector2i ev1) {
	if (ev0.x == ev1.x) {
		return ev0.x;
	} else if (ev0.x == ev1.y) {
		return ev0.x;
	} else if (ev0.y == ev1.x) {
		return ev0.y;
	} else if (ev0.y == ev1.y) {
		return ev0.y;
	} else {
		return DMesh3::InvalidID;
	}
}

static int find_edge_other_v(Vector2i ev, int v) {
	if (ev.x == v) {
		return ev.y;
	} else if (ev.y == v) {
		return ev.x;
	} else {
		return DMesh3::InvalidID;
	}
}

class EdgeLoop;
typedef std::shared_ptr<EdgeLoop> EdgeLoopPtr;

/// <summary>
/// Sequential set of vertices/edges in a mesh, that form a closed loop.
///
/// If all you have are the vertices, use EdgeLoop.VertexLoopToEdgeLoop() to construct an EdgeLoop
/// </summary>
class EdgeLoop {
public:
	DMesh3Ptr Mesh;

	std::vector<int> Vertices;
	std::vector<int> Edges;
	std::vector<int> BowtieVertices; // this may not be initialized!

	EdgeLoop(DMesh3Ptr mesh) {
		Mesh = mesh;
	}
	EdgeLoop(DMesh3Ptr mesh, std::vector<int> vertices, std::vector<int> edges) {
		Mesh = mesh;
		Vertices = vertices;
		Edges = edges;
	}

	EdgeLoop(EdgeLoopPtr copy) {
		Mesh = copy->Mesh;
		Vertices = copy->Vertices;
		Edges = copy->Edges;
		if (copy->BowtieVertices.size()) {
			BowtieVertices = copy->BowtieVertices;
		}
	}

	/// <summary>
	/// construct EdgeLoop from a list of edges of mesh
	/// </summary>
	EdgeLoopPtr FromEdges(DMesh3Ptr mesh, std::list<int> edges) {
		Edges.clear();
		for (int edge : edges) {
			Edges.push_back(edge);
		}

		std::vector<int> Vertices;
		Vertices.resize(Edges.size());
		Vector2i start_ev = mesh->GetEdgeV(Edges[0]);
		Vector2i prev_ev = start_ev;
		for (int i = 1; i < Edges.size(); ++i) {
			Vector2i next_ev = mesh->GetEdgeV(Edges[i % Edges.size()]);
			Vertices[i] = find_shared_edge_v(prev_ev, next_ev);
			prev_ev = next_ev;
		}
		Vertices[0] = find_edge_other_v(start_ev, Vertices[1]);
		return std::make_shared<EdgeLoop>(mesh, Vertices, Edges);
	}

	/// <summary>
	/// construct EdgeLoop from a list of vertices of mesh
	/// </summary>
	EdgeLoopPtr FromVertices(DMesh3Ptr mesh, std::list<int> vertices) {
		int NV = vertices.size();
		std::vector<int> Vertices;
		Vertices.clear();
		for (int vertex : vertices) {
			Vertices.push_back(vertex);
		}
		int NE = NV;
		std::vector<int> Edges;
		Edges.resize(NE);
		for (int i = 0; i < NE; ++i) {
			Edges[i] = mesh->FindEdge(Vertices[i], Vertices[(i + 1) % NE]);
			if (Edges[i] == DMesh3::InvalidID) {
				throw std::logic_error("EdgeLoop.FromVertices: vertices are not connected by edge!");
			}
		}
		return std::make_shared<EdgeLoop>(mesh, Vertices, Edges);
	}

	/// <summary>
	/// construct EdgeLoop from a list of vertices of mesh
	/// if loop is a boundary edge, we can correct orientation if requested
	/// </summary>
	EdgeLoopPtr FromVertices(DMesh3Ptr mesh, std::list<int> vertices, bool bAutoOrient = true) {
		std::vector<int> Vertices;
		for (int vertex : vertices) {
			Vertices.push_back(vertex);
		}

		if (bAutoOrient) {
			int a = Vertices[0], b = Vertices[1];
			int eid = mesh->FindEdge(a, b);
			if (mesh->IsBoundaryEdge(eid)) {
				Vector2i ev = mesh->GetOrientedBoundaryEdgeV(eid);
				if (ev.x == b && ev.y == a) {
					std::reverse(Vertices.begin(), Vertices.end());
				}
			}
		}

		std::vector<int> Edges;
		Edges.resize(Vertices.size());
		for (int i = 0; i < Edges.size(); ++i) {
			int a = Vertices[i], b = Vertices[(i + 1) % Vertices.size()];
			Edges[i] = mesh->FindEdge(a, b);
			if (Edges[i] == DMesh3::InvalidID) {
				throw std::logic_error("EdgeLoop.FromVertices: invalid edge [" + std::to_string(a) + ',' + std::to_string(b) + ']');
			}
		}

		return std::make_shared<EdgeLoop>(mesh, Vertices, Edges);
	}

	int VertexCount() {
		return Vertices.size();
	}
	int EdgeCount() {
		return Edges.size();
	}

	Vector3d GetVertex(int i) {
		return Mesh->GetVertex(Vertices[i]);
	}

	AABB GetBounds() {
		AABB box;
		for (int i = 0; i < Vertices.size(); ++i) {
			box.Contain(Mesh->GetVertex(Vertices[i]));
		}
		return box;
	}

	//DCurve3 ToCurve(DMesh3Ptr sourceMesh)
	//{
	//	if (sourceMesh == null)
	//		sourceMesh = Mesh;
	//	DCurve3 curve = ExtractLoopV(sourceMesh, Vertices);
	//	curve.Closed = true;
	//	return curve;
	//}

	//	/// <summary>
	//	/// if this is a border edge-loop, we can check that it is oriented correctly, and
	//	/// if not, reverse it.
	//	/// Returns true if we reversed orientation.
	//	/// </summary>
	//	public bool CorrectOrientation()
	//	{
	//		int a = Vertices[0], b = Vertices[1];
	//		int eid = Mesh.FindEdge(a, b);
	//		if (Mesh.IsBoundaryEdge(eid)) {
	//			Vector2i ev = Mesh.GetOrientedBoundaryEdgeV(eid);
	//			if (ev.a == b && ev.b == a) {
	//				Reverse();
	//				return true;
	//			}
	//		}
	//		return false;
	//	}

	void Reverse() {
		std::reverse(Vertices.begin(), Vertices.end());
		std::reverse(Edges.begin(), Edges.end());
	}

	/// <summary>
	/// check if all edges of this loop are internal edges (ie none on boundary)
	/// </summary>
	/// <returns></returns>
	bool IsInternalLoop() {
		int NV = Vertices.size();
		for (int i = 0; i < NV; ++i) {
			int eid = Mesh->FindEdge(Vertices[i], Vertices[(i + 1) % NV]);
			if (eid == DMesh3::InvalidID) {
				continue;
			}
			if (Mesh->IsBoundaryEdge(eid)) {
				return false;
			}
		}
		return true;
	}

	/// <summary>
	/// Check if all edges of this loop are boundary edges.
	/// If testMesh != null, will check that mesh instead of internal Mesh
	/// </summary>
	bool IsBoundaryLoop(DMesh3Ptr testMesh) {
		DMesh3Ptr useMesh = !testMesh ? testMesh : Mesh;

		int NV = Vertices.size();
		for (int i = 0; i < NV; ++i) {
			int eid = useMesh->FindEdge(Vertices[i], Vertices[(i + 1) % NV]);
			if (eid == DMesh3::InvalidID) {
				continue;
			}
			if (useMesh->IsBoundaryEdge(eid) == false) {
				return false;
			}
		}
		return true;
	}

	/// <summary>
	/// find index of vertex vID in Vertices list, or -1 if not found
	/// </summary>
	size_t FindVertexIndex(int vID) {
		size_t N = Vertices.size();
		for (size_t i = 0; i < N; ++i) {
			if (Vertices[i] == vID) {
				return i;
			}
		}
		return -1;
	}

	/// <summary>
	/// find index of vertices of loop that is closest to point v
	/// </summary>
	size_t FindNearestVertex(Vector3d v) {
		size_t iNear = -1;
		double fNearSqr = std::numeric_limits<double>::max();
		size_t N = Vertices.size();
		for (size_t i = 0; i < N; ++i) {
			Vector3d lv = Mesh->GetVertex(Vertices[i]);
			double d2 = (lv - v).squaredNorm();
			if (d2 < fNearSqr) {
				fNearSqr = d2;
				iNear = i;
			}
		}
		return iNear;
	}

	// count # of vertices in loop that are within tol of v
	// final param returns last encountered index within tolerance, or -1 if return is 0
	int CountWithinTolerance(Vector3d v, double tol, int &last_in_tol) {
		last_in_tol = -1;
		int count = 0;
		int N = Vertices.size();
		for (int i = 0; i < N; ++i) {
			Vector3d lv = Mesh->GetVertex(Vertices[i]);
			if ((lv - v).norm() < tol) {
				count++;
				last_in_tol = i;
			}
		}
		return count;
	}

	// for ((i++) % N)-type loops, but where we might be using (i--)
	static int WrapSignedIndex(int val, int mod) {
		while (val < 0) {
			val += mod;
		}
		return val % mod;
	}

	// Check if Loop2 is the same set of positions on another mesh.
	// Does not require the indexing to be the same
	// Currently doesn't handle loop-reversal
	bool IsSameLoop(EdgeLoopPtr Loop2, bool bReverse2 = false, double tolerance = 1e-08) {
		// find a duplicate starting vertex
		int N = Vertices.size();
		int N2 = Loop2->Vertices.size();
		if (N != N2)
			return false;

		DMesh3Ptr Mesh2 = Loop2->Mesh;

		int start_i = 0, start_j = -1;

		// try to find a unique same-vertex on each loop. Do not
		// use vertices that have duplicate positions.
		bool bFoundGoodStart = false;
		while (!bFoundGoodStart && start_i < N) {
			Vector3d start_v = Mesh->GetVertex(start_i);
			int count = Loop2->CountWithinTolerance(start_v, tolerance, start_j);
			if (count == 1) {
				bFoundGoodStart = true;
			} else {
				start_i++;
			}
		}
		if (!bFoundGoodStart)
			return false; // no within-tolerance duplicate vtx to start at

		for (int ii = 0; ii < N; ++ii) {
			int i = (start_i + ii) % N;
			int j = (bReverse2) ?
							  WrapSignedIndex(start_j - ii, N2) :
							  (start_j + ii) % N2;
			Vector3d v = Mesh->GetVertex(Vertices[i]);
			Vector3d v2 = Mesh2->GetVertex(Loop2->Vertices[j]);
			if ((v2 - v).norm() > tolerance)
				return false;
		}

		return true;
	}

	/// <summary>
	/// stores vertices [starti, starti+1, ... starti+count-1] in span, and returns span, or null if invalid range
	/// </summary>
	std::vector<int> GetVertexSpan(int starti, int count, std::vector<int> span, bool reverse = false) {
		int N = Vertices.size();
		if (starti < 0 || starti >= N || count > N - 1)
			return std::vector<int>();
		if (reverse) {
			for (int k = 0; k < count; ++k)
				span[count - k - 1] = Vertices[(starti + k) % N];
		} else {
			for (int k = 0; k < count; ++k)
				span[k] = Vertices[(starti + k) % N];
		}
		return span;
	}

	/// <summary>
	/// Convert a vertex loop to an edge loop. This should be somewhere else...
	/// </summary>
	std::vector<int> VertexLoopToEdgeLoop(DMesh3Ptr mesh, std::vector<int> vertex_loop) {
		int NV = vertex_loop.size();
		std::vector<int> edges;
		edges.resize(NV);
		for (int i = 0; i < NV; ++i) {
			int v0 = vertex_loop[i];
			int v1 = vertex_loop[(i + 1) % NV];
			edges[i] = mesh->FindEdge(v0, v1);
		}
		return edges;
	}
};
} // namespace g3
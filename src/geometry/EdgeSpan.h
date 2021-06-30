#pragma once

#include "g3types.h"
#include "EdgeLoop.h"
#include "DCurve3.h"

#include <cfloat>
#include <vector>
#include <list>

#include "core/config/engine.h"

namespace g3 {
/// <summary>
/// An EdgeSpan is a continous set of edges in a Mesh that is *not* closed
/// (that would be an EdgeLoop)
/// </summary>
class EdgeSpan {
public:
	DMesh3Ptr Mesh;

	std::vector<int> Vertices;
	std::vector<int> Edges;

	std::vector<int> BowtieVertices; // this may not be initialized!

	EdgeSpan(DMesh3Ptr mesh) {
		Mesh = mesh;
	}

	EdgeSpan(DMesh3Ptr mesh, std::vector<int> vertices, std::vector<int> edges) {
		Mesh = mesh;
		Vertices = vertices;
		Edges = edges;
	}

	/// <summary>
	/// construct EdgeSpan from a list of edges of mesh
	/// </summary>
	static EdgeSpan FromEdges(DMesh3Ptr mesh, std::list<int> edges) {
		std::vector<int> Edges;
		for (int edge : edges) {
			Edges.push_back(edge);
        }
		std::vector<int> Vertices;
        Vertices.resize(Edges.size() + 1);
		Vector2i start_ev = mesh->GetEdgeV(Edges[0]);
		Vector2i prev_ev = start_ev;
		if (Edges.size() > 1) {
			for (int i = 1; i < Edges.size(); ++i) {
				Vector2i next_ev = mesh->GetEdgeV(Edges[i]);
				Vertices[i] = find_shared_edge_v(prev_ev, next_ev);
				prev_ev = next_ev;
			}
			Vertices[0] = find_edge_other_v(start_ev, Vertices[1]);
			Vertices[Vertices.size() - 1] = find_edge_other_v(prev_ev, Vertices[Vertices.size() - 2]);
		} else {
			Vertices[0] = start_ev[0];
			Vertices[1] = start_ev[1];
		}
		return EdgeSpan(mesh, Vertices, Edges, false);
	}

	/// <summary>
	/// construct EdgeSpan from a list of vertices of mesh
	/// </summary>
	static EdgeSpan FromVertices(DMesh3Ptr mesh, std::list<int> vertices) {
		int NV = vertices.size();
		std::vector<int> Vertices;
		Vertices.resize(NV);
        std::list<int>::iterator it = vertices.begin();
		for (int i = 0; i < NV; ++i) {
            std::advance(it, 1);
			Vertices[i] = *it;
		}
		int NE = NV - 1;
		std::vector<int> Edges;
		Edges.resize(NE);
		for (int i = 0; i < NE; ++i) {
			Edges[i] = mesh->FindEdge(Vertices[i], Vertices[i+1]);
			ERR_CONTINUE_MSG(Edges[i] == DMesh3::InvalidID, "EdgeSpan.FromVertices: vertices are not connected by edge!");
		}
		return EdgeSpan(mesh, Vertices, Edges);
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
		for (int i = 0; i < Vertices.size(); ++i)
			box.Contain(Mesh->GetVertex(Vertices[i]));
		return box;
	}

	DCurve3Ptr ToCurve(DMesh3Ptr sourceMesh = nullptr) {
		if (sourceMesh == nullptr) {
			sourceMesh = Mesh;
		}
		DCurve3Ptr curve = ExtractLoopV(sourceMesh, Vertices);
		curve->SetClosed(false);
		return curve;
	}

	bool IsInternalSpan() {
		int NV = Vertices.size();
		for (int i = 0; i < NV - 1; ++i) {
			int eid = Mesh->FindEdge(Vertices[i], Vertices[i + 1]);
			if (eid == DMesh3::InvalidID) {
                continue;
            }
			if (Mesh->IsBoundaryEdge(eid)) {
				return false;
            }
		}
		return true;
	}

	bool IsBoundarySpan(DMesh3Ptr testMesh = nullptr) {
		
		DMesh3Ptr useMesh = Mesh;
		if (testMesh) {
			useMesh = testMesh;
		}

		int NV = Vertices.size();
		for (int i = 0; i < NV - 1; ++i) {
			int eid = useMesh->FindEdge(Vertices[i], Vertices[i + 1]);
			ERR_CONTINUE(eid == DMesh3::InvalidID);
			if (useMesh->IsBoundaryEdge(eid) == false) {
				return false;
			}
		}
		return true;
	}

	int FindNearestVertex(Vector3d v) {
		int iNear = -1;
		double fNearSqr = DBL_MAX;
		int N = Vertices.size();
		for (int i = 0; i < N; ++i) {
			Vector3d lv = Mesh->GetVertex(Vertices[i]);
			double d2 = v.DistanceSquared(lv);
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

	// Check if Spanw is the same set of positions on another mesh.
	// Does not require the indexing to be the same
	bool IsSameSpan(EdgeSpan Spanw, bool bReverse2 = false, double tolerance = Wm5::ZERO_TOLERANCE) {
		// [RMS] this is much easier than for a loop, because it has to have
		//   same endpoints. But don't have time right now.
		ERR_FAIL_MSG_V("todo!", false);
	}

	// /// <summary>
	// /// Exhaustively check that verts and edges of this EdgeSpan are consistent. Not for production use.
	// /// </summary>
	// bool CheckValidity(FailMode eFailMode = FailMode.Throw)
	// {
	//     bool is_ok = true;
	//     Action<bool> CheckOrFailF = (b) => { is_ok = is_ok && b; };
	//     if (eFailMode == FailMode.DebugAssert) {
	//         CheckOrFailF = (b) => { Debug.Assert(b); is_ok = is_ok && b; };
	//     } else if (eFailMode == FailMode.gDevAssert) {
	//         CheckOrFailF = (b) => { Util.gDevAssert(b); is_ok = is_ok && b; };
	//     } else if (eFailMode == FailMode.Throw) {
	//         CheckOrFailF = (b) => { if (b == false) throw new Exception("EdgeSpan.CheckValidity: check failed"); };
	//     }

	//     CheckOrFailF(Vertices.Length == Edges.Length + 1);
	//     for (int ei = 0; ei < Edges.Length; ++ei) {
	//         Vector2i ev = Mesh.GetEdgeV(Edges[ei]);
	//         CheckOrFailF(Mesh.IsVertex(ev.a));
	//         CheckOrFailF(Mesh.IsVertex(ev.b));
	//         CheckOrFailF(Mesh.FindEdge(ev.a, ev.b) != DMesh3.InvalidID);
	//         CheckOrFailF(Vertices[ei] == ev.a || Vertices[ei] == ev.b);
	//         CheckOrFailF(Vertices[ei + 1] == ev.a || Vertices[ei + 1] == ev.b);
	//     }
	//     for (int vi = 0; vi < Vertices.Length-1; ++vi) {
	//         int a = Vertices[vi], b = Vertices[vi + 1];
	//         CheckOrFailF(Mesh.IsVertex(a));
	//         CheckOrFailF(Mesh.IsVertex(b));
	//         CheckOrFailF(Mesh.FindEdge(a, b) != DMesh3.InvalidID);
	//         if (vi < Vertices.Length - 2) {
	//             int n = 0, edge_before_b = Edges[vi], edge_after_b = Edges[vi + 1];
	//             for (int nbr_e : Mesh->VtxEdgesItr(b)) {
	//                 if (nbr_e == edge_before_b || nbr_e == edge_after_b)
	//                     n++;
	//             }
	//             CheckOrFailF(n == 2);
	//         }
	//     }
	//     return true;
	// }

	/// <summary>
	/// Convert vertex span to list of edges. This should be somewhere else.
	/// </summary>
	static std::vector<int> VerticesToEdges(DMesh3Ptr mesh, std::vector<int> vertex_span) {
		int NV = vertex_span.size();
		std::vector<int> edges;
		edges.resize(NV - 1);
		for (int i = 0; i < NV - 1; ++i) {
			int v0 = vertex_span[i];
			int v1 = vertex_span[(i + 1)];
			edges[i] = mesh->FindEdge(v0, v1);
		}
		return edges;
	}
};
} // namespace g3
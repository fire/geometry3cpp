#pragma once

#include "EdgeLoop.h"
#include "EdgeSpan.h"

#include <list>
#include <vector>
#include <algorithm>

namespace g3 {

enum class MeshBoundaryLoopsException {
	UnclosedLoop,
	BowtieFailure,
	RepeatedEdge,
};

class MeshBoundaryLoops;
typedef std::shared_ptr<MeshBoundaryLoops> MeshBoundaryLoopsPtr;

/// Extract boundary EdgeLoops from Mesh. Can also extract EdgeSpans for open areas,
/// however default behavior is to ignore these. Set .SpanBehavior to configure.
class MeshBoundaryLoops : public EdgeLoop {
public:
	DMesh3Ptr Mesh;
	std::list<EdgeLoopPtr> Loops;

	std::list<EdgeSpan> Spans; // spans are unclosed loops
	bool SawOpenSpans = false; // will be set to true if we find any open spans
	bool FellBackToSpansOnFailure = false; // set to true if we had to add spans to recover from failure
	// currently this happens if we cannot extract simple loops from a loop with bowties

	// What should we do if we encounter open spans. Mainly a result of EdgeFilter, but can also
	// happen on meshes w/ crazy bowties
	enum class SpanBehaviors {
		Ignore,
		ThrowException,
		Compute
	};
	SpanBehaviors SpanBehavior = SpanBehaviors::Compute;

	// What should we do if we encounter an unrecoverable failure while walking a loop
	enum class FailureBehaviors {
		ThrowException, // die, and you clean up
		ConvertToOpenSpan // keep un-closed loop as a span
	};
	FailureBehaviors FailureBehavior = FailureBehaviors::ConvertToOpenSpan;

	//// if enabled, only edges where this returns true are considered
	//Func<int, bool> EdgeFilterF = null;

	//// if we throw an exception, we will try to set FailureBowties, so that client
	//// can try repairing these vertices
	//std::list<int> FailureBowties = null;

	MeshBoundaryLoops(DMesh3Ptr mesh, bool bAutoCompute = true) :
			EdgeLoop(mesh) {
		if (bAutoCompute) {
			Compute();
		}
	}

	size_t Count() {
		return Loops.size();
	}
	int SpanCount() {
		return Spans.size();
	}

	/// <summary>
	/// Index of Loop with largest vertex count
	/// </summary>
	int MaxVerticesLoopIndex() {
		int j = 0;
		for (int i = 1; i < Loops.size(); ++i) {
            std::list<EdgeLoopPtr>::iterator it = Loops.begin();
            std::advance(it, i);
            std::list<EdgeLoopPtr>::iterator itj = Loops.begin();
            std::advance(itj, j);
			if ((*it)->Vertices.size() > (*itj)->Vertices.size()){
                j = i;            
            }
		}
		return j;
	}

	/// <summary>
	/// find pair (loop_index,in_loop_index) of vertex vID in EdgeLoops, or Vector2i.Max if not found
	/// </summary>
	Vector2i FindVertexIndex(int vID) {
		int N = Loops.size();
		for (int li = 0; li < N; ++li) {
            std::list<EdgeLoopPtr>::iterator ili = Loops.begin();
            std::advance(ili, li);
			int idx = (*ili)->FindVertexIndex(vID);
			if (idx >= 0)
				return Vector2i(li, idx);
		}
		return Vector2i::Max;
	}

	/// <summary>
	/// find the loop index that contains a vertex, or return -1
	/// </summary>
	int FindLoopContainingVertex(int vid) {
		int N = Loops.size();
		for (int li = 0; li < N; ++li) {
            std::list<EdgeLoopPtr>::iterator it = Loops.begin();
            std::advance(it, li);

            std::vector<int>::iterator found;
            found = std::find((*it)->Vertices.begin(), (*it)->Vertices.end());
			if (found != (*it)->Vertices.end()) {
				return *found;
            }
		}
		return -1;
	}

	/// <summary>
	/// find the loop index that contains an edge, or return -1
	/// </summary>
	int FindLoopContainingEdge(int eid) {
		int N = Loops.size();
        std::list<EdgeLoopPtr>::iterator it = Loops.begin();
        std::advance(it, li);
		for (int li = 0; li < N; ++li) {
			if (Loops[li].Edges.Contains(eid)) {
				return li;
            }
		}
		return -1;
	}

	/// <summary>
	/// Find the set of boundary EdgeLoops. Note that if we encounter topological
	/// issues, we will throw MeshBoundaryLoopsException w/ more info (if possible)
	/// </summary>
	bool Compute() {
		// This algorithm assumes that triangles are oriented consistently,
		// so closed boundary-loop can be followed by walking edges in-order

		Loops.clear();
		Spans.clear();

		// early-out if we don't actually have boundaries
		if (Mesh->CachedIsClosed()) {
			return true;
		}

		constexpr int NE = Mesh->MaxEdgeID();

		// Temporary memory used to indicate when we have "used" an edge.
		std::vector<bool> used_edge;
		used_edge.resize(NE);

		// current loop is stored here, cleared after each loop extracted
		std::list<int> loop_edges; // [RMS] not sure we need this...
		std::list<int> loop_verts;
		std::list<int> bowties;

		// Temp buffer for reading back all boundary edges of a vertex.
		// probably always small but in pathological cases it could be large...
		std::vector<int> all_e;
		all_e.resize(16);

		// [TODO] might make sense to precompute some things here, like num_be for each bdry vtx?

		// process all edges of mesh
		for (int eid = 0; eid < NE; ++eid) {
			if (!Mesh->IsEdge(eid))
				continue;
			if (used_edge[eid] == true)
				continue;
			if (Mesh->IsBoundaryEdge(eid) == false)
				continue;

			if (EdgeFilterF != nullptr && EdgeFilterF(eid) == false) {
				used_edge[eid] = true;
				continue;
			}

			// ok this is start of a boundary chain
			int eStart = eid;
			used_edge[eStart] = true;
			loop_edges.push_back(eStart);

			int eCur = eid;

			// follow the chain in order of oriented edges
			bool bClosed = false;
			bool bIsOpenSpan = false;
			while (!bClosed) {
				Vector2i ev = Mesh.GetOrientedBoundaryEdgeV(eCur);
				int cure_a = ev.a, cure_b = ev.b;
				if (bIsOpenSpan) {
					cure_a = ev.b;
					cure_b = ev.a;
				} else {
					loop_verts.push_back(cure_a);
				}

				int e0 = -1, e1 = 1;
				int bdry_nbrs = Mesh.VtxBoundaryEdges(cure_b, e0, e1);

				// have to filter this list, if we are filtering. this is ugly.
				if (EdgeFilterF != nullptr) {
					if (bdry_nbrs > 2) {
						if (bdry_nbrs >= all_e.Length)
                            std::vector<int> new_all_e;
                            new_all_e.resize(bdry_nbrs);
							all_e = new_all_e;
						// we may repeat this below...irritating...
						int num_be = Mesh->VtxAllBoundaryEdges(cure_b, all_e);
						num_be = BufferUtil.CountValid(all_e, EdgeFilterF, num_be);
					} else {
						if (EdgeFilterF(e0) == false)
							bdry_nbrs--;
						if (EdgeFilterF(e1) == false)
							bdry_nbrs--;
					}
				}

				if (bdry_nbrs < 2) { // hit an 'endpoint' vertex (should only happen when Filter is on...)
					if (SpanBehavior == SpanBehaviors.ThrowException)
						throw new MeshBoundaryLoopsException("MeshBoundaryLoops.Compute: found open span at vertex " + cure_b){ UnclosedLoop = true };
					if (bIsOpenSpan) {
						bClosed = true;
						continue;
					} else {
						bIsOpenSpan = true; // begin open span
						eCur = loop_edges[0]; // restart at other end of loop
						loop_edges.Reverse(); // do this so we can push to front
						continue;
					}
				}

				int eNext = -1;

				if (bdry_nbrs > 2) {
					// found "bowtie" vertex...things just got complicated!

					if (cure_b == loop_verts[0]) {
						// The "end" of the current edge is the same as the start vertex.
						// This means we can close the loop here. Might as well!
						eNext = -2; // sentinel value used below

					} else {
						// try to find an unused outgoing edge that is oriented properly.
						// This could create sub-loops, we will handle those later
						if (bdry_nbrs >= all_e.Length)
							all_e = new int[2 * bdry_nbrs];
						int num_be = Mesh.VtxAllBoundaryEdges(cure_b, all_e);
						ERR_FAIL_COND(num_be != bdry_nbrs);

						if (EdgeFilterF != nullptr) {
							num_be = BufferUtil.FilterInPlace(all_e, EdgeFilterF, num_be);
						}

						// Try to pick the best "turn left" vertex.
						eNext = find_left_turn_edge(eCur, cure_b, all_e, num_be, used_edge);
						if (eNext == -1) {
							if (FailureBehavior == FailureBehaviors.ThrowException || SpanBehavior == SpanBehaviors.ThrowException)
								throw new MeshBoundaryLoopsException("MeshBoundaryLoops.Compute: cannot find valid outgoing edge at bowtie vertex " + cure_b){ BowtieFailure = true };

							// ok, we are stuck. all we can do now is terminate this loop and keep it as a span
							if (bIsOpenSpan) {
								bClosed = true;
							} else {
								bIsOpenSpan = true;
								bClosed = true;
							}
							continue;
						}
					}

					if (bowties.Contains(cure_b) == false)
						bowties.Add(cure_b);

				} else {
					// walk forward to next available edge
					Debug.Assert(e0 == eCur || e1 == eCur);
					eNext = (e0 == eCur) ? e1 : e0;
				}

				if (eNext == -2) {
					// found a bowtie vert that is the same as start-of-loop, so we
					// are just closing it off explicitly
					bClosed = true;
				} else if (eNext == eStart) {
					// found edge at start of loop, so loop is done.
					bClosed = true;
				} else if (used_edge[eNext] != false) {
					// disaster case - the next edge is already used, but it is not the start of our loop
					// All we can do is convert to open span and terminate
					if (FailureBehavior == FailureBehaviors.ThrowException || SpanBehavior == SpanBehaviors.ThrowException)
						throw new MeshBoundaryLoopsException("MeshBoundaryLoops.Compute: encountered repeated edge " + eNext){ RepeatedEdge = true };
					bIsOpenSpan = true;
					bClosed = true;

				} else {
					// push onto accumulated list
					if (used_edge[eNext] != false) {
						continue;
					}
					loop_edges.Add(eNext);
					used_edge[eNext] = true;
					eCur = eNext;
				}
			}

			if (bIsOpenSpan) {
				SawOpenSpans = true;
				if (SpanBehavior == SpanBehaviors.Compute) {
					std::reverse(loop_edges.begin(), loop_edges.end()); // orient properly
					EdgeSpan span = EdgeSpan::FromEdges(Mesh, loop_edges);
					Spans.Add(span);
				}
			} else if (bowties.Count > 0) {
				// if we saw a bowtie vertex, we might need to break up this loop,
				// so call extract_subloops
				Subloops subloops = extract_subloops(loop_verts, loop_edges, bowties);
				for (var loop : subloops.Loops) {
					Loops.Add(loop);
				}
				if (subloops.Spans.Count > 0) {
					FellBackToSpansOnFailure = true;
					foreach (var span in subloops.Spans)
						Spans.Add(span);
				}
			} else {
				// clean simple loop, convert to EdgeLoop instance
				EdgeLoop loop = new EdgeLoop(Mesh);
				loop.Vertices = loop_verts.ToArray();
				loop.Edges = loop_edges.ToArray();
				Loops.Add(loop);
			}

			// reset these lists
			loop_edges.Clear();
			loop_verts.Clear();
			bowties.Clear();
		}

		return true;
	}

	// [TODO] cache this in a dictionary? we will not need very many, but we will
	//   need each multiple times!
	Vector3d get_vtx_normal(int vid) {
		Vector3d n = Vector3d.Zero;
		for (int ti : Mesh->VtxTrianglesItr(vid)) {
			n += Mesh->GetTriNormal(ti);
		}
		n.normalize();
		return n;
	}

	// ok, bdry_edges[0...bdry_edges_count] contains the boundary edges coming out of bowtie_v.
	// We want to pick the best one to continue the loop that came in to bowtie_v on incoming_e.
	// If the loops are all sane, then we will get the smallest loops by "turning left" at bowtie_v.
	// So, we compute the tangent plane at bowtie_v, and then the signed angle for each
	// viable edge in this plane.
	//
	// [TODO] handle degenerate edges. what do we do then? Currently will only chose
	//  degenerate edge if there are no other options (I think...)
	int find_left_turn_edge(int incoming_e, int bowtie_v, std::vector<int> bdry_edges, int bdry_edges_count, std::vector<bool> used_edges) {
		// compute normal and edge [a,bowtie]
		Vector3d n = get_vtx_normal(bowtie_v);
		int other_v = Mesh.edge_other_v(incoming_e, bowtie_v);
		Vector3d ab = Mesh.GetVertex(bowtie_v) - Mesh.GetVertex(other_v);

		// our winner
		int best_e = -1;
		double best_angle = double.MaxValue;

		for (int i = 0; i < bdry_edges_count; ++i) {
			int bdry_eid = bdry_edges[i];
			if (used_edges[bdry_eid] == true)
				continue; // this edge is already used
			Vector2i bdry_ev = Mesh.GetOrientedBoundaryEdgeV(bdry_eid);
			if (bdry_ev.a != bowtie_v)
				continue; // have to be able to chain to end of current edge, orientation-wise

			// compute projected angle
			Vector3d bc = Mesh.GetVertex(bdry_ev.b) - Mesh.GetVertex(bowtie_v);
			float fAngleS = MathUtil.PlaneAngleSignedD((Vector3f)ab, (Vector3f)bc, (Vector3f)n);

			// turn left!
			if (best_angle == double.MaxValue || fAngleS < best_angle) {
				best_angle = fAngleS;
				best_e = bdry_eid;
			}
		}

		// [RMS] w/ bowtie vertices and open spans, this does happen
		//Debug.Assert(best_e != -1);

		return best_e;
	}

	struct Subloops {
		std::list<EdgeLoop> Loops;
		std::list<EdgeSpan> Spans;
	};

	// This is called when loopV contains one or more "bowtie" vertices.
	// These vertices *might* be duplicated in loopV (but not necessarily)
	// If they are, we have to break loopV into subloops that don't contain duplicates.
	//
	// The list bowties contains all the possible duplicates
	// (all v in bowties occur in loopV at least once)
	//
	// Currently loopE is not used, and the returned EdgeLoop objects do not have their Edges
	// arrays initialized. Perhaps to improve in future.
	//
	// An unhandled case to think about is where we have a sequence [..A..B..A..B..] where
	// A and B are bowties. In this case there are no A->A or B->B subloops. What should
	// we do here??
	Subloops
	extract_subloops(std::list<int> loopV, std::list<int> loopE, std::list<int> bowties) {
		Subloops subs;

		// figure out which bowties we saw are actually duplicated in loopV
		std::list<int> dupes;
		foreach (int bv in bowties) {
			if (count_in_list(loopV, bv) > 1)
				dupes.Add(bv);
		}

		// we might not actually have any duplicates, if we got luck. Early out in that case
		if (dupes.Count == 0) {
			subs.Loops.Add(new EdgeLoop(Mesh){
					Vertices = loopV.ToArray(), Edges = loopE.ToArray(), BowtieVertices = bowties.ToArray() });
			return subs;
		}

		// This loop extracts subloops until we have dealt with all the
		// duplicate vertices in loopV
		while (dupes.Count > 0) {
			// Find shortest "simple" loop, ie a loop from a bowtie to itself that
			// does not contain any other bowties. This is an independent loop.
			// We're doing a lot of extra work here if we only have one element in dupes...
			int bi = 0, bv = 0;
			int start_i = -1, end_i = -1;
			int bv_shortest = -1;
			int shortest = int.MaxValue;
			for (; bi < dupes.Count; ++bi) {
				bv = dupes[bi];
				if (is_simple_bowtie_loop(loopV, dupes, bv, out start_i, out end_i)) {
					int len = count_span(loopV, start_i, end_i);
					if (len < shortest) {
						bv_shortest = bv;
						shortest = len;
					}
				}
			}

			// failed to find a simple loop. Not sure what to do in this situation.
			// If we don't want to throw, all we can do is convert the remaining
			// loop to a span and return.
			// (Or should we keep it as a loop and set flag??)
			if (bv_shortest == -1) {
				if (FailureBehavior == FailureBehaviors.ThrowException) {
					FailureBowties = dupes;
					throw new MeshBoundaryLoopsException("MeshBoundaryLoops.Compute: Cannot find a valid simple loop");
				}
				EdgeSpan span = EdgeSpan(Mesh);
				std::list<int> verts;
				for (int i = 0; i < loopV.Count; ++i) {
					if (loopV[i] != -1)
						verts.Add(loopV[i]);
				}
				span.Vertices = verts.ToArray();
				span.Edges = EdgeSpan.VerticesToEdges(Mesh, span.Vertices);
				span.BowtieVertices = bowties.ToArray();
				subs.Spans.Add(span);
				return subs;
			}

			if (bv != bv_shortest) {
				bv = bv_shortest;
				// running again just to get start_i and end_i...
				is_simple_bowtie_loop(loopV, dupes, bv, out start_i, out end_i);
			}

			Debug.Assert(loopV[start_i] == bv && loopV[end_i] == bv);

			EdgeLoop loop = new EdgeLoop(Mesh);
			loop.Vertices = extract_span(loopV, start_i, end_i, true);
			loop.Edges = EdgeLoop.VertexLoopToEdgeLoop(Mesh, loop.Vertices);
			loop.BowtieVertices = bowties.ToArray();
			subs.Loops.Add(loop);

			// If there are no more duplicates of this bowtie, we can treat
			// it like a regular vertex now
			if (count_in_list(loopV, bv) < 2)
				dupes.Remove(bv);
		}

		// Should have one loop left that contains duplicates.
		// Extract this as a separate loop
		int nLeft = 0;
		for (int i = 0; i < loopV.Count; ++i) {
			if (loopV[i] != -1)
				nLeft++;
		}
		if (nLeft > 0) {
			EdgeLoop loop = EdgeLoop(Mesh);
            loop.Vertices.clear();
            loop.Vertices.resize(nLeft);
			int vi = 0;
			for (int i = 0; i < loopV.Count; ++i) {
				if (loopV[i] != -1)
					loop.Vertices[vi++] = loopV[i];
			}
			loop.Edges = EdgeLoop.VertexLoopToEdgeLoop(Mesh, loop.Vertices);
			loop.BowtieVertices = bowties.ToArray();
			subs.Loops.Add(loop);
		}

		return subs;
	}

	/*
		* In all the functions below, the list loopV is assumed to possibly
		* contain "removed" vertices indicated by -1. These are ignored.
		*/

	// Check if the loop from bowtieV to bowtieV inside loopV contains any other bowtie verts.
	// Also returns start and end indices in loopV of "clean" loop
	// Note that start may be < end, if the "clean" loop wraps around the end
	bool is_simple_bowtie_loop(std::list<int> loopV, std::list<int> bowties, int bowtieV, int &start_i, int &end_i) {
		// find two indices of bowtie vert
		start_i = find_index(loopV, 0, bowtieV);
		end_i = find_index(loopV, start_i + 1, bowtieV);

		if (is_simple_path(loopV, bowties, bowtieV, start_i, end_i)) {
			return true;

		} else if (is_simple_path(loopV, bowties, bowtieV, end_i, start_i)) {
			int tmp = start_i;
			start_i = end_i;
			end_i = tmp;
			return true;

		} else
			return false; // not a simple bowtie loop!
	}

	// check if forward path from loopV[i1] to loopV[i2] contains any bowtie verts other than bowtieV
	bool is_simple_path(std::list<int> loopV, std::list<int> bowties, int bowtieV, int i1, int i2) {
		int N = loopV.Count;
		for (int i = i1; i != i2; i = (i + 1) % N) {
			int vi = loopV[i];
			if (vi == -1)
				continue; // skip removed vertices
			if (vi != bowtieV && bowties.Contains(vi))
				return false;
		}
		return true;
	}

	// Read out the span from loop[i0] to loop [i1-1] into an array.
	// If bMarkInvalid, then these values are set to -1 in loop
	std::vector<int> extract_span(std::list<int> loop, int i0, int i1, bool bMarkInvalid) {
		int num = count_span(loop, i0, i1);
		std::vector<int> a = new int[num];
		int ai = 0;
		int N = loop.size();
		for (int i = i0; i != i1; i = (i + 1) % N) {
			if (loop[i] != -1) {                
                std::list<EdgeLoopPtr>::iterator it = Loops.begin();
                std::advance(it, i);
				a[ai++] = *it;
				if (bMarkInvalid)
					loop[i] = -1;
			}
		}
		return a;
	}

	// count number of valid vertices in l between loop[i0] and loop[i1-1]
	int count_span(std::list<int> l, int i0, int i1) {
		int c = 0;
		int N = l.Count;
		for (int i = i0; i != i1; i = (i + 1) % N) {
			if (l[i] != -1)
				c++;
		}
		return c;
	}

	// find the index of item in loop, starting at start index
	int find_index(std::list<int> loop, int start, int item) {
		for (int i = start; i < loop.Count; ++i)
			if (loop[i] == item)
				return i;
		return -1;
	}

	// count number of times item appears in loop
	int count_in_list(std::list<int> loop, int item) {
		int c = 0;
		for (int i = 0; i < loop.Count; ++i) {
			if (loop[i] == item) {
				c++;
			}
		}
		return c;
	}
};
} // namespace g3
#pragma once
#include "core/object/ref_counted.h"
#include "scene/resources/mesh.h"

#include "scene/resources/mesh.h"

#include "MeshBoundaryLoop.h"
#include "profile_util.h"
#include "src/geometry/g3types.h"
#include "src/mesh/DMesh3Builder.h"
#include "src/mesh/Remesher.h"
#include "src/spatial/BasicProjectionTargets.h"
#include <DMesh3.h>
#include <DMeshAABBTree3.h>
#include <MeshQueries.h>
#include <MeshSubdivider.h>
#include <VectorUtil.h>
#include <refcount_vector.h>
#include <small_list_set.h>
#include <algorithm>
#include <limits>
#include <list>

#include "../../godot/scene/resources/surface_tool.h"
#include "geometry3_process.h"
#include "src/geometry/MeshBoundaryLoop.h"
#include "src/mesh/DMesh3.h"
#include "src/mesh/MeshConstraints.h"
#include "src/spatial/DCurveProject.h"

namespace g3 {
void PreserveAllBoundaryEdges(g3::MeshConstraintsPtr cons,
		g3::DMesh3Ptr p_mesh) {
	if (!p_mesh) {
		return;
	}
	if (p_mesh->EdgeCount() <= 1) {
		return;
	}
	int32_t max_edge_id = p_mesh->MaxEdgeID();
	for (int edge_id = 0; edge_id < max_edge_id; ++edge_id) {
		if (p_mesh->IsEdge(edge_id) && p_mesh->IsBoundaryEdge(edge_id)) {
			cons->SetOrUpdateEdgeConstraint(edge_id,
					g3::EdgeConstraint::FullyConstrained());

			Vector2i ev = p_mesh->GetEdgeV(edge_id);
			VertexConstraint vc = VertexConstraint::Pinned();
			cons->SetOrUpdateVertexConstraint(ev.x(), vc);
			cons->SetOrUpdateVertexConstraint(ev.y(), vc);
		}
	}
}

// https://github.com/gradientspace/geometry3Sharp/blob/master/mesh/MeshConstraintUtil.cs
/// <summary>
/// Remove 'fin' triangles that have only one connected triangle.
/// Removing one fin can create another, by default will keep iterating
/// until all fins removed (in a not very efficient way!).
/// Pass bRepeatToConvergence=false to only do one pass.
/// [TODO] if we are repeating, construct face selection from numbers of first
// list and iterate over that on future passes!
/// </summary>
size_t RemoveFinTriangles(g3::DMesh3Ptr mesh,
		bool bRepeatToConvergence = true) {
	size_t nRemoved = 0;
	std::list<int> to_remove;
	while (true) {
		for (int tid : mesh->TriangleIndices()) {
			Index3i nbrs = mesh->GetTriNeighbourTris(tid);
			int c = ((nbrs.x() != DMesh3::InvalidID) ? 1 : 0) +
					((nbrs.y() != DMesh3::InvalidID) ? 1 : 0) +
					((nbrs.z() != DMesh3::InvalidID) ? 1 : 0);
			if (c <= 1) {
				to_remove.push_back(tid);
			}
		}
		if (!to_remove.size()) {
			return nRemoved;
		}
		nRemoved += to_remove.size();

		for (int tid : to_remove) {
			mesh->RemoveTriangle(tid, false, true);
		}
		to_remove.clear();
		if (!bRepeatToConvergence) {
			break;
		}
	}
	return nRemoved;
}

// for all vertices in loopV, constrain to target
// for all edges in loopV, disable flips and constrain to target
static void ConstrainVtxLoopToMesh(MeshConstraintsPtr cons, DMesh3Ptr mesh,
		std::list<int> loopV,
		IProjectionTargetPtr target,
		int setID = -1) {
	VertexConstraint vc = VertexConstraint(target);
	int N = loopV.size();
	for (int i = 0; i < N; ++i) {
		std::list<int>::iterator it = loopV.begin();
		std::advance(it, i);
		cons->SetOrUpdateVertexConstraint(*it, vc);
	}

	EdgeConstraint ec = EdgeConstraint(EdgeRefineFlags::NoFlip, target);
	ec.TrackingSetID = setID;
	for (int i = 0; i < N; ++i) {
		std::list<int>::iterator it = loopV.begin();
		std::advance(it, i);
		int v0 = *it;
		it = loopV.begin();
		std::advance(it, (i + 1) % N);
		int v1 = *it;

		int eid = mesh->FindEdge(v0, v1);
		if (eid == DMesh3::InvalidID) {
			return;
		}
		cons->SetOrUpdateEdgeConstraint(eid, ec);
	}
}
static void ConstrainVtxLoopTo(Remesher r, std::vector<int> loopV,
		IProjectionTargetPtr target, int setID = -1) {
	if (r.Constraints() == nullptr) {
		r.SetExternalConstraints(std::make_shared<MeshConstraints>());
	}
	std::list<int> loop_list_verts;
	for (int index : loopV) {
		loop_list_verts.push_back(index);
	}
	ConstrainVtxLoopToMesh(r.Constraints(), r.Mesh(), loop_list_verts, target);
}

// https://github.com/gradientspace/geometry3Sharp/blob/master/mesh/MeshConstraintUtil.cs
void PreserveBoundaryLoops(g3::MeshConstraintsPtr cons, g3::DMesh3Ptr mesh) {
  std::list<g3::EdgeLoopPtr> loops = MeshBoundaryLoops(mesh).Loops;
  for (EdgeLoopPtr loop : loops) {
    DCurve3Ptr loopC = ExtractLoopV(mesh, loop->Vertices);
    IProjectionTargetPtr target =
        std::make_shared<DCurveProjectionTarget>(loopC);
    std::list<int> loop_verts;
    for (int index : loop->Vertices) {
      loop_verts.push_back(index);
    }
    ConstrainVtxLoopToMesh(cons, mesh, loop_verts, target);
  }
}

static void EdgeLengthStats(DMesh3Ptr mesh, double &minEdgeLen,
		double &maxEdgeLen, double &avgEdgeLen,
		int samples = 0) {
	minEdgeLen = std::numeric_limits<double>::max();
	maxEdgeLen = std::numeric_limits<double>::max();
	avgEdgeLen = 0;
	int avg_count = 0;
	int MaxID = mesh->MaxEdgeID();

	// if we are only taking some samples, use a prime-modulo-loop instead of
	// random
	int nPrime = (samples == 0) ? 1 : nPrime = 31337;
	int max_count = (samples == 0) ? MaxID : samples;

	Vector3d a, b;
	int eid = 0;
	int count = 0;
	do {
		if (mesh->IsEdge(eid)) {
			mesh->GetEdgeV(eid, a, b);
			double len = (b - a).norm();
			if (len < minEdgeLen)
				minEdgeLen = len;
			if (len > maxEdgeLen)
				maxEdgeLen = len;
			avgEdgeLen += len;
			avg_count++;
		}
		eid = (eid + nPrime) % MaxID;
	} while (eid != 0 && count++ < max_count);

	avgEdgeLen /= (double)avg_count;
}

//https://github.com/gradientspace/geometry3Sharp/blob/master/mesh/MeshConstraintUtil.cs
void preserve_group_region_border_loops(DMesh3Ptr mesh) {
  int set_id = 1;
  // int[][] group_tri_sets = FaceGroupUtil.FindTriangleSetsByGroup(mesh);
  // for (int[] tri_list : group_tri_sets) {
  //   MeshRegionBoundaryLoops loops = new MeshRegionBoundaryLoops(mesh,
  //   tri_list); foreach (EdgeLoop loop in loops) {
  //     MeshConstraintUtil.ConstrainVtxLoopTo(
  //         r, loop.Vertices, new DCurveProjectionTarget(loop.ToCurve()),
  //         set_id++);
  //   }
  // }
}

Array geometry3_process(Array p_mesh) {
	uint32_t ticks = OS::get_singleton()->get_ticks_msec();
	g3::DMesh3Ptr g3_mesh = std::make_shared<DMesh3>();

	::Vector<::Vector3> vertex_array = p_mesh[Mesh::ARRAY_VERTEX];
	::Vector<::Vector3> normal_array = p_mesh[Mesh::ARRAY_NORMAL];
	::Vector<::Color> color_array = p_mesh[Mesh::ARRAY_COLOR];
	::Vector<::Vector2> uv1_array = p_mesh[Mesh::ARRAY_TEX_UV];
	::Vector<int32_t> bones_array = p_mesh[Mesh::ARRAY_BONES];
	::Vector<float> weights_array = p_mesh[Mesh::ARRAY_WEIGHTS];
	if (normal_array.size()) {
		g3_mesh->EnableVertexNormals(Vector3f(0.0, 0.0, 1.0));
	}
	if (color_array.size()) {
		g3_mesh->EnableVertexColors(Vector3f(1.0, 1.0, 1.0));
	}
	if (uv1_array.size()) {
		g3_mesh->EnableVertexUVs(Vector2f());
	}
	for (int32_t vert_i = 0; vert_i < vertex_array.size(); vert_i++) {
		::Vector3 vert = vertex_array[vert_i];
		::Vector3 normal = ::Vector3(0.0, 0.0, 1.0);
		if (normal_array.size()) {
			normal = normal_array[vert_i];
		}
		::Color color = ::Color(1.0f, 1.0f, 1.0f);
		if (color_array.size()) {
			color = color_array[vert_i];
		}
		::Vector2 uv1;
		if (uv1_array.size()) {
			uv1 = uv1_array[vert_i];
		}
		VectoriDynamic bones;
		VectorfDynamic weights;
		if (bones_array.size() && bones_array.size() == vertex_array.size() * 8
		&&
		    bones_array.size() == weights_array.size()) {
		  const int count = 8;
		  bones.resize(count);
		  for (int32_t i = 0; i < count; i++) {
		    bones[i] = bones_array[vert_i * count + i];
		  }
		  weights.resize(count);
		  for (int32_t i = 0; i < count; i++) {
		    weights[i] = weights_array[vert_i * count + i];
		  }
		} else if (bones_array.size() &&
		           bones_array.size() == vertex_array.size() * 4 &&
		           weights_array.size() == bones_array.size()) {
		  const int count = 4;
		  bones.resize(count);
		  for (int32_t i = 0; i < count; i++) {
		    bones[i] = bones_array[vert_i * count + i];
		  }
		  weights.resize(count);
		  for (int32_t i = 0; i < count; i++) {
		    weights[i] = weights_array[vert_i * count + i];
		  }
		}
		// NewVertexInfo info = NewVertexInfo(Vector3d(vert.x, vert.y, vert.z),
		//                                    Vector3f(normal.x, normal.y,
		//                                    normal.z), Vector3f(color.r, color.g,
		//                                    color.b), Vector2f(uv1.x, uv1.y),
		//                                    bones, weights);
		// TODO Add bone process in geom3 2021-01-23 fire
		// Handle collapse case with multiple bones and every other case.
		NewVertexInfo info = NewVertexInfo(Vector3d(vert.x, vert.y, vert.z),
				Vector3f(normal.x, normal.y, normal.z),
				Vector3f(color.r, color.g, color.b),
				Vector2f(uv1.x, uv1.y));
		g3_mesh->AppendVertex(info);
	}
	::Vector<int32_t> index_array = p_mesh[Mesh::ARRAY_INDEX];
	for (int32_t index_i = 0; index_i < index_array.size(); index_i += 3) {
		ERR_CONTINUE(index_i >= index_array.size());
		Index3i new_tri =
				Index3i(index_array[index_i + 0], index_array[index_i + 1],
						index_array[index_i + 2]);
		g3_mesh->AppendTriangle(new_tri);
	}
	Remesher r(g3_mesh);
	g3::MeshConstraintsPtr cons = std::make_shared<MeshConstraints>();
	// PreserveAllBoundaryEdges(cons, g3_mesh);
	r.SmoothType = Remesher::SmoothTypes::Uniform;
	r.SetExternalConstraints(cons);
	r.SetProjectionTarget(MeshProjectionTarget::AutoPtr(g3_mesh, true));
	PreserveBoundaryLoops(cons, g3_mesh);
	// http://www.gradientspace.com/tutorials/2018/7/5/remeshing-and-constraints
	int iterations = 2;
	r.EnableParallelSmooth = true; // TODO Implement parallel smooth 2021-06-29 FIRE
	r.PreventNormalFlips = true;
	double avg_edge_len = 0.0;
	double min_edge_len = 0.0;
	double max_edge_len = 0.0;
	EdgeLengthStats(g3_mesh, min_edge_len, max_edge_len, avg_edge_len);
	print_line(vformat("avg edge len %.2f", avg_edge_len));
  double target_edge_len = avg_edge_len;
  target_edge_len /= 0.5f;
	print_line(vformat("target edge len %.2f", target_edge_len));
	r.SetTargetEdgeLength(target_edge_len);
	r.SmoothSpeedT = 1.0f / iterations;
	r.Precompute();
	for (int k = 0; k < iterations; ++k) {
		r.BasicRemeshPass();
		print_line("remesh pass " + itos(k));
	}
	print_line("remesh done");
	RemoveFinTriangles(g3_mesh, true);
	std::cout << g3_mesh->MeshInfoString();
	vertex_array.clear();
	index_array.clear();
	uv1_array.clear();
	normal_array.clear();
	color_array.clear();

	for (int32_t vertex_i = 0; vertex_i < g3_mesh->MaxVertexID(); vertex_i++) {
		NewVertexInfo v = g3_mesh->GetVertexAll(vertex_i);

		::Vector3 new_vert;
		Vector3f v_float = v.v.cast<float>();
		new_vert.x = v_float.x();
		new_vert.y = v_float.y();
		new_vert.z = v_float.z();
		vertex_array.push_back(new_vert);

		::Vector3 new_norm;
		Vector3f n_float = v.n.cast<float>();
		new_norm.x = n_float.x();
		new_norm.y = n_float.y();
		new_norm.z = n_float.z();
		normal_array.push_back(new_norm);

		// ::Vector<int> new_bones;
		// // VectoriDynamic bones_float = v.bones;
		// for (int bone_i = 0; bone_i < bones_float.size(); bone_i++) {
		//   new_bones.push_back(bones_float[bone_i]);
		// }
		// bones_array.append_array(new_bones);

		// ::Vector<float> new_weights;
		// VectorfDynamic weights_float = v.weights;
		// for (int weight_i = 0; weight_i < weights_float.size(); weight_i++) {
		//   new_weights.push_back(weights_float[weight_i]);
		// }
		// weights_array.append_array(new_weights);

		::Vector2 new_uv1;
		Vector2f uv_float = v.uv.cast<float>();
		new_uv1.x = uv_float.x();
		new_uv1.y = uv_float.y();
		uv1_array.push_back(new_uv1);
	}
	for (int tri_i = 0; tri_i < g3_mesh->MaxTriangleID(); tri_i++) {
		Index3i tri_index = g3_mesh->GetTriangle(tri_i);
		index_array.push_back(tri_index.x());
		index_array.push_back(tri_index.y());
		index_array.push_back(tri_index.z());
	}

	Array mesh;
	mesh.resize(ArrayMesh::ARRAY_MAX);
	mesh[Mesh::ARRAY_VERTEX] = vertex_array;
	mesh[Mesh::ARRAY_INDEX] = index_array;
	if (uv1_array.size()) {
		mesh[Mesh::ARRAY_TEX_UV] = uv1_array;
	}
	if (normal_array.size()) {
		mesh[Mesh::ARRAY_NORMAL] = normal_array;
	}
	if (color_array.size()) {
		mesh[Mesh::ARRAY_COLOR] = color_array;
	}
	print_line(vformat("remesh took %.2fs",
			(OS::get_singleton()->get_ticks_msec() - ticks) / 1000.0f));
	return mesh;
}
} // namespace g3

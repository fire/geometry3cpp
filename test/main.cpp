#include <iostream>

#include "profile_util.h"
#include <BasicProjectionTargets.h>
#include <DMesh3.h>
#include <DMeshAABBTree3.h>
#include <MeshQueries.h>
#include <MeshSubdivider.h>
#include <OBJReader.h>
#include <OBJWriter.h>
#include <Remesher.h>
#include <VectorUtil.h>
#include <refcount_vector.h>
#include <small_list_set.h>

using namespace g3;

void PreserveAllBoundaryEdges(MeshConstraintsPtr cons, DMesh3* p_mesh)
{
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
				EdgeConstraint::FullyConstrained());

			Index2i ev = p_mesh->GetEdgeV(edge_id);
			VertexConstraint vc = VertexConstraint::Pinned();
			cons->SetOrUpdateVertexConstraint(ev.x(), vc);
			cons->SetOrUpdateVertexConstraint(ev.y(), vc);
		}
	}
}

int main(int argc, char** argv) {
	OBJReader reader;
	DMesh3Builder builder;

	std::ifstream input("sailor.obj");

	BlockTimer read_timer("read", true);
	ReadOptions options = ReadOptions::Defaults();
	options.ReadMaterials = true;
	reader.Read(input, options, builder);
	read_timer.Stop();
	std::cout << "read " << builder.Meshes.size() << " meshes, took "
		<< read_timer.ToString() << std::endl;
	auto mesh1 = builder.Meshes[0];
	std::cout << mesh1->MeshInfoString();
	// TODO 2021-01-21 Color, UV, Groups aren't input // fire

	DMeshAABBTree3 spatialTest(mesh1, true);
	spatialTest.Build();
	spatialTest.TestCoverage();
	BlockTimer remesh_timer("remesh", true);
	Remesher r(mesh1);
	MeshConstraintsPtr cons;
	PreserveAllBoundaryEdges(cons, mesh1.get());
	// https://github.com/gradientspace/geometry3Sharp/blob/master/mesh/MeshConstraintUtil.cs
	//void PreserveBoundaryLoops(MeshConstraints cons, DMesh3 mesh) {
	//  // MeshBoundaryLoops loops = new MeshBoundaryLoops(mesh);
	//  for (int32_t loop_i = 0;;) {
	//    DCurve3 loopC = MeshUtil.ExtractLoopV(mesh, loop.Vertices);
	//    DCurveProjectionTarget target = new DCurveProjectionTarget(loopC);
	//    ConstrainVtxLoopTo(cons, mesh, loop.Vertices, target);
	//  }
	//}
	// https://github.com/gradientspace/geometry3Sharp/blob/master/mesh/MeshConstraintUtil.cs
	// preserve group-region-border-loops
	// int set_id = 1;
	// int[][] group_tri_sets = FaceGroupUtil.FindTriangleSetsByGroup(mesh);
	// foreach (int[] tri_list in group_tri_sets) {
	//     MeshRegionBoundaryLoops loops = new MeshRegionBoundaryLoops(mesh,
	//     tri_list); foreach (EdgeLoop loop in loops) {
	//         MeshConstraintUtil.ConstrainVtxLoopTo(r, loop.Vertices,
	//             new DCurveProjectionTarget(loop.ToCurve()), set_id++);
	//     }
	//  }
	r.SetExternalConstraints(cons);
	r.SetProjectionTarget(MeshProjectionTarget::AutoPtr(mesh1, true));
	// http://www.gradientspace.com/tutorials/2018/7/5/remeshing-and-constraints
	int iterations = 5;
	r.SmoothSpeedT /= iterations;
	r.EnableParallelSmooth = true;
	double avg_edge_len = 0.0;
	for (int32_t edge_i = 1; edge_i < mesh1->EdgeCount(); edge_i++) {
		double edge_len = (mesh1->GetEdgePoint(edge_i - 1, edge_i - 1) -
			mesh1->GetEdgePoint(edge_i - 1, edge_i))
			.norm();
		avg_edge_len += edge_len;
		avg_edge_len /= 2.0;
	}
	std::cout << "avg edge len " << avg_edge_len << std::endl;
	double target_edge_len = avg_edge_len;
	target_edge_len = Clamp(target_edge_len, 0.008, 1.0); // Meters
	std::cout << "target edge len " << target_edge_len << std::endl;
	r.SetTargetEdgeLength(target_edge_len);
	r.Precompute();
	for (int k = 0; k < iterations; ++k) {
		r.BasicRemeshPass();
		std::cout << "remesh pass " << k << std::endl;
	}
	// https://github.com/gradientspace/geometry3Sharp/blob/master/mesh/MeshConstraintUtil.cs
	// RemoveFinTriangles
	// /// <summary>
	// /// Remove 'fin' triangles that have only one connected triangle.
	// /// Removing one fin can create another, by default will keep iterating
	// /// until all fins removed (in a not very efficient way!).
	// /// Pass bRepeatToConvergence=false to only do one pass.
	// /// [TODO] if we are repeating, construct face selection from nbrs of first
	// list and iterate over that on future passes!
	// /// </summary>
	// public static int RemoveFinTriangles(DMesh3 mesh, Func<DMesh3, int, bool>
	// removeF = null, bool bRepeatToConvergence = true)
	// {
	//     MeshEditor editor = new MeshEditor(mesh);
	//     int nRemoved = 0;
	//     List<int> to_remove = new List<int>();
	//     repeat:
	//     foreach ( int tid in mesh.TriangleIndices()) {
	//         Index3i nbrs = mesh.GetTriNeighbourTris(tid);
	//         int c = ((nbrs.a != DMesh3.InvalidID)?1:0) + ((nbrs.b !=
	//         DMesh3.InvalidID)?1:0) + ((nbrs.c != DMesh3.InvalidID)?1:0); if (c
	//         <= 1) {
	//             if (removeF == null || removeF(mesh, tid) == true )
	//                 to_remove.Add(tid);
	//         }
	//     }
	//     if (to_remove.Count == 0)
	//         return nRemoved;
	//     nRemoved += to_remove.Count;
	//     RemoveTriangles(mesh, to_remove, true);
	//     to_remove.Clear();
	//     if (bRepeatToConvergence)
	//         goto repeat;
	//     return nRemoved;
	// }
	remesh_timer.Stop();
	std::cout << "remesh took " << remesh_timer.ToString() << std::endl;
	std::cout << mesh1->MeshInfoString();

	std::ofstream output("output_sailor.obj");
	std::vector<WriteMesh> write_meshes;
	write_meshes.push_back(WriteMesh(mesh1));
	OBJWriter writer;
	writer.Write(output, write_meshes, WriteOptions::Defaults());
}

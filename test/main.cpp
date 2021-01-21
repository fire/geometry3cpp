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
	if (mesh1->EdgeCount() <= 1) {
		return 1;
	}
	double avg_edge_len = 0.0;
	for (int32_t edge_i = 1; edge_i < mesh1->EdgeCount(); edge_i++) {
		double edge_len = (mesh1->GetEdgePoint(edge_i - 1, edge_i - 1) - mesh1->GetEdgePoint(edge_i - 1, edge_i)).norm();
		avg_edge_len += edge_len;
		avg_edge_len /= 2.0;

	}
	std::cout << "avg edge len " << avg_edge_len << std::endl;
	Remesher r(mesh1); 
	MeshConstraintsPtr cons = std::make_shared<MeshConstraints>();
	{
		int32_t max_edge_id = mesh1->MaxEdgeID();
		for (int edge_id = 0; edge_id < max_edge_id; ++edge_id) {
			if (mesh1->IsEdge(edge_id) && mesh1->IsBoundaryEdge(edge_id)) {
				cons->SetOrUpdateEdgeConstraint(edge_id, EdgeConstraint::FullyConstrained());

				Index2i ev = mesh1->GetEdgeV(edge_id);
				VertexConstraint vc = VertexConstraint::Pinned();
				cons->SetOrUpdateVertexConstraint(ev.x(), vc);
				cons->SetOrUpdateVertexConstraint(ev.y(), vc);
			}
		}
	}
	r.SetExternalConstraints(cons);
	r.SetProjectionTarget(MeshProjectionTarget::AutoPtr(mesh1, true));
	// http://www.gradientspace.com/tutorials/2018/7/5/remeshing-and-constraints
	int iterations = 5;
	r.SmoothSpeedT /= iterations;
	r.EnableParallelSmooth = true;
	double target_edge_len = avg_edge_len;
	target_edge_len = Clamp(target_edge_len, 0.008, 1.0); // Meters
	std::cout << "target edge len " << target_edge_len  << std::endl;
	r.SetTargetEdgeLength(target_edge_len);
	r.Precompute();
	for (int k = 0; k < iterations; ++k) {
		r.BasicRemeshPass();
		std::cout << "remesh pass " << k << std::endl;
	}
	remesh_timer.Stop();
	std::cout << "remesh took " << remesh_timer.ToString() << std::endl;
	std::cout << mesh1->MeshInfoString();

	std::ofstream output("output_sailor.obj");
	std::vector<WriteMesh> write_meshes;
	write_meshes.push_back(WriteMesh(mesh1));
	OBJWriter writer;
	writer.Write(output, write_meshes, WriteOptions::Defaults());
}

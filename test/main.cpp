#include <iostream>

#include "profile_util.h"
#include <BasicProjectionTargets.h>
#include <DMesh3.h>
#include <DMeshAABBTree3.h>
#include <MeshQueries.h>
#include <OBJReader.h>
#include <OBJWriter.h>
#include <Remesher.h>
#include <VectorUtil.h>
#include <refcount_vector.h>
#include <small_list_set.h>
#include <MeshSubdivider.h>

using namespace g3;

int main(int argc, char **argv) {
  DMesh3 mesh;
  mesh.AppendVertex(Vector3d(0, 0, 0));
  mesh.AppendVertex(Vector3d(1, 0, 0));
  mesh.AppendVertex(Vector3d(0, 1, 0));
  mesh.AppendVertex(Vector3d(1, 1, 0));
  int t0 = mesh.AppendTriangle(Vector3i(0, 1, 2));
  int t1 = mesh.AppendTriangle(Vector3i(0, 3, 1));

  bool bValid = mesh.CheckValidity();
  mesh.CompactInPlace();
  bValid = mesh.CheckValidity();

  for (Vector3d v : mesh.Vertices()) {
    std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
  }
  for (Index3i t : mesh.Triangles()) {
    std::cout << t[0] << " " << t[1] << " " << t[2] << std::endl;
  }

  for (int eid : mesh.BoundaryEdgeIndices()) {
    Index4i ev = mesh.GetEdge(eid);
    if (ev[3] != InvalidID)
      std::cout << "BAD BOUNDARY EDGE ITR!!";
  }

  for (int vid : mesh.VertexIndices()) {
    std::cout << "vertex " << vid << " tris: ";
    for (int tid : mesh.VtxTrianglesItr(vid))
      std::cout << tid << " ";
    std::cout << std::endl;
  }

  DMeshAABBTree3 spatial(std::make_shared<DMesh3>(mesh), true);
  spatial.Build();
  spatial.TestCoverage();
  double fNearDistSqr;

  int n0 = 0, n1 = 0, fails = 0;
  Wml::Mathd::IntervalRandom(-10, 10, 31337);
  for (int k = 0; k < 1; ++k) {
    double x = Wml::Mathd::IntervalRandom(-10, 10);
    double y = Wml::Mathd::IntervalRandom(-10, 10);
    double z = Wml::Mathd::IntervalRandom(-10, 10);
    int near_tid = spatial.FindNearestTriangle(Vector3d(x, y, z), fNearDistSqr);
    if (near_tid == 0)
      n0++;
    else
      n1++;

    int brute_near_tid =
        MeshQueries::FindNearestTriangle_LinearSearch(mesh, Vector3d(x, y, z));
    if (near_tid != brute_near_tid)
      fails++;
  }
  std::cout << "near0 " << n0 << " near1 " << n1 << " fails " << fails
            << std::endl;

  OBJReader reader;
  DMesh3Builder builder;

  std::ifstream input("sailor.obj");

  BlockTimer read_timer("read", true);
  reader.Read(input, ReadOptions::Defaults(), builder);
  read_timer.Stop();
  std::cout << "read " << builder.Meshes.size() << " meshes, took "
            << read_timer.ToString() << std::endl;
  auto mesh1 = builder.Meshes[0];
  std::cout << mesh1->MeshInfoString();

  double cur_len =
      (mesh1->GetEdgePoint(0, 0) - mesh1->GetEdgePoint(0, 1)).norm();

  DMeshAABBTree3 spatialTest(mesh1, true);
  spatialTest.Build();
  spatialTest.TestCoverage();

  BlockTimer remesh_timer("remesh", true);
  Remesher r(mesh1);
  r.SetProjectionTarget(MeshProjectionTarget::AutoPtr(mesh1, true));
  r.SetTargetEdgeLength(0.01);
  r.SmoothSpeedT = 1.0;
  r.Precompute();
  for (int k = 0; k < 5; ++k) {
  	r.BasicRemeshPass();
	std::cout << "remesh pass " << k << std::endl;
  }

  remesh_timer.Stop();
  std::cout << "remesh took " << remesh_timer.ToString() << std::endl;

  std::ofstream output("output_sailor.obj");
  std::vector<WriteMesh> write_meshes;
  write_meshes.push_back(WriteMesh(mesh1));
  OBJWriter writer;
  writer.Write(output, write_meshes, WriteOptions::Defaults());

  Vector3d axis(0, 0, 1);
  Matrix3d matrix;
  g3::ComputeAlignZAxisMatrix(axis, matrix, false);
}

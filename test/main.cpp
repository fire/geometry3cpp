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

int main(int argc, char **argv) {
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

  double cur_len =
      (mesh1->GetEdgePoint(0, 0) - mesh1->GetEdgePoint(0, 1)).norm();

  DMeshAABBTree3 spatialTest(mesh1, true);
  spatialTest.Build();
  spatialTest.TestCoverage();
  BlockTimer remesh_timer("remesh", true);
  MeshSubdivider s;
  s.Split1to4(*mesh1);
  if (false) {
    Remesher r(mesh1);
    r.SetProjectionTarget(MeshProjectionTarget::AutoPtr(mesh1, true));
    // http://www.gradientspace.com/tutorials/2018/7/5/remeshing-and-constraints
    // TODO 2021-01-21 Fix all boundary edges // fire
    // "half the edge length", approximately
    int iterations = 25;
    r.SmoothSpeedT = 1.0;
    r.SmoothSpeedT /= iterations;
    r.EnableParallelSmooth = true;
    r.PreventNormalFlips = true;
    for (int k = 0; k < iterations; ++k) {
      r.BasicRemeshPass();
      std::cout << "remesh pass " << k << std::endl;
    }
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

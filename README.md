# geometry3cpp

Open-Source (Boost-license) C++11 library for geometric computing. 

geometry3cpp is an in-progress port of [geometry3Sharp](https://github.com/gradientspace/geometry3Sharp), the gradientspace C# library for geometric computing. (Except the parts of that library that were ported from the C++ WildMagic/GTEngine libraries, we just use those directly)

**WORK IN PROGRESS - HARDLY TESTED - BUYER BEWARE**


# Goals

g3cpp is intended to be a general-purpose high-level geometric computing package, with a focus on triangle mesh processing. It includes *other* math libraries which provide most of the low-level vector math stuff, solvers, etc. 

I would like the library to be a header-only library. *However* the WildMagic5 dependency is currently structured to produce a compiled library, and GTEngine is also not fully header-only, there are few .cpp files. My intention is to eventually refactor these libraries so that they are fully header-only as well. 

*Note: I would welcome comments about whether a header-only approach is suitable/desirable for production use.*

# Current State

The dependencies contain an enormous amount of functionality, much more than geometry3Sharp, at the lower level. At the mesh level the following classes have been ported

* **DMesh3** - dynamic mesh, fully ported
* **DMeshAABBTree3** - AABB bounding volume hierarchy for DMesh3, only nearest-point and generic traversal currently ported
* **Remesher** - majority ported, no parallel smoothing/projection currently

**Old Code** This repository has evolved from an earlier attempt (pre-geometry3Sharp) at a mesh processing library. This previous code has not been fully deleted/updated, but will be. 

**goemetry3_tests** project has some basic sample code in **main.cpp**. The CMake generators add this as a separate project, and it generates a standalone executable linking to g3Cpp. If you want to just try something out, I would recommend starting here.


# Dependencies

All dependencies are included in the repository, for convenience.

1) [**Eigen** ](https://eigen.tuxfamily.org/), the C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms (according to their website). **MPL2 License**. Header-only, source included in */external/Eigen* for ease of compiling, but you could probably get CMake to look somewhere else. The Vector types (Vector2d, Vector3d, etc) in g3cpp are Eigen vectors.

2) **WildMagic5** from [GeometricTools](https://www.geometrictools.com/), written by David Eberly. **Boost License**. Only the LibCore and LibMathematics components. Vector math, Geometric intersection and distance tests in 2D and 3D, containment fitters, geometric approximations fitters, Computational Geometry algorithms, Numerical methods, rational number types, 1/2/3D interpolation methods. It's amazing, I've been using WildMagic for 10+ years, since version 2. The source is included in the */external/* subdirectory. This library is no longer maintained and so I have made various local changes to ease porting from the C# version and make it easier to pass vector types between Eigen and Wm5. 

3) [**libigl**](https://libigl.github.io/), the C++ header-only mesh/geometry processing library. libigl is built on Eigen and has implementations of most standard geometry processing algorithms/techniques. Source is included in */external/libigl*. **MPL2 License** is used for the core library, and this is all that geometry3cpp will call. However the code includes calls to various other libaries, including CGAL, LGPL/GPL-licensed code, etc. These will not be included in your binaries unless you explicitly call them via the **igl::copyleft::** namespace. 

# Building

```
mkdir build
cmake ..
cd build
# os specific build like make or geometry3cpp.sln
# execute the tests
```

# libigl interop

Since libigl also uses Eigen, many things are compatible. The main interop required is in passing meshes between the libraries. libigl uses Nx3 Eigen matrices for vertices and triangles (Eigen::MatrixXd and MatrixXi, respectively). [More details in their tutorial](https://libigl.github.io/tutorial/#mesh-representation). g3cpp provides functions to convert to/from DMesh3 as follows:

    Eigen::MatrixXd V; Eigen::MatrixXi F;   // your libigl mesh
    DMesh3 mesh(V,F);                       // convert to DMesh3
    mesh.ToIGLMesh(V,F);                    // convert DMesh3 to libigl mesh. Resizes V and F.

Note that ToIGLMesh() properly handles gaps in the index spaces (ie if vertices or triangles were deleted), so there will not be a 1-1 correspondence between DMesh3 vert/tri indices and libigl row indices unless the DMesh3 is compact.

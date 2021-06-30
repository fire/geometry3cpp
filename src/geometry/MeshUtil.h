#pragma once

#include <DMesh3.h>
#include <GeometryInterfaces.h>
#include <g3types.h>

#include <VectorUtil.h>
#include <Wm5ContBox3.h>
#include <Wm5PolyhedralMassProperties.h>

#include <atomic_util.h>
#include <algorithm>
#include <limits>

namespace g3 {

double Area(const IPackedMesh *pMesh);
double Volume(const IPackedMesh *pMesh);

Box3f AxisBoundingBox(const IPackedMesh *pMesh);
Box3f OrientedBoundingBox(const IPackedMesh *pMesh, bool bFast = true);

void Translate(IPackedMesh *pMesh, const Vector3f &vTranslation);
void Scale(IPackedMesh *pMesh, const Vector3f &vScale);
inline void Scale(IPackedMesh *pMesh, float fScale) {
	Scale(pMesh, Vector3f(fScale, fScale, fScale));
}

inline Vector3d UniformSmooth(DMesh3Ptr mesh, int vID, double t) {
	Vector3d v = mesh->GetVertex(vID);
	Vector3d c;
	mesh->VtxOneRingCentroid(vID, c);
	return (1.0 - t) * v + (t)*c;
}

template <class Real>
const Real *get_value(const Real *pBuffer, unsigned int nIndex,
		unsigned int nStride = 3) {
	return &pBuffer[nIndex * nStride];
}

double Area(const IPackedMesh *pMesh) {
	// was using this in tbb version but we removed tbb...
	// atomic_accumulator<double> area(0.0);

	const float *pVertices = pMesh->GetPositionsBuffer();
	const unsigned int *pTriangles = pMesh->GetIndicesBuffer();
	unsigned int nTriangles = pMesh->GetTriangleCount();

	double area = 0;
	for (unsigned int ti = 0; ti < nTriangles; ++ti) {
		auto pTri = get_value(pTriangles, (unsigned int)ti);
		double a = Area(f2d(Vector3f(get_value(pVertices, pTri[0]))),
				f2d(Vector3f(get_value(pVertices, pTri[1]))),
				f2d(Vector3f(get_value(pVertices, pTri[2]))));
		area += a;
	}

	return area;
}

class PackedMeshVertexSource : public Wml::VertexSource<double> {
public:
	const float *pVertices;
	virtual ~PackedMeshVertexSource() = default;
	virtual Vector3<double> operator[](unsigned int i) const {
		i = 3 * i;
		return Vector3<double>(pVertices[i], pVertices[i + 1],
				pVertices[i + 2]);
	}
};

double Volume(const IPackedMesh *pMesh) {
	PackedMeshVertexSource tmp;
	tmp.pVertices = pMesh->GetPositionsBuffer();

	double fMass;
	Vector3d vCoM;
	Wml::Matrix3d vMoI;
	Wml::ComputeMassProperties(&tmp, (int)pMesh->GetTriangleCount(),
			(const int *)pMesh->GetIndicesBuffer(), false,
			fMass, vCoM, vMoI);
	return fMass;
}

Box3f AxisBoundingBox(const IPackedMesh *pMesh) {
	int NV = pMesh->GetVertexCount();
	std::vector<Vector3f> pts;
	pts.reserve(NV);
	const float *p = pMesh->GetPositionsBuffer();
	for (int k = 0; k < NV; ++k)
		pts.push_back(Vector3f(p[3 * k], p[3 * k + 1], p[3 * k + 2]));
	return Wml::ContAlignedBox(NV, &pts[0]);
}

Box3f OrientedBoundingBox(const IPackedMesh *pMesh, bool bFast) {
	int NV = pMesh->GetVertexCount();
	std::vector<Vector3f> pts;
	pts.reserve(NV);
	const float *p = pMesh->GetPositionsBuffer();
	for (int k = 0; k < NV; ++k)
		pts.push_back(Vector3f(p[3 * k], p[3 * k + 1], p[3 * k + 2]));
	return Wml::ContOrientedBox(NV, &pts[0]);
}

void Translate(IPackedMesh *pMesh, const Vector3f &vTranslation) {
	// AAAAHHHH
	float *pVertices = const_cast<float *>(pMesh->GetPositionsBuffer());
	unsigned int nVertices = pMesh->GetVertexCount();

	for (unsigned int vi = 0; vi < nVertices; ++vi) {
		unsigned int k = vi * 3;
		pVertices[k] += vTranslation[0];
		pVertices[k + 1] += vTranslation[1];
		pVertices[k + 2] += vTranslation[2];
	}

	pMesh->updateTimeStamp();
}

void Scale(IPackedMesh *pMesh, const Vector3f &vScale) {
	// AAAAHHHH
	float *pVertices = const_cast<float *>(pMesh->GetPositionsBuffer());
	unsigned int nVertices = pMesh->GetVertexCount();

	for (unsigned int vi = 0; vi < nVertices; ++vi) {
		unsigned int k = vi * 3;
		pVertices[k] *= vScale[0];
		pVertices[k + 1] *= vScale[1];
		pVertices[k + 2] *= vScale[2];
	}

	pMesh->updateTimeStamp();
}

static g3::DCurve3Ptr ExtractLoopV(g3::DMesh3Ptr mesh, std::vector<int> vertices) {
	g3::DCurve3Ptr curve = std::make_shared<g3::DCurve3>();
	for (int vid : vertices) {
		curve->AppendVertex(mesh->GetVertex(vid));
	}
	curve->SetClosed(true);
	return curve;
}
};
} // namespace g3

#pragma once

#include <GeometryInterfaces.h>
#include <g3types.h>
#include "VectorUtil.h"
#include "g3Debug.h"

namespace g3 {

/*
	 * Simple mesh class where data is stored in internal buffers.
	 * Colors are float-RGB (enforced internally w/ asserts, etc)
	 */
class g3External PackedMesh : public IPackedMesh {
public:
	virtual ~PackedMesh() {}

	/*
		* IPackedMesh interface
		*/

	virtual unsigned int GetVertexCount() const override;
	virtual unsigned int GetTriangleCount() const override;

	virtual bool HasPositions() const override;
	virtual const float *GetPositionsBuffer() const override;

	virtual bool HasNormals() const override;
	virtual const float *GetNormalsBuffer() const override;

	virtual bool HasColorsFloat() const override;
	virtual const float *GetColorsFloatBuffer() const override;

	virtual bool HasIndices() const override;
	virtual const unsigned int *GetIndicesBuffer() const override;

	virtual bool HasTriColors() const override;
	virtual const float *GetTriColorsBuffer() const override;

	virtual bool HasTriGroups() const override;
	virtual const unsigned int *GetTriGroupsBuffer() const override;

	/*
		 * construction
		 */
	void ResizeVertices(unsigned int nVertices, bool bNormals, bool bColors);
	void ResizeTriangles(unsigned int nTriangles);

	void CopyVertices(const std::vector<float> &vPositions,
			const std::vector<float> *pNormals = nullptr,
			const std::vector<float> *pColors = nullptr);
	void CopyVertices(const float *pPositions, size_t nVertices,
			const float *pNormals = nullptr,
			const float *pColors = nullptr);

	void CopyTriangles(const std::vector<unsigned int> &vIndices,
			const std::vector<unsigned int> *pGroups = nullptr);
	void CopyTriangles(const unsigned int *pIndices, size_t nTriangles,
			const unsigned int *pGroups = nullptr);
	void CopyTriangles(const int *pIndices, size_t nTriangles,
			const int *pGroups = nullptr);
	void CopyTriColors(const std::vector<float> *pColors);
	void CopyTriColors(const float *pColors, size_t nTriangles);
	void CopyTriGroups(const std::vector<unsigned int> *pGroups);
	void CopyTriGroups(const unsigned int *pGroups, size_t nTriangles);

	void SetConstantTriGroup(unsigned int nGroupID);

	/*
		 * Utility
		 */
	float *GetPosition(unsigned int k) {
		return &m_vPositions[3 * k];
	}
	const float *GetPosition(unsigned int k) const {
		return &m_vPositions[3 * k];
	}
	void SetPosition(unsigned int k, float x, float y, float z) {
		float *f = &m_vPositions[3 * k];
		f[0] = x;
		f[1] = y;
		f[2] = z;
	}

	float *GetNormal(unsigned int k) {
		return &m_vNormals[3 * k];
	}
	const float *GetNormal(unsigned int k) const {
		return &m_vNormals[3 * k];
	}
	float *GetColor(unsigned int k) {
		return &m_vColors[3 * k];
	}
	const float *GetColor(unsigned int k) const {
		return &m_vColors[3 * k];
	}
	float *GetTriColor(unsigned int k) {
		return &m_vTriColors[3 * k];
	}
	const float *GetTriColor(unsigned int k) const {
		return &m_vTriColors[3 * k];
	}
	unsigned int *GetTriGroup(unsigned int k) {
		return &m_vTriGroups[k];
	}
	unsigned int GetTriGroup(unsigned int k) const {
		return m_vTriGroups[k];
	}

	void GetTriangle(unsigned int k, unsigned int &i0, unsigned int &i1, unsigned int &i2) const;
	void GetTriangle(unsigned int k, unsigned int *iTri) const;
	void GetTriangle(unsigned int k, Vector3 &v0, Vector3 &v1, Vector3 &v2) const;
	void GetTriangle(unsigned int k, Vector3 *pTri) const;

	void SetTriangle(unsigned int k, int a, int b, int c) {
		unsigned int *t = &m_vIndices[3 * k];
		t[0] = a;
		t[1] = b;
		t[2] = c;
	}

	Vector3 GetTriangleNormal(unsigned int k) const;

	void EstimateNormals();

protected:
	std::vector<float> m_vPositions;
	std::vector<float> m_vNormals;
	std::vector<float> m_vColors;
	std::vector<unsigned int> m_vIndices;
	std::vector<float> m_vTriColors;
	std::vector<unsigned int> m_vTriGroups;
};

unsigned int PackedMesh::GetVertexCount() const {
	return (unsigned int)m_vPositions.size() / 3;
}

unsigned int PackedMesh::GetTriangleCount() const {
	return (unsigned int)m_vIndices.size() / 3;
}

bool PackedMesh::HasPositions() const {
	return !m_vPositions.empty();
}

const float *PackedMesh::GetPositionsBuffer() const {
	return HasPositions() ? &m_vPositions[0] : nullptr;
}

bool PackedMesh::HasNormals() const {
	return HasPositions() && m_vNormals.size() == m_vPositions.size();
}

const float *PackedMesh::GetNormalsBuffer() const {
	return HasNormals() ? &m_vNormals[0] : nullptr;
}

bool PackedMesh::HasColorsFloat() const {
	return HasPositions() && m_vPositions.size() == m_vColors.size();
}

const float *PackedMesh::GetColorsFloatBuffer() const {
	return HasColorsFloat() ? &m_vColors[0] : nullptr;
}

bool PackedMesh::HasIndices() const {
	return !m_vIndices.empty();
}

const unsigned int *PackedMesh::GetIndicesBuffer() const {
	return HasIndices() ? &m_vIndices[0] : nullptr;
}

bool PackedMesh::HasTriColors() const {
	return !m_vTriColors.empty();
}
const float *PackedMesh::GetTriColorsBuffer() const {
	return HasTriColors() ? &m_vTriColors[0] : nullptr;
}

bool PackedMesh::HasTriGroups() const {
	return !m_vTriGroups.empty();
}
const unsigned int *PackedMesh::GetTriGroupsBuffer() const {
	return HasTriGroups() ? &m_vTriGroups[0] : nullptr;
}

void PackedMesh::ResizeVertices(unsigned int nVertices, bool bNormals,
		bool bColors) {
	m_vPositions.resize(nVertices * 3);
	if (bNormals)
		m_vNormals.resize(m_vPositions.size());
	if (bColors)
		m_vColors.resize(m_vColors.size());
	updateTimeStamp();
}

void PackedMesh::ResizeTriangles(unsigned int nTriangles) {
	m_vIndices.resize(nTriangles * 3);
	updateTimeStamp();
}

void PackedMesh::CopyVertices(const std::vector<float> &vPositions,
		const std::vector<float> *pNormals,
		const std::vector<float> *pColors) {
	m_vPositions = vPositions;
	if (pNormals) {
		gDevAssert(pNormals->size() == m_vPositions.size());
		m_vNormals = *pNormals;
	}
	if (pColors) {
		gDevAssert(pColors->size() == m_vPositions.size());
		m_vColors = *pColors;
	}
	updateTimeStamp();
}

void PackedMesh::CopyVertices(const float *pPositions, size_t nVertices,
		const float *pNormals, const float *pColors) {
	m_vPositions.resize(nVertices * 3);
	memcpy(&m_vPositions[0], pPositions, nVertices * 3 * sizeof(float));

	if (pNormals) {
		m_vNormals.resize(nVertices * 3);
		memcpy(&m_vNormals[0], pNormals, nVertices * 3 * sizeof(float));
	}
	if (pColors) {
		m_vColors.resize(nVertices * 3);
		memcpy(&m_vColors[0], pColors, nVertices * 3 * sizeof(float));
	}
	updateTimeStamp();
}

void PackedMesh::CopyTriangles(const std::vector<unsigned int> &vIndices,
		const std::vector<unsigned int> *pGroups) {
	m_vIndices = vIndices;
	if (pGroups)
		m_vTriGroups = *pGroups;
	updateTimeStamp();
}

void PackedMesh::CopyTriangles(const unsigned int *pIndices,
		size_t nTriangles,
		const unsigned int *pGroups) {
	m_vIndices.resize(nTriangles * 3);
	memcpy(&m_vIndices[0], pIndices, nTriangles * 3 * sizeof(unsigned int));
	if (pGroups) {
		m_vTriGroups.resize(nTriangles);
		memcpy(&m_vTriGroups[0], pGroups, nTriangles * sizeof(unsigned int));
	}
	updateTimeStamp();
}

void PackedMesh::CopyTriangles(const int *pIndices, size_t nTriangles,
		const int *pGroups) {
	CopyTriangles((unsigned int *)pIndices, nTriangles, (unsigned int *)pGroups);
	updateTimeStamp();
}

void PackedMesh::CopyTriColors(const std::vector<float> *pColors) {
	gDevAssert(pColors->size() == GetTriangleCount() * 3);
	m_vTriColors = *pColors;
	updateTimeStamp();
}
void PackedMesh::CopyTriColors(const float *pColors, size_t nTriangles) {
	gDevAssert(nTriangles == GetTriangleCount());
	m_vTriColors.resize(nTriangles * 3);
	memcpy(&m_vTriColors[0], pColors, nTriangles * 3 * sizeof(float));
	updateTimeStamp();
}

void PackedMesh::CopyTriGroups(const std::vector<unsigned int> *pGroups) {
	gDevAssert(pGroups->size() == GetTriangleCount());
	m_vTriGroups = *pGroups;
	updateTimeStamp();
}
void PackedMesh::CopyTriGroups(const unsigned int *pGroups, size_t nTriangles) {
	gDevAssert(nTriangles == GetTriangleCount());
	m_vTriGroups.resize(nTriangles);
	memcpy(&m_vTriGroups[0], pGroups, nTriangles * sizeof(unsigned int));
	updateTimeStamp();
}

void PackedMesh::SetConstantTriGroup(unsigned int nGroupID) {
	m_vTriGroups.resize(GetTriangleCount());
	std::fill_n(m_vTriGroups.begin(), GetTriangleCount(), nGroupID);
	updateTimeStamp();
}

void PackedMesh::GetTriangle(unsigned int k, unsigned int &i0, unsigned int &i1,
		unsigned int &i2) const {
	const unsigned int *t = &m_vIndices[3 * k];
	i0 = t[0];
	i1 = t[1];
	i2 = t[2];
}
void PackedMesh::GetTriangle(unsigned int k, unsigned int *iTri) const {
	const unsigned int *t = &m_vIndices[3 * k];
	iTri[0] = t[0];
	iTri[1] = t[1];
	iTri[2] = t[2];
}

void PackedMesh::GetTriangle(unsigned int k, Vector3f &v0, Vector3f &v1,
		Vector3f &v2) const {
	unsigned int i0, i1, i2;
	GetTriangle(k, i0, i1, i2);
	v0 = Vector3f(GetPosition(i0));
	v1 = Vector3f(GetPosition(i1));
	v2 = Vector3f(GetPosition(i2));
}
void PackedMesh::GetTriangle(unsigned int k, Vector3f *pTri) const {
	unsigned int i0, i1, i2;
	GetTriangle(k, i0, i1, i2);
	pTri[0] = Vector3f(GetPosition(i0));
	pTri[1] = Vector3f(GetPosition(i1));
	pTri[2] = Vector3f(GetPosition(i2));
}

Vector3f PackedMesh::GetTriangleNormal(unsigned int k) const {
	unsigned int i0, i1, i2;
	GetTriangle(k, i0, i1, i2);
	return Normal(Vector3f(GetPosition(i0)),
			Vector3f(GetPosition(i1)),
			Vector3f(GetPosition(i2)));
}

void PackedMesh::EstimateNormals() {
	m_vNormals.resize(0);
	m_vNormals.resize(m_vPositions.size(), 0.0f);
	float *pNormals = &m_vNormals[0];

	unsigned int nTriangles = GetTriangleCount();
	std::vector<Vector3f> vTriNormals(nTriangles);
	for (unsigned int k = 0; k < nTriangles; ++k)
		vTriNormals[k] = GetTriangleNormal(k);

	for (unsigned int k = 0; k < nTriangles; ++k) {
		unsigned int nTri[3];
		GetTriangle(k, nTri);
		array3f_add(pNormals, nTri[0], vTriNormals[k].data());
		array3f_add(pNormals, nTri[1], vTriNormals[k].data());
		array3f_add(pNormals, nTri[2], vTriNormals[k].data());
	}

	unsigned int nVertices = GetVertexCount();
	for (unsigned int k = 0; k < nVertices; ++k) {
		array3f_normalize(pNormals, k);
	}

	updateTimeStamp();
}

} // namespace g3

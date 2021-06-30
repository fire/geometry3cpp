#pragma once

#include "g3Debug.h"
#include <GeometryInterfaces.h>
#include <g3types.h>
#include <string.h>

namespace g3 {

class PackedLines : public IPackedLines {
public:
	virtual ~PackedLines() {}

	// /*
	// 	* IPackedLines interface
	// 	*/
	// virtual LinesType GetLinesType() const;

	// virtual unsigned int GetVertexCount() const;
	// virtual unsigned int GetLineCount() const;

	// virtual bool HasPositions() const;
	// virtual const float *GetPositionsBuffer() const;

	// virtual bool HasColorsFloat() const;
	// virtual const float *GetColorsFloatBuffer() const;

	// virtual bool HasIndices() const;
	// virtual const unsigned int *GetIndicesBuffer() const;

	// /*
	// 	 * construction
	// 	 */
	// void ResizeVertices(unsigned int nVertices, bool bNormals, bool bColors);
	// void ResizeLines(unsigned int nTriangles);

	// void CopyVertices(const std::vector<float> &vPositions,
	// 		const std::vector<float> *pColors = nullptr);
	// void CopyVertices(const float *pPositions, size_t nVertices,
	// 		const float *pColors = nullptr);

	// void CopyLines(const std::vector<unsigned int> &vIndices);
	// void CopyLines(const unsigned int *pIndices, size_t nLines);
	// void CopyLines(const int *pIndices, size_t nLines);

	// /*
	// 	 * Utility
	// 	 */
	// const float *GetPosition(unsigned int k) const {
	// 	return &m_vPositions[3 * k];
	// }
	// const float *GetColor(unsigned int k) const {
	// 	return &m_vColors[3 * k];
	// }
	// void GetSegmentIndices(unsigned int k, unsigned int &i0, unsigned int &i1) const;
	// void GetSegment(unsigned int k, Vector3 &v0, Vector3 &v1) const;

	PackedLines(LinesType eType) { m_eType = eType; }

	LinesType GetLinesType() const { return m_eType; }

	unsigned int GetVertexCount() const {
		return (unsigned int)m_vPositions.size() / 3;
	}

	unsigned int GetLineCount() const {
		switch (m_eType) {
			default:
			case Segments:
				return (unsigned int)m_vPositions.size() / 2;
			case Strip:
				return (m_vPositions.empty()) ? 0 : (unsigned int)m_vPositions.size() - 1;
			case Loop:
				return (unsigned int)m_vPositions.size();
			case IndexedSegments:
				return (unsigned int)m_vIndices.size() / 2;
		}
	}

	bool HasPositions() const { return !m_vPositions.empty(); }
	const float *GetPositionsBuffer() const {
		return HasPositions() ? &m_vPositions[0] : nullptr;
	}

	bool HasColorsFloat() const {
		return HasPositions() && m_vPositions.size() == m_vColors.size();
	}
	const float *GetColorsFloatBuffer() const {
		return HasColorsFloat() ? &m_vColors[0] : nullptr;
	}

	bool HasIndices() const {
		return (m_eType == IndexedSegments) && (!m_vIndices.empty());
	}
	const unsigned int *GetIndicesBuffer() const {
		return HasIndices() ? &m_vIndices[0] : nullptr;
	}

	void ResizeVertices(unsigned int nVertices, bool bNormals,
			bool bColors) {
		m_vPositions.resize(nVertices * 3);
		if (bColors)
			m_vColors.resize(m_vColors.size());
		updateTimeStamp();
	}

	void ResizeLines(unsigned int nLines) {
		m_vIndices.resize(nLines * 2);
		updateTimeStamp();
	}

	void CopyVertices(const std::vector<float> &vPositions,
			const std::vector<float> *pColors) {
		m_vPositions = vPositions;
		if (pColors) {
			gDevAssert(pColors->size() == m_vPositions.size());
			m_vColors = *pColors;
		}
		updateTimeStamp();
	}

	void CopyVertices(const float *pPositions, size_t nVertices,
			const float *pColors) {
		m_vPositions.resize(nVertices * 3);
		memcpy(&m_vPositions[0], pPositions, nVertices * 3 * sizeof(float));

		if (pColors) {
			m_vColors.resize(nVertices * 3);
			memcpy(&m_vColors[0], pColors, nVertices * 3 * sizeof(float));
		}
		updateTimeStamp();
	}

	void CopyLines(const std::vector<unsigned int> &vIndices) {
		gDevAssert(m_eType == Segments || m_eType == IndexedSegments);
		m_vIndices = vIndices;
		m_eType = IndexedSegments;
		updateTimeStamp();
	}

	void CopyLines(const unsigned int *pIndices, size_t nLines) {
		gDevAssert(m_eType == Segments || m_eType == IndexedSegments);
		m_vIndices.resize(nLines * 2);
		memcpy(&m_vIndices[0], pIndices, nLines * 2 * sizeof(unsigned int));
		m_eType = IndexedSegments;
		updateTimeStamp();
	}

	void CopyLines(const int *pIndices, size_t nLines) {
		CopyLines((unsigned int *)pIndices, nLines);
		updateTimeStamp();
	}

	void GetSegmentIndices(unsigned int k, unsigned int &i0,
			unsigned int &i1) const {
		i0 = k;
		i1 = (k + 1) % m_vPositions.size();
		if (m_eType == Segments) {
			i0 = 2 * k;
			i1 = 2 * k + 1;
		} else if (m_eType == IndexedSegments) {
			i0 = m_vIndices[2 * k];
			i1 = m_vIndices[2 * k + 1];
		}
	}

	void GetSegment(unsigned int k, Vector3f &v0,
			Vector3f &v1) const {
		unsigned int i0, i1;
		GetSegmentIndices(k, i0, i1);
		v0 = Vector3f(&m_vPositions[i0]);
		v1 = Vector3f(&m_vPositions[i1]);
	}

protected:
	LinesType m_eType;

	std::vector<float> m_vPositions;
	std::vector<float> m_vColors;

	std::vector<unsigned int> m_vIndices;
};

} // namespace g3
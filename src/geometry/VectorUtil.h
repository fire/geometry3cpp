#pragma once

#include <g3types.h>
#include <math.h>

/*
 * VectorUtil contains a set of utility functions for matrix/vector math.
 */

namespace g3 {

/*
	 * if bToZAxis is false, compute matrix that rotates Z axis into vAlignWith
	 * if bToZAxis is true, compute matrix that rotates vAlignWith into Z axis
	 */
// template <class real_t>
// void ComputeAlignZAxisMatrix(const Vector3 &vAlignWith,
// 		Basis &matrix, bool bToZAxis = false);

// template <class real_t>
// void ComputeAlignAxisMatrix(const Vector3 &vInitial,
// 		const Vector3 &vAlignWith, Basis &matrix);

// //! compute vectors in a plane perpendicular to vIn
// template <class real_t>
// void ComputePerpVectors(const Vector3 &vIn,
// 		Vector3 &vOut1, Vector3 &vOut2,
// 		bool bInIsNormalized = false);

// //! compute tangent vectors in plane perp to vNormal, using non-orthogonal vEstX as estimate of vOut1
// template <class real_t>
// void ComputePerpVectors(const Vector3 &vNormal, const Vector3 &vEstX,
// 		Vector3 &vOut1, Vector3 &vOut2,
// 		bool bInputIsNormalized = false);

// template <class real_t>
// void ToGLMatrix(const Basis &matrix, real_t glMatrix[16]);

// template <class real_t>
// Vector2<real_t> ToUV(const Vector3 &vec, int nUIndex, int nVIndex);
// template <class real_t>
// Vector3 To3D(const Vector2<real_t> &vec, int nUIndex, int nVIndex);

// template <class real_t>
// real_t VectorAngleR(const Vector2<real_t> &v1, const Vector2<real_t> &v2);
// template <class real_t>
// real_t VectorAngleR(const Vector3 &v1, const Vector3 &v2);

// template <class real_t>
// real_t VectorAngleD(const Vector2<real_t> &v1, const Vector2<real_t> &v2);
// template <class real_t>
// real_t VectorAngleD(const Vector3 &v1, const Vector3 &v2);

// template <class real_t>
// real_t VectorCot(const Vector3 &v1, const Vector3 &v2);

// template <class real_t>
// Vector2<real_t> Lerp(const Vector2<real_t> &v1, const Vector2<real_t> &v2, real_t t);
// template <class real_t>
// Vector3 Lerp(const Vector3 &v1, const Vector3 &v2, real_t t);

// template <class real_t>
// void BarycentricCoords(const Vector3 &vTriVtx1,
// 		const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3,
// 		const Vector3 &vVertex,
// 		real_t &fBary1, real_t &fBary2, real_t &fBary3);

// template <class real_t>
// real_t Area(const Vector3 &vTriVtx1,
// 		const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3);

// template <class real_t>
// void BarycentricCoords(const Vector2<real_t> &vTriVtx1,
// 		const Vector2<real_t> &vTriVtx2,
// 		const Vector2<real_t> &vTriVtx3,
// 		const Vector2<real_t> &vVertex,
// 		real_t &fBary1, real_t &fBary2, real_t &fBary3);

// template <class real_t>
// real_t Area(const Vector2<real_t> &vTriVtx1,
// 		const Vector2<real_t> &vTriVtx2,
// 		const Vector2<real_t> &vTriVtx3);

// template <class real_t>
// Vector3 Normal(const Vector3 &vTriVtx1,
// 		const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3, real_t *pArea = nullptr);

// template <class real_t>
// Vector3 InterpNormal(const Vector3 &vTriVtx1,
// 		const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3,
// 		const Vector3 &vTriNorm1,
// 		const Vector3 &vTriNorm2,
// 		const Vector3 &vTriNorm3,
// 		const Vector3 &vPointInTri);

// //! This metric is from Texture Mapping Progressive Meshes, Sander et al, Siggraph 2001
// template <class real_t>
// void StretchMetric1(const Vector3 &vTriVtx1,
// 		const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3,
// 		const Vector2<real_t> &vVtxParam1,
// 		const Vector2<real_t> &vVtxParam2,
// 		const Vector2<real_t> &vVtxParam3,
// 		real_t &MaxSV, real_t &MinSV, real_t &L2Norm, real_t &LInfNorm);

// template <class real_t>
// void StretchMetric3(const Vector3 &vTriVtx1,
// 		const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3,
// 		const Vector3 &vVtxParam1,
// 		const Vector3 &vVtxParam2,
// 		const Vector3 &vVtxParam3,
// 		real_t &MaxSV, real_t &MinSV, real_t &L2Norm, real_t &LInfNorm);

// template <class real_t>
// bool IsObtuse(const Vector2<real_t> &v1, const Vector2<real_t> &v2, const Vector2<real_t> &v3);
// template <class real_t>
// bool IsObtuse(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3);

// /*
// 	 * inline utilities
// 	 */

// inline bool IsFinite(const Vector3d &v) {
// 	return isfinite(v.x()) && isfinite(v.y()) && isfinite(v.z());
// }

// inline Vector2f d2f(const Vector2d &v) {
// 	return Vector2f((float)v[0], (float)v[1]);
// }
// inline Vector3f d2f(const Vector3d &v) {
// 	return Vector3f((float)v[0], (float)v[1], (float)v[2]);
// }
// inline Vector2d f2d(const Vector2f &v) {
// 	return Vector2d((double)v[0], (double)v[1]);
// }
// inline Vector3d f2d(const Vector3f &v) {
// 	return Vector3d((double)v[0], (double)v[1], (double)v[2]);
// }

// inline Matrix2f d2f(const Matrix2d &v) {
// 	return v.cast<float>();
// }
// inline Matrix3f d2f(const Matrix3d &v) {
// 	return v.cast<float>();
// }
// inline Matrix3d f2d(const Matrix3f &v) {
// 	return v.cast<double>();
// }
// inline Matrix2d f2d(const Matrix2f &v) {
// 	return v.cast<double>();
// }

//inline AxisAlignedBox2f d2f(const AxisAlignedBox2d & v) {
//	return AxisAlignedBox2f((float)v.Min[0], (float)v.Max[0], (float)v.Min[1], (float)v.Max[1] );
//}
//inline AxisAlignedBox3f d2f(const AABB & v) {
//	return AxisAlignedBox3f((float)v.Min[0], (float)v.Max[0], (float)v.Min[1], (float)v.Max[1], (float)v.Min[2], (float)v.Max[2] );
//}
//inline AxisAlignedBox2d f2d(const AxisAlignedBox2f & v) {
//	return AxisAlignedBox2d((double)v.Min[0], (double)v.Max[0], (double)v.Min[1], (double)v.Max[1]);
//}
//inline AABB f2d(const AxisAlignedBox3f & v) {
//	return AABB((double)v.Min[0], (double)v.Max[0], (double)v.Min[1], (double)v.Max[1], (double)v.Min[2], (double)v.Max[2]);
//}

// template <class real_t>
// inline real_t Clamp(const real_t &fValue, const real_t &fMin, const real_t &fMax) {
// 	if (fValue < fMin)
// 		return fMin;
// 	else if (fValue > fMax)
// 		return fMax;
// 	else
// 		return fValue;
// }

// inline void array3f_add(float *pBuffer, unsigned int nIndex, const float *pAdd) {
// 	pBuffer[3 * nIndex] += pAdd[0];
// 	pBuffer[3 * nIndex + 1] += pAdd[1];
// 	pBuffer[3 * nIndex + 2] += pAdd[2];
// }
// inline void array3f_normalize(float *pBuffer, unsigned int nIndex, float fEpsilon = 0.0f) {
// 	auto v = Vector3f(&pBuffer[3 * nIndex]).normalized();
// 	pBuffer[3 * nIndex] = v[0];
// 	pBuffer[3 * nIndex + 1] = v[1];
// 	pBuffer[3 * nIndex + 2] = v[2];
// }
// inline void vectorf_push(std::vector<float> &v, const Vector3f &p) {
// 	v.push_back(p[0]);
// 	v.push_back(p[1]);
// 	v.push_back(p[2]);
// }
// inline void vectori_push(std::vector<unsigned int> &v, const Vector2i &p) {
// 	v.push_back(p[0]);
// 	v.push_back(p[1]);
// }
// inline void vectori_push(std::vector<unsigned int> &v, const Vector3i &p) {
// 	v.push_back(p[0]);
// 	v.push_back(p[1]);
// 	v.push_back(p[2]);
// }

template <typename DerivedA, typename DerivedB>
inline bool EpsilonEqual(const Eigen::MatrixBase<DerivedA> &m1, const Eigen::MatrixBase<DerivedB> &m2, double eps) {
	return (m1 - m2).cwiseAbs().maxCoeff() <= eps;
}

/// <summary>
/// compute vector in direction of triangle normal (cross-product). No normalization.
/// </summary>
/// <returns>The normal direction.</returns>
template <typename DerivedA, typename DerivedB, typename DerivedC>
inline Eigen::MatrixBase<DerivedC> FastNormalDirection(
		const Eigen::MatrixBase<DerivedA> &v1,
		const Eigen::MatrixBase<DerivedB> &v2,
		const Eigen::MatrixBase<DerivedC> &v3) {
	return (v2 - v1).cross(v3 - v1);
}

template <typename DerivedA, typename T>
inline bool Contains(const Eigen::MatrixBase<DerivedA> &m1, T value) {
	int n = (int)m1.array().size();
	for (int k = 0; k < n; ++k) {
		if (m1.array()[k] == value)
			return true;
	}
	return false;
}

template <typename T1, typename T2>
inline bool ContainsKey(const std::map<T1, T2> &dict, T1 value) {
	return dict.find(value) != dict.end();
}

template <typename T>
inline bool Contains(const std::vector<T> &vec, T value) {
	return std::find(vec.begin(), vec.end(), value) != vec.end();
}
template <typename T>
inline bool Remove(std::vector<T> &vec, T value) {
	auto itr = std::find(vec.begin(), vec.end(), value);
	if (itr == vec.end())
		return false;
	vec.erase(itr);
	return true;
}

// template <class real_t>
// void ComputeAlignZAxisMatrix(const Vector3 &vAlignWith,
// 		Basis &matrix, bool bInvert) {
// 	// compute cosine of angle between vectors
// 	real_t axisDot = vAlignWith.dot(Vector3::UnitZ());

// 	// compute rotation axis
// 	Vector3 axisCross(Vector3::UnitZ().cross(vAlignWith));

// 	real_t fInverter = (bInvert) ? (real_t)-1 : (real_t)1;

// 	// apply rotation if necessary
// 	if (axisCross.squaredNorm() > Math<real_t>::EPSILON) {
// 		// compute normalized axis and angle, then create rotation around axis
// 		axisCross.normalize();
// 		real_t fAngle = Math<real_t>::ACos(axisDot / vAlignWith.norm());
// 		matrix = Wml::Basis(axisCross, fAngle * fInverter);

// 	} else if (axisDot < (real_t)0) {
// 		axisCross = Vector3::UnitX();
// 		real_t fAngle = (real_t)180 * Math<real_t>::DEG_TO_RAD * fInverter;
// 		matrix = Wml::Basis(axisCross, fAngle);
// 	} else {
// 		matrix = Basis::Identity();
// 	}
// }

// template <class real_t>
// void ComputeAlignAxisMatrix(const Vector3 &vInitial,
// 		const Vector3 &vAlignWith,
// 		Basis &matrix) {
// 	// compute cosine of angle between vectors
// 	real_t axisDot = vAlignWith.dot(vInitial);

// 	// compute rotation axis
// 	Vector3 axisCross(vInitial.cross(vAlignWith));

// 	// apply rotation if necessary
// 	if (axisCross.squaredNorm() > Math<real_t>::EPSILON) {
// 		// compute normalized axis and angle, then create rotation around axis
// 		axisCross.normalize();
// 		real_t fAngle = Math<real_t>::ACos(axisDot / vAlignWith.norm());
// 		matrix = Wml::Basis(axisCross, fAngle);

// 	} else if (axisDot < (real_t)0) {
// 		// find some perpendicular vectors
// 		Vector3 vPerp1, vPerp2;
// 		ComputePerpVectors(vInitial, vPerp1, vPerp2);

// 		matrix = Wml::Basis(vPerp1, (real_t)180 * Math<real_t>::DEG_TO_RAD);
// 	} else {
// 		matrix = Basis::Identity();
// 	}
// }

// template <class real_t>
// void ComputePerpVectors(const Vector3 &vIn, Vector3 &vOut1,
// 		Vector3 &vOut2, bool bInIsNormalized) {
// 	Vector3 vPerp(vIn);
// 	if (!bInIsNormalized)
// 		vPerp.normalize();

// 	if (Math<real_t>::FAbs(vPerp.x()) >= Math<real_t>::FAbs(vPerp.y()) &&
// 			Math<real_t>::FAbs(vPerp.x()) >= Math<real_t>::FAbs(vPerp.z())) {
// 		vOut1.x() = -vPerp.y();
// 		vOut1.y() = vPerp.x();
// 		vOut1.z() = (real_t)0.0;
// 	} else {
// 		vOut1.x() = (real_t)0.0;
// 		vOut1.y() = vPerp.z();
// 		vOut1.z() = -vPerp.y();
// 	}

// 	vOut1.normalize();
// 	vOut2 = vPerp.cross(vOut1);
// }

// template <class real_t>
// void ComputePerpVectors(const Vector3 &vNormal,
// 		const Vector3 &vEstX, Vector3 &vOut1,
// 		Vector3 &vOut2, bool bInputIsNormalized) {
// 	Vector3 n(vNormal);
// 	Vector3 tan2(vEstX);
// 	if (!bInputIsNormalized) {
// 		n.normalize();
// 		tan2.normalize();
// 	}
// 	Vector3 tan1 = n.cross(tan2.cross(n));
// 	tan1.normalize();
// 	tan2 = n.cross(tan1);

// 	vOut1 = tan2;
// 	vOut2 = tan1;
// }

// template <class real_t>
// real_t VectorAngleR(const Vector2<real_t> &v1, const Vector2<real_t> &v2) {
// 	real_t fDot = Clamp(v1.dot(v2), (real_t)-1.0, (real_t)1.0);
// 	return (real_t)acos(fDot);
// }

real_t VectorAngleR(const Vector3 &v1, const Vector3 &v2) {
	real_t fDot = Clamp(v1.dot(v2), (real_t)-1.0, (real_t)1.0);
	return (real_t)acos(fDot);
}

// template <class real_t>
// real_t VectorAngleD(const Vector2<real_t> &v1, const Vector2<real_t> &v2) {
// 	real_t fDot = Clamp(v1.dot(v2), (real_t)-1.0, (real_t)1.0);
// 	return (real_t)acos(fDot) * Wml::Math<real_t>::RAD_TO_DEG;
// }

// template <class real_t>
// real_t VectorAngleD(const Vector3 &v1, const Vector3 &v2) {
// 	real_t fDot = Clamp(v1.dot(v2), (real_t)-1.0, (real_t)1.0);
// 	return (real_t)acos(fDot) * Wml::Math<real_t>::RAD_TO_DEG;
// }

// template <class real_t>
// real_t VectorCot(const Vector3 &v1, const Vector3 &v2) {
// 	real_t fDot = v1.dot(v2);

// 	real_t lensqr1 = v1.squaredNorm();
// 	real_t lensqr2 = v2.squaredNorm();
// 	real_t d = Clamp(lensqr1 * lensqr2 - fDot * fDot, (real_t)0.0, (real_t)std::numeric_limits<real_t>::max());
// 	if (d < std::numeric_limits<real_t>::epsilon()) {
// 		return 0;
// 	} else {
// 		return fDot / sqrt(d);
// 	}
// }

// template <class real_t>
// Vector2<real_t> Lerp(const Vector2<real_t> &v1, const Vector2<real_t> &v2,
// 		real_t t) {
// 	return (1 - t) * v1 + (t)*v2;
// }
// template <class real_t>
// Vector3 Lerp(const Vector3 &v1, const Vector3 &v2,
// 		real_t t) {
// 	return (1 - t) * v1 + (t)*v2;
// }

// template <class real_t>
// Vector2<real_t> ToUV(const Vector3 &vec, int nUIndex, int nVIndex) {
// 	return Vector2<real_t>(vec[nUIndex], vec[nVIndex]);
// }

// template <class real_t>
// Vector3 To3D(const Vector2<real_t> &vec, int nUIndex, int nVIndex) {
// 	Vector3 tmp = Vector3::Zero();
// 	tmp[nUIndex] = vec.x();
// 	tmp[nVIndex] = vec.y();
// 	return tmp;
// }

// template <class real_t>
// void ToGLMatrix(const Basis &matrix, real_t glMatrix[16]) {
// 	for (int r = 0; r < 3; ++r)
// 		for (int c = 0; c < 4; ++c)
// 			glMatrix[c * 4 + r] = (c < 3) ? matrix(r, c) : 0;
// 	glMatrix[3] = glMatrix[7] = glMatrix[11] = 0;
// 	glMatrix[15] = 1;
// }

// template <class real_t>
// void BarycentricCoords(const Vector3 &vTriVtx1,
// 		const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3,
// 		const Vector3 &vVertex, real_t &fBary1,
// 		real_t &fBary2, real_t &fBary3) {
// 	Vector3 kV02 = vTriVtx1 - vTriVtx3;
// 	Vector3 kV12 = vTriVtx2 - vTriVtx3;
// 	Vector3 kPV2 = vVertex - vTriVtx3;

// 	real_t fM00 = kV02.dot(kV02);
// 	real_t fM01 = kV02.dot(kV12);
// 	real_t fM11 = kV12.dot(kV12);
// 	real_t fR0 = kV02.dot(kPV2);
// 	real_t fR1 = kV12.dot(kPV2);
// 	real_t fDet = fM00 * fM11 - fM01 * fM01;
// 	//    lgASSERT( Math<real_t>::FAbs(fDet) > (real_t)0.0 );
// 	real_t fInvDet = ((real_t)1.0) / fDet;

// 	fBary1 = (fM11 * fR0 - fM01 * fR1) * fInvDet;
// 	fBary2 = (fM00 * fR1 - fM01 * fR0) * fInvDet;
// 	fBary3 = (real_t)1.0 - fBary1 - fBary2;
// }

// template <class real_t>
// real_t Area(const Vector3 &vTriVtx1, const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3) {
// 	Vector3 edge1(vTriVtx2 - vTriVtx1);
// 	Vector3 edge2(vTriVtx3 - vTriVtx1);
// 	Vector3 vCross(edge1.cross(edge2));

// 	return (real_t)0.5 * vCross.norm();
// }

// template <class real_t>
// void BarycentricCoords(const Vector2<real_t> &vTriVtx1,
// 		const Vector2<real_t> &vTriVtx2,
// 		const Vector2<real_t> &vTriVtx3,
// 		const Vector2<real_t> &vVertex, real_t &fBary1,
// 		real_t &fBary2, real_t &fBary3) {
// 	Vector2<real_t> kV02 = vTriVtx1 - vTriVtx3;
// 	Vector2<real_t> kV12 = vTriVtx2 - vTriVtx3;
// 	Vector2<real_t> kPV2 = vVertex - vTriVtx3;

// 	real_t fM00 = kV02.dot(kV02);
// 	real_t fM01 = kV02.dot(kV12);
// 	real_t fM11 = kV12.dot(kV12);
// 	real_t fR0 = kV02.dot(kPV2);
// 	real_t fR1 = kV12.dot(kPV2);
// 	real_t fDet = fM00 * fM11 - fM01 * fM01;
// 	//    lgASSERT( Math<real_t>::FAbs(fDet) > (real_t)0.0 );
// 	real_t fInvDet = ((real_t)1.0) / fDet;

// 	fBary1 = (fM11 * fR0 - fM01 * fR1) * fInvDet;
// 	fBary2 = (fM00 * fR1 - fM01 * fR0) * fInvDet;
// 	fBary3 = (real_t)1.0 - fBary1 - fBary2;
// }

// template <class real_t>
// real_t Area(const Vector2<real_t> &vTriVtx1, const Vector2<real_t> &vTriVtx2,
// 		const Vector2<real_t> &vTriVtx3) {
// 	Vector2<real_t> edge1(vTriVtx2 - vTriVtx1);
// 	Vector2<real_t> edge2(vTriVtx3 - vTriVtx1);
// 	real_t fDot = edge1.dot(edge2);
// 	return (real_t)0.5 *
// 		   sqrt(edge1.squaredNorm() * edge2.squaredNorm() - fDot * fDot);
// }

// template <class real_t>
// Vector3 Normal(const Vector3 &vTriVtx1,
// 		const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3, real_t *pArea) {
// 	Vector3 edge1(vTriVtx2 - vTriVtx1);
// 	Vector3 edge2(vTriVtx3 - vTriVtx1);
// 	if (pArea) {
// 		real_t fDot = edge1.dot(edge2);
// 		*pArea = (real_t)0.5 *
// 				 sqrt(edge1.squaredNorm() * edge2.squaredNorm() - fDot * fDot);
// 	}
// 	edge1.normalize();
// 	edge2.normalize();
// 	Vector3 vCross(edge1.cross(edge2));
// 	vCross.normalize();
// 	return vCross;
// }

// template <class real_t>
// Vector3
// InterpNormal(const Vector3 &vTriVtx1, const Vector3 &vTriVtx2,
// 		const Vector3 &vTriVtx3, const Vector3 &vTriNorm1,
// 		const Vector3 &vTriNorm2, const Vector3 &vTriNorm3,
// 		const Vector3 &vPointInTri) {
// 	real_t fBary[3];
// 	BarycentricCoords(vTriVtx1, vTriVtx2, vTriVtx3, vPointInTri, fBary[0],
// 			fBary[1], fBary[2]);
// 	Vector3 vNormal(fBary[0] * vTriNorm1 + fBary[1] * vTriNorm1 +
// 						  fBary[2] * vTriNorm1);
// 	vNormal.normalize();
// 	return vNormal;
// }

// template <class real_t>
// void StretchMetric1(const Vector3 &q1, const Vector3 &q2,
// 		const Vector3 &q3, const Vector2<real_t> &p1,
// 		const Vector2<real_t> &p2, const Vector2<real_t> &p3,
// 		real_t &MaxSV, real_t &MinSV, real_t &L2Norm,
// 		real_t &LInfNorm) {
// 	real_t s1 = p1.x();
// 	real_t t1 = p1.y();
// 	real_t s2 = p2.x();
// 	real_t t2 = p2.y();
// 	real_t s3 = p3.x();
// 	real_t t3 = p3.y();

// 	real_t A = (real_t)0.5 * ((s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1));
// 	if (A > 0) {
// 		Vector3 Ss =
// 				(q1 * (t2 - t3) + q2 * (t3 - t1) + q3 * (t1 - t2)) / (2 * A);
// 		Vector3 St =
// 				(q1 * (s3 - s2) + q2 * (s1 - s3) + q3 * (s2 - s1)) / (2 * A);

// 		real_t a = Ss.dot(Ss);
// 		real_t b = Ss.dot(St);
// 		real_t c = St.dot(St);

// 		real_t discrim = (real_t)sqrt((a - c) * (a - c) + 4 * b * b);

// 		MaxSV = (real_t)sqrt((real_t)0.5 * ((a + c) + discrim));
// 		MinSV = (real_t)sqrt((real_t)0.5 * ((a + c) - discrim));

// 		L2Norm = (real_t)sqrt((real_t)0.5 * (a + c));
// 		LInfNorm = MaxSV;
// 	} else {
// 		MaxSV = MinSV = L2Norm = LInfNorm = std::numeric_limits<real_t>::max();
// 	}
// }

// template <class real_t>
// void StretchMetric3(const Vector3 &q1, const Vector3 &q2,
// 		const Vector3 &q3, const Vector3 &p1_3D,
// 		const Vector3 &p2_3D, const Vector3 &p3_3D,
// 		real_t &MaxSV, real_t &MinSV, real_t &L2Norm,
// 		real_t &LInfNorm) {
// 	// compute plane containing p1/2/3
// 	Vector3 e1(p2_3D - p1_3D);
// 	e1.normalize();
// 	Vector3 e2(p3_3D - p1_3D);
// 	e2.normalize();
// 	Vector3 n(e1.cross(e2));
// 	n.normalize();
// 	e2 = n.cross(e1);
// 	e2.normalize();

// 	Vector2<real_t> p1(Vector2<real_t>::Zero());
// 	Vector2<real_t> p2((p2_3D - p1_3D).dot(e1), (p2_3D - p1_3D).dot(e2));
// 	Vector2<real_t> p3((p3_3D - p1_3D).dot(e1), (p3_3D - p1_3D).dot(e2));

// 	real_t s1 = p1.x();
// 	real_t t1 = p1.y();
// 	real_t s2 = p2.x();
// 	real_t t2 = p2.y();
// 	real_t s3 = p3.x();
// 	real_t t3 = p3.y();

// 	real_t A = (real_t)0.5 * ((s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1));
// 	if (A > 0) {
// 		Vector3 Ss =
// 				(q1 * (t2 - t3) + q2 * (t3 - t1) + q3 * (t1 - t2)) / (2 * A);
// 		Vector3 St =
// 				(q1 * (s3 - s2) + q2 * (s1 - s3) + q3 * (s2 - s1)) / (2 * A);

// 		real_t a = Ss.dot(Ss);
// 		real_t b = Ss.dot(St);
// 		real_t c = St.dot(St);

// 		real_t discrim = (real_t)sqrt((a - c) * (a - c) + 4 * b * b);

// 		MaxSV = (real_t)sqrt((real_t)0.5 * ((a + c) + discrim));
// 		MinSV = (real_t)sqrt((real_t)0.5 * ((a + c) - discrim));

// 		L2Norm = (real_t)sqrt((real_t)0.5 * (a + c));
// 		LInfNorm = MaxSV;
// 	} else {
// 		MaxSV = MinSV = L2Norm = LInfNorm = std::numeric_limits<real_t>::max();
// 	}
// }

// template <class real_t>
// bool IsObtuse(const Vector2<real_t> &v1, const Vector2<real_t> &v2,
// 		const Vector2<real_t> &v3) {
// 	// from http://mathworld.wolfram.com/ObtuseTriangle.html
// 	real_t a2 = (v1 - v2).squaredNorm();
// 	real_t b2 = (v1 - v3).squaredNorm();
// 	real_t c2 = (v2 - v3).squaredNorm();
// 	return (a2 + b2 < c2) || (b2 + c2 < a2) || (c2 + a2 < b2);
// }
// template <class real_t>
// bool IsObtuse(const Vector3 &v1, const Vector3 &v2,
// 		const Vector3 &v3) {
// 	real_t a2 = (v1 - v2).squaredNorm();
// 	real_t b2 = (v1 - v3).squaredNorm();
// 	real_t c2 = (v2 - v3).squaredNorm();
// 	return (a2 + b2 < c2) || (b2 + c2 < a2) || (c2 + a2 < b2);
// }

// template void ToGLMatrix(const Matrix3<float> &matrix,
// 		float glMatrix[16]);
// template void ToGLMatrix(const Matrix3<double> &matrix,
// 		double glMatrix[16]);

// template void
// ComputeAlignZAxisMatrix(const Vector3<float> &vAlignWith,
// 		Matrix3<float> &matrix, bool bInvert);
// template void
// ComputeAlignZAxisMatrix(const Vector3<double> &vAlignWith,
// 		Matrix3<double> &matrix, bool bInvert);

// template void
// ComputeAlignAxisMatrix(const Vector3<float> &vInitial,
// 		const Vector3<float> &vAlignWith,
// 		Matrix3<float> &matrix);
// template void
// ComputeAlignAxisMatrix(const Vector3<double> &vInitial,
// 		const Vector3<double> &vAlignWith,
// 		Matrix3<double> &matrix);

// template void ComputePerpVectors(const Vector3<float> &vIn,
// 		Vector3<float> &vOut1,
// 		Vector3<float> &vOut2,
// 		bool bInIsNormalized);
// template void ComputePerpVectors(const Vector3<double> &vIn,
// 		Vector3<double> &vOut1,
// 		Vector3<double> &vOut2,
// 		bool bInIsNormalized);

// template void ComputePerpVectors(const Vector3<float> &vNormal,
// 		const Vector3<float> &vEstX,
// 		Vector3<float> &vOut1,
// 		Vector3<float> &vOut2,
// 		bool bInputIsNormalized);
// template void ComputePerpVectors(const Vector3<double> &vNormal,
// 		const Vector3<double> &vEstX,
// 		Vector3<double> &vOut1,
// 		Vector3<double> &vOut2,
// 		bool bInputIsNormalized);

// template float VectorAngleR(const Vector2<float> &v1,
// 		const Vector2<float> &v2);
// template double VectorAngleR(const Vector2<double> &v1,
// 		const Vector2<double> &v2);
// template float VectorAngleR(const Vector3<float> &v1,
// 		const Vector3<float> &v2);
// template double VectorAngleR(const Vector3<double> &v1,
// 		const Vector3<double> &v2);

// template float VectorAngleD(const Vector2<float> &v1,
// 		const Vector2<float> &v2);
// template double VectorAngleD(const Vector2<double> &v1,
// 		const Vector2<double> &v2);
// template float VectorAngleD(const Vector3<float> &v1,
// 		const Vector3<float> &v2);
// template double VectorAngleD(const Vector3<double> &v1,
// 		const Vector3<double> &v2);

// template float VectorCot(const Vector3<float> &v1,
// 		const Vector3<float> &v2);
// template double VectorCot(const Vector3<double> &v1,
// 		const Vector3<double> &v2);

// template Vector2<float> Lerp(const Vector2<float> &v1,
// 		const Vector2<float> &v2, float t);
// template Vector2<double>
// Lerp(const Vector2<double> &v1, const Vector2<double> &v2, double t);
// template Vector3<float> Lerp(const Vector3<float> &v1,
// 		const Vector3<float> &v2, float t);
// template Vector3<double>
// Lerp(const Vector3<double> &v1, const Vector3<double> &v2, double t);

// template Vector2<float> ToUV(const Vector3<float> &vec, int nUIndex,
// 		int nVIndex);
// template Vector2<double> ToUV(const Vector3<double> &vec,
// 		int nUIndex, int nVIndex);
// template Vector3<float> To3D(const Vector2<float> &vec, int nUIndex,
// 		int nVIndex);
// template Vector3<double> To3D(const Vector2<double> &vec,
// 		int nUIndex, int nVIndex);

// template void
// BarycentricCoords(const Vector3<float> &TriVtx1, const Vector3<float> &TriVtx2,
// 		const Vector3<float> &TriVtx3, const Vector3<float> &vVertex,
// 		float &fWeight1, float &fWeight2, float &fWeight3);
// template void BarycentricCoords(const Vector3<double> &TriVtx1,
// 		const Vector3<double> &TriVtx2,
// 		const Vector3<double> &TriVtx3,
// 		const Vector3<double> &vVertex,
// 		double &fWeight1, double &fWeight2,
// 		double &fWeight3);
// template void
// BarycentricCoords(const Vector2<float> &TriVtx1, const Vector2<float> &TriVtx2,
// 		const Vector2<float> &TriVtx3, const Vector2<float> &vVertex,
// 		float &fWeight1, float &fWeight2, float &fWeight3);
// template void BarycentricCoords(const Vector2<double> &TriVtx1,
// 		const Vector2<double> &TriVtx2,
// 		const Vector2<double> &TriVtx3,
// 		const Vector2<double> &vVertex,
// 		double &fWeight1, double &fWeight2,
// 		double &fWeight3);

// template float Area(const Vector3<float> &TriVtx1,
// 		const Vector3<float> &TriVtx2,
// 		const Vector3<float> &TriVtx3);
// template double Area(const Vector3<double> &TriVtx1,
// 		const Vector3<double> &TriVtx2,
// 		const Vector3<double> &TriVtx3);
// template float Area(const Vector2<float> &TriVtx1,
// 		const Vector2<float> &TriVtx2,
// 		const Vector2<float> &TriVtx3);
// template double Area(const Vector2<double> &TriVtx1,
// 		const Vector2<double> &TriVtx2,
// 		const Vector2<double> &TriVtx3);

// template Vector3<float> Normal(const Vector3<float> &TriVtx1,
// 		const Vector3<float> &TriVtx2,
// 		const Vector3<float> &TriVtx3,
// 		float *pArea);
// template Vector3<double> Normal(const Vector3<double> &TriVtx1,
// 		const Vector3<double> &TriVtx2,
// 		const Vector3<double> &TriVtx3,
// 		double *pArea);

// template Vector3<float>
// InterpNormal(const Vector3<float> &vTriVtx1, const Vector3<float> &vTriVtx2,
// 		const Vector3<float> &vTriVtx3, const Vector3<float> &vTriNorm1,
// 		const Vector3<float> &vTriNorm2, const Vector3<float> &vTriNorm3,
// 		const Vector3<float> &vPointInTri);
// template Vector3<double>
// InterpNormal(const Vector3<double> &vTriVtx1, const Vector3<double> &vTriVtx2,
// 		const Vector3<double> &vTriVtx3, const Vector3<double> &vTriNorm1,
// 		const Vector3<double> &vTriNorm2, const Vector3<double> &vTriNorm3,
// 		const Vector3<double> &vPointInTri);

// template void
// StretchMetric1(const Vector3<float> &q1, const Vector3<float> &q2,
// 		const Vector3<float> &q3, const Vector2<float> &p1,
// 		const Vector2<float> &p2, const Vector2<float> &p3, float &MaxSV,
// 		float &MinSV, float &L2Norm, float &LInfNorm);
// template void
// StretchMetric1(const Vector3<double> &q1, const Vector3<double> &q2,
// 		const Vector3<double> &q3, const Vector2<double> &p1,
// 		const Vector2<double> &p2, const Vector2<double> &p3,
// 		double &MaxSV, double &MinSV, double &L2Norm, double &LInfNorm);

// template void
// StretchMetric3(const Vector3<float> &q1, const Vector3<float> &q2,
// 		const Vector3<float> &q3, const Vector3<float> &p1,
// 		const Vector3<float> &p2, const Vector3<float> &p3, float &MaxSV,
// 		float &MinSV, float &L2Norm, float &LInfNorm);
// template void
// StretchMetric3(const Vector3<double> &q1, const Vector3<double> &q2,
// 		const Vector3<double> &q3, const Vector3<double> &p1,
// 		const Vector3<double> &p2, const Vector3<double> &p3,
// 		double &MaxSV, double &MinSV, double &L2Norm, double &LInfNorm);

// template bool IsObtuse(const Vector2<float> &v1,
// 		const Vector2<float> &v2,
// 		const Vector2<float> &v3);
// template bool IsObtuse(const Vector2<double> &v1,
// 		const Vector2<double> &v2,
// 		const Vector2<double> &v3);
// template bool IsObtuse(const Vector3<float> &v1,
// 		const Vector3<float> &v2,
// 		const Vector3<float> &v3);
// template bool IsObtuse(const Vector3<double> &v1,
// 		const Vector3<double> &v2,
// 		const Vector3<double> &v3);
// } // namespace g3

#pragma once

#include <corecrt_math.h>
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
template <class Real>
void ComputeAlignZAxisMatrix(const Vector3<Real> &vAlignWith,
		Matrix3<Real> &matrix, bool bToZAxis = false);

template <class Real>
void ComputeAlignAxisMatrix(const Vector3<Real> &vInitial,
		const Vector3<Real> &vAlignWith, Matrix3<Real> &matrix);

//! compute vectors in a plane perpendicular to vIn
template <class Real>
void ComputePerpVectors(const Vector3<Real> &vIn,
		Vector3<Real> &vOut1, Vector3<Real> &vOut2,
		bool bInIsNormalized = false);

//! compute tangent vectors in plane perp to vNormal, using non-orthogonal vEstX as estimate of vOut1
template <class Real>
void ComputePerpVectors(const Vector3<Real> &vNormal, const Vector3<Real> &vEstX,
		Vector3<Real> &vOut1, Vector3<Real> &vOut2,
		bool bInputIsNormalized = false);

template <class Real>
void ToGLMatrix(const Matrix3<Real> &matrix, Real glMatrix[16]);

template <class Real>
Vector2<Real> ToUV(const Vector3<Real> &vec, int nUIndex, int nVIndex);
template <class Real>
Vector3<Real> To3D(const Vector2<Real> &vec, int nUIndex, int nVIndex);

template <class Real>
Real VectorAngleR(const Vector2<Real> &v1, const Vector2<Real> &v2);
template <class Real>
Real VectorAngleR(const Vector3<Real> &v1, const Vector3<Real> &v2);

template <class Real>
Real VectorAngleD(const Vector2<Real> &v1, const Vector2<Real> &v2);
template <class Real>
Real VectorAngleD(const Vector3<Real> &v1, const Vector3<Real> &v2);

template <class Real>
Vector2<Real> Lerp(const Vector2<Real> &v1, const Vector2<Real> &v2, Real t);
template <class Real>
Vector3<Real> Lerp(const Vector3<Real> &v1, const Vector3<Real> &v2, Real t);

template <class Real>
void BarycentricCoords(const Vector3<Real> &vTriVtx1,
		const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3,
		const Vector3<Real> &vVertex,
		Real &fBary1, Real &fBary2, Real &fBary3);

template <class Real>
Real Area(const Vector3<Real> &vTriVtx1,
		const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3);

template <class Real>
void BarycentricCoords(const Vector2<Real> &vTriVtx1,
		const Vector2<Real> &vTriVtx2,
		const Vector2<Real> &vTriVtx3,
		const Vector2<Real> &vVertex,
		Real &fBary1, Real &fBary2, Real &fBary3);

template <class Real>
Real Area(const Vector2<Real> &vTriVtx1,
		const Vector2<Real> &vTriVtx2,
		const Vector2<Real> &vTriVtx3);

template <class Real>
Vector3<Real> Normal(const Vector3<Real> &vTriVtx1,
		const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3, Real *pArea = nullptr);

template <class Real>
Vector3<Real> InterpNormal(const Vector3<Real> &vTriVtx1,
		const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3,
		const Vector3<Real> &vTriNorm1,
		const Vector3<Real> &vTriNorm2,
		const Vector3<Real> &vTriNorm3,
		const Vector3<Real> &vPointInTri);

//! This metric is from Texture Mapping Progressive Meshes, Sander et al, Siggraph 2001
template <class Real>
void StretchMetric1(const Vector3<Real> &vTriVtx1,
		const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3,
		const Vector2<Real> &vVtxParam1,
		const Vector2<Real> &vVtxParam2,
		const Vector2<Real> &vVtxParam3,
		Real &MaxSV, Real &MinSV, Real &L2Norm, Real &LInfNorm);

template <class Real>
void StretchMetric3(const Vector3<Real> &vTriVtx1,
		const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3,
		const Vector3<Real> &vVtxParam1,
		const Vector3<Real> &vVtxParam2,
		const Vector3<Real> &vVtxParam3,
		Real &MaxSV, Real &MinSV, Real &L2Norm, Real &LInfNorm);

template <class Real>
bool IsObtuse(const Vector2<Real> &v1, const Vector2<Real> &v2, const Vector2<Real> &v3);
template <class Real>
bool IsObtuse(const Vector3<Real> &v1, const Vector3<Real> &v2, const Vector3<Real> &v3);

/*
	 * inline utilities
	 */

inline bool IsFinite(const Vector3d &v) {
	return isfinite(v.x()) && isfinite(v.y()) && isfinite(v.z());
}

inline Vector2f d2f(const Vector2d &v) {
	return Vector2f((float)v[0], (float)v[1]);
}
inline Vector3f d2f(const Vector3d &v) {
	return Vector3f((float)v[0], (float)v[1], (float)v[2]);
}
inline Vector2d f2d(const Vector2f &v) {
	return Vector2d((double)v[0], (double)v[1]);
}
inline Vector3d f2d(const Vector3f &v) {
	return Vector3d((double)v[0], (double)v[1], (double)v[2]);
}

inline Matrix2f d2f(const Matrix2d &v) {
	return v.cast<float>();
}
inline Matrix3f d2f(const Matrix3d &v) {
	return v.cast<float>();
}
inline Matrix3d f2d(const Matrix3f &v) {
	return v.cast<double>();
}
inline Matrix2d f2d(const Matrix2f &v) {
	return v.cast<double>();
}

//inline AxisAlignedBox2f d2f(const AxisAlignedBox2d & v) {
//	return AxisAlignedBox2f((float)v.Min[0], (float)v.Max[0], (float)v.Min[1], (float)v.Max[1] );
//}
//inline AxisAlignedBox3f d2f(const AxisAlignedBox3d & v) {
//	return AxisAlignedBox3f((float)v.Min[0], (float)v.Max[0], (float)v.Min[1], (float)v.Max[1], (float)v.Min[2], (float)v.Max[2] );
//}
//inline AxisAlignedBox2d f2d(const AxisAlignedBox2f & v) {
//	return AxisAlignedBox2d((double)v.Min[0], (double)v.Max[0], (double)v.Min[1], (double)v.Max[1]);
//}
//inline AxisAlignedBox3d f2d(const AxisAlignedBox3f & v) {
//	return AxisAlignedBox3d((double)v.Min[0], (double)v.Max[0], (double)v.Min[1], (double)v.Max[1], (double)v.Min[2], (double)v.Max[2]);
//}

template <class Real>
inline Real Clamp(const Real &fValue, const Real &fMin, const Real &fMax) {
	if (fValue < fMin)
		return fMin;
	else if (fValue > fMax)
		return fMax;
	else
		return fValue;
}

inline void array3f_add(float *pBuffer, unsigned int nIndex, const float *pAdd) {
	pBuffer[3 * nIndex] += pAdd[0];
	pBuffer[3 * nIndex + 1] += pAdd[1];
	pBuffer[3 * nIndex + 2] += pAdd[2];
}
inline void array3f_normalize(float *pBuffer, unsigned int nIndex, float fEpsilon = 0.0f) {
	auto v = Vector3f(&pBuffer[3 * nIndex]).normalized();
	pBuffer[3 * nIndex] = v[0];
	pBuffer[3 * nIndex + 1] = v[1];
	pBuffer[3 * nIndex + 2] = v[2];
}
inline void vectorf_push(std::vector<float> &v, const Vector3f &p) {
	v.push_back(p[0]);
	v.push_back(p[1]);
	v.push_back(p[2]);
}
inline void vectori_push(std::vector<unsigned int> &v, const Vector2i &p) {
	v.push_back(p[0]);
	v.push_back(p[1]);
}
inline void vectori_push(std::vector<unsigned int> &v, const Vector3i &p) {
	v.push_back(p[0]);
	v.push_back(p[1]);
	v.push_back(p[2]);
}

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

template <class Real>
void ComputeAlignZAxisMatrix(const Vector3<Real> &vAlignWith,
		Matrix3<Real> &matrix, bool bInvert) {
	// compute cosine of angle between vectors
	Real axisDot = vAlignWith.dot(Vector3<Real>::UnitZ());

	// compute rotation axis
	Vector3<Real> axisCross(Vector3<Real>::UnitZ().cross(vAlignWith));

	Real fInverter = (bInvert) ? (Real)-1 : (Real)1;

	// apply rotation if necessary
	if (axisCross.squaredNorm() > Math<Real>::EPSILON) {
		// compute normalized axis and angle, then create rotation around axis
		axisCross.normalize();
		Real fAngle = Math<Real>::ACos(axisDot / vAlignWith.norm());
		matrix = Wml::Matrix3<Real>(axisCross, fAngle * fInverter);

	} else if (axisDot < (Real)0) {
		axisCross = Vector3<Real>::UnitX();
		Real fAngle = (Real)180 * Math<Real>::DEG_TO_RAD * fInverter;
		matrix = Wml::Matrix3<Real>(axisCross, fAngle);
	} else {
		matrix = Matrix3<Real>::Identity();
	}
}

template <class Real>
void ComputeAlignAxisMatrix(const Vector3<Real> &vInitial,
		const Vector3<Real> &vAlignWith,
		Matrix3<Real> &matrix) {
	// compute cosine of angle between vectors
	Real axisDot = vAlignWith.dot(vInitial);

	// compute rotation axis
	Vector3<Real> axisCross(vInitial.cross(vAlignWith));

	// apply rotation if necessary
	if (axisCross.squaredNorm() > Math<Real>::EPSILON) {
		// compute normalized axis and angle, then create rotation around axis
		axisCross.normalize();
		Real fAngle = Math<Real>::ACos(axisDot / vAlignWith.norm());
		matrix = Wml::Matrix3<Real>(axisCross, fAngle);

	} else if (axisDot < (Real)0) {
		// find some perpendicular vectors
		Vector3<Real> vPerp1, vPerp2;
		ComputePerpVectors(vInitial, vPerp1, vPerp2);

		matrix = Wml::Matrix3<Real>(vPerp1, (Real)180 * Math<Real>::DEG_TO_RAD);
	} else {
		matrix = Matrix3<Real>::Identity();
	}
}

template <class Real>
void ComputePerpVectors(const Vector3<Real> &vIn, Vector3<Real> &vOut1,
		Vector3<Real> &vOut2, bool bInIsNormalized) {
	Vector3<Real> vPerp(vIn);
	if (!bInIsNormalized)
		vPerp.normalize();

	if (Math<Real>::FAbs(vPerp.x()) >= Math<Real>::FAbs(vPerp.y()) &&
			Math<Real>::FAbs(vPerp.x()) >= Math<Real>::FAbs(vPerp.z())) {
		vOut1.x() = -vPerp.y();
		vOut1.y() = vPerp.x();
		vOut1.z() = (Real)0.0;
	} else {
		vOut1.x() = (Real)0.0;
		vOut1.y() = vPerp.z();
		vOut1.z() = -vPerp.y();
	}

	vOut1.normalize();
	vOut2 = vPerp.cross(vOut1);
}

template <class Real>
void ComputePerpVectors(const Vector3<Real> &vNormal,
		const Vector3<Real> &vEstX, Vector3<Real> &vOut1,
		Vector3<Real> &vOut2, bool bInputIsNormalized) {
	Vector3<Real> n(vNormal);
	Vector3<Real> tan2(vEstX);
	if (!bInputIsNormalized) {
		n.normalize();
		tan2.normalize();
	}
	Vector3<Real> tan1 = n.cross(tan2.cross(n));
	tan1.normalize();
	tan2 = n.cross(tan1);

	vOut1 = tan2;
	vOut2 = tan1;
}

template <class Real>
Real VectorAngleR(const Vector2<Real> &v1, const Vector2<Real> &v2) {
	Real fDot = Clamp(v1.dot(v2), (Real)-1.0, (Real)1.0);
	return (Real)acos(fDot);
}

template <class Real>
Real VectorAngleR(const Vector3<Real> &v1, const Vector3<Real> &v2) {
	Real fDot = Clamp(v1.dot(v2), (Real)-1.0, (Real)1.0);
	return (Real)acos(fDot);
}

template <class Real>
Real VectorAngleD(const Vector2<Real> &v1, const Vector2<Real> &v2) {
	Real fDot = Clamp(v1.dot(v2), (Real)-1.0, (Real)1.0);
	return (Real)acos(fDot) * Wml::Math<Real>::RAD_TO_DEG;
}

template <class Real>
Real VectorAngleD(const Vector3<Real> &v1, const Vector3<Real> &v2) {
	Real fDot = Clamp(v1.dot(v2), (Real)-1.0, (Real)1.0);
	return (Real)acos(fDot) * Wml::Math<Real>::RAD_TO_DEG;
}

template <class Real>
Vector2<Real> Lerp(const Vector2<Real> &v1, const Vector2<Real> &v2,
		Real t) {
	return (1 - t) * v1 + (t)*v2;
}
template <class Real>
Vector3<Real> Lerp(const Vector3<Real> &v1, const Vector3<Real> &v2,
		Real t) {
	return (1 - t) * v1 + (t)*v2;
}

template <class Real>
Vector2<Real> ToUV(const Vector3<Real> &vec, int nUIndex, int nVIndex) {
	return Vector2<Real>(vec[nUIndex], vec[nVIndex]);
}

template <class Real>
Vector3<Real> To3D(const Vector2<Real> &vec, int nUIndex, int nVIndex) {
	Vector3<Real> tmp = Vector3<Real>::Zero();
	tmp[nUIndex] = vec.x();
	tmp[nVIndex] = vec.y();
	return tmp;
}

template <class Real>
void ToGLMatrix(const Matrix3<Real> &matrix, Real glMatrix[16]) {
	for (int r = 0; r < 3; ++r)
		for (int c = 0; c < 4; ++c)
			glMatrix[c * 4 + r] = (c < 3) ? matrix(r, c) : 0;
	glMatrix[3] = glMatrix[7] = glMatrix[11] = 0;
	glMatrix[15] = 1;
}

template <class Real>
void BarycentricCoords(const Vector3<Real> &vTriVtx1,
		const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3,
		const Vector3<Real> &vVertex, Real &fBary1,
		Real &fBary2, Real &fBary3) {
	Vector3<Real> kV02 = vTriVtx1 - vTriVtx3;
	Vector3<Real> kV12 = vTriVtx2 - vTriVtx3;
	Vector3<Real> kPV2 = vVertex - vTriVtx3;

	Real fM00 = kV02.dot(kV02);
	Real fM01 = kV02.dot(kV12);
	Real fM11 = kV12.dot(kV12);
	Real fR0 = kV02.dot(kPV2);
	Real fR1 = kV12.dot(kPV2);
	Real fDet = fM00 * fM11 - fM01 * fM01;
	//    lgASSERT( Math<Real>::FAbs(fDet) > (Real)0.0 );
	Real fInvDet = ((Real)1.0) / fDet;

	fBary1 = (fM11 * fR0 - fM01 * fR1) * fInvDet;
	fBary2 = (fM00 * fR1 - fM01 * fR0) * fInvDet;
	fBary3 = (Real)1.0 - fBary1 - fBary2;
}

template <class Real>
Real Area(const Vector3<Real> &vTriVtx1, const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3) {
	Vector3<Real> edge1(vTriVtx2 - vTriVtx1);
	Vector3<Real> edge2(vTriVtx3 - vTriVtx1);
	Vector3<Real> vCross(edge1.cross(edge2));

	return (Real)0.5 * vCross.norm();
}

template <class Real>
void BarycentricCoords(const Vector2<Real> &vTriVtx1,
		const Vector2<Real> &vTriVtx2,
		const Vector2<Real> &vTriVtx3,
		const Vector2<Real> &vVertex, Real &fBary1,
		Real &fBary2, Real &fBary3) {
	Vector2<Real> kV02 = vTriVtx1 - vTriVtx3;
	Vector2<Real> kV12 = vTriVtx2 - vTriVtx3;
	Vector2<Real> kPV2 = vVertex - vTriVtx3;

	Real fM00 = kV02.dot(kV02);
	Real fM01 = kV02.dot(kV12);
	Real fM11 = kV12.dot(kV12);
	Real fR0 = kV02.dot(kPV2);
	Real fR1 = kV12.dot(kPV2);
	Real fDet = fM00 * fM11 - fM01 * fM01;
	//    lgASSERT( Math<Real>::FAbs(fDet) > (Real)0.0 );
	Real fInvDet = ((Real)1.0) / fDet;

	fBary1 = (fM11 * fR0 - fM01 * fR1) * fInvDet;
	fBary2 = (fM00 * fR1 - fM01 * fR0) * fInvDet;
	fBary3 = (Real)1.0 - fBary1 - fBary2;
}

template <class Real>
Real Area(const Vector2<Real> &vTriVtx1, const Vector2<Real> &vTriVtx2,
		const Vector2<Real> &vTriVtx3) {
	Vector2<Real> edge1(vTriVtx2 - vTriVtx1);
	Vector2<Real> edge2(vTriVtx3 - vTriVtx1);
	Real fDot = edge1.dot(edge2);
	return (Real)0.5 *
		   sqrt(edge1.squaredNorm() * edge2.squaredNorm() - fDot * fDot);
}

template <class Real>
Vector3<Real> Normal(const Vector3<Real> &vTriVtx1,
		const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3, Real *pArea) {
	Vector3<Real> edge1(vTriVtx2 - vTriVtx1);
	Vector3<Real> edge2(vTriVtx3 - vTriVtx1);
	if (pArea) {
		Real fDot = edge1.dot(edge2);
		*pArea = (Real)0.5 *
				 sqrt(edge1.squaredNorm() * edge2.squaredNorm() - fDot * fDot);
	}
	edge1.normalize();
	edge2.normalize();
	Vector3<Real> vCross(edge1.cross(edge2));
	vCross.normalize();
	return vCross;
}

template <class Real>
Vector3<Real>
InterpNormal(const Vector3<Real> &vTriVtx1, const Vector3<Real> &vTriVtx2,
		const Vector3<Real> &vTriVtx3, const Vector3<Real> &vTriNorm1,
		const Vector3<Real> &vTriNorm2, const Vector3<Real> &vTriNorm3,
		const Vector3<Real> &vPointInTri) {
	Real fBary[3];
	BarycentricCoords(vTriVtx1, vTriVtx2, vTriVtx3, vPointInTri, fBary[0],
			fBary[1], fBary[2]);
	Vector3<Real> vNormal(fBary[0] * vTriNorm1 + fBary[1] * vTriNorm1 +
						  fBary[2] * vTriNorm1);
	vNormal.normalize();
	return vNormal;
}

template <class Real>
void StretchMetric1(const Vector3<Real> &q1, const Vector3<Real> &q2,
		const Vector3<Real> &q3, const Vector2<Real> &p1,
		const Vector2<Real> &p2, const Vector2<Real> &p3,
		Real &MaxSV, Real &MinSV, Real &L2Norm,
		Real &LInfNorm) {
	Real s1 = p1.x();
	Real t1 = p1.y();
	Real s2 = p2.x();
	Real t2 = p2.y();
	Real s3 = p3.x();
	Real t3 = p3.y();

	Real A = (Real)0.5 * ((s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1));
	if (A > 0) {
		Vector3<Real> Ss =
				(q1 * (t2 - t3) + q2 * (t3 - t1) + q3 * (t1 - t2)) / (2 * A);
		Vector3<Real> St =
				(q1 * (s3 - s2) + q2 * (s1 - s3) + q3 * (s2 - s1)) / (2 * A);

		Real a = Ss.dot(Ss);
		Real b = Ss.dot(St);
		Real c = St.dot(St);

		Real discrim = (Real)sqrt((a - c) * (a - c) + 4 * b * b);

		MaxSV = (Real)sqrt((Real)0.5 * ((a + c) + discrim));
		MinSV = (Real)sqrt((Real)0.5 * ((a + c) - discrim));

		L2Norm = (Real)sqrt((Real)0.5 * (a + c));
		LInfNorm = MaxSV;
	} else {
		MaxSV = MinSV = L2Norm = LInfNorm = std::numeric_limits<Real>::max();
	}
}

template <class Real>
void StretchMetric3(const Vector3<Real> &q1, const Vector3<Real> &q2,
		const Vector3<Real> &q3, const Vector3<Real> &p1_3D,
		const Vector3<Real> &p2_3D, const Vector3<Real> &p3_3D,
		Real &MaxSV, Real &MinSV, Real &L2Norm,
		Real &LInfNorm) {
	// compute plane containing p1/2/3
	Vector3<Real> e1(p2_3D - p1_3D);
	e1.normalize();
	Vector3<Real> e2(p3_3D - p1_3D);
	e2.normalize();
	Vector3<Real> n(e1.cross(e2));
	n.normalize();
	e2 = n.cross(e1);
	e2.normalize();

	Vector2<Real> p1(Vector2<Real>::Zero());
	Vector2<Real> p2((p2_3D - p1_3D).dot(e1), (p2_3D - p1_3D).dot(e2));
	Vector2<Real> p3((p3_3D - p1_3D).dot(e1), (p3_3D - p1_3D).dot(e2));

	Real s1 = p1.x();
	Real t1 = p1.y();
	Real s2 = p2.x();
	Real t2 = p2.y();
	Real s3 = p3.x();
	Real t3 = p3.y();

	Real A = (Real)0.5 * ((s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1));
	if (A > 0) {
		Vector3<Real> Ss =
				(q1 * (t2 - t3) + q2 * (t3 - t1) + q3 * (t1 - t2)) / (2 * A);
		Vector3<Real> St =
				(q1 * (s3 - s2) + q2 * (s1 - s3) + q3 * (s2 - s1)) / (2 * A);

		Real a = Ss.dot(Ss);
		Real b = Ss.dot(St);
		Real c = St.dot(St);

		Real discrim = (Real)sqrt((a - c) * (a - c) + 4 * b * b);

		MaxSV = (Real)sqrt((Real)0.5 * ((a + c) + discrim));
		MinSV = (Real)sqrt((Real)0.5 * ((a + c) - discrim));

		L2Norm = (Real)sqrt((Real)0.5 * (a + c));
		LInfNorm = MaxSV;
	} else {
		MaxSV = MinSV = L2Norm = LInfNorm = std::numeric_limits<Real>::max();
	}
}

template <class Real>
bool IsObtuse(const Vector2<Real> &v1, const Vector2<Real> &v2,
		const Vector2<Real> &v3) {
	// from http://mathworld.wolfram.com/ObtuseTriangle.html
	Real a2 = (v1 - v2).squaredNorm();
	Real b2 = (v1 - v3).squaredNorm();
	Real c2 = (v2 - v3).squaredNorm();
	return (a2 + b2 < c2) || (b2 + c2 < a2) || (c2 + a2 < b2);
}
template <class Real>
bool IsObtuse(const Vector3<Real> &v1, const Vector3<Real> &v2,
		const Vector3<Real> &v3) {
	Real a2 = (v1 - v2).squaredNorm();
	Real b2 = (v1 - v3).squaredNorm();
	Real c2 = (v2 - v3).squaredNorm();
	return (a2 + b2 < c2) || (b2 + c2 < a2) || (c2 + a2 < b2);
}

template void ToGLMatrix(const Matrix3<float> &matrix,
		float glMatrix[16]);
template void ToGLMatrix(const Matrix3<double> &matrix,
		double glMatrix[16]);

template void
ComputeAlignZAxisMatrix(const Vector3<float> &vAlignWith,
		Matrix3<float> &matrix, bool bInvert);
template void
ComputeAlignZAxisMatrix(const Vector3<double> &vAlignWith,
		Matrix3<double> &matrix, bool bInvert);

template void
ComputeAlignAxisMatrix(const Vector3<float> &vInitial,
		const Vector3<float> &vAlignWith,
		Matrix3<float> &matrix);
template void
ComputeAlignAxisMatrix(const Vector3<double> &vInitial,
		const Vector3<double> &vAlignWith,
		Matrix3<double> &matrix);

template void ComputePerpVectors(const Vector3<float> &vIn,
		Vector3<float> &vOut1,
		Vector3<float> &vOut2,
		bool bInIsNormalized);
template void ComputePerpVectors(const Vector3<double> &vIn,
		Vector3<double> &vOut1,
		Vector3<double> &vOut2,
		bool bInIsNormalized);

template void ComputePerpVectors(const Vector3<float> &vNormal,
		const Vector3<float> &vEstX,
		Vector3<float> &vOut1,
		Vector3<float> &vOut2,
		bool bInputIsNormalized);
template void ComputePerpVectors(const Vector3<double> &vNormal,
		const Vector3<double> &vEstX,
		Vector3<double> &vOut1,
		Vector3<double> &vOut2,
		bool bInputIsNormalized);

template float VectorAngleR(const Vector2<float> &v1,
		const Vector2<float> &v2);
template double VectorAngleR(const Vector2<double> &v1,
		const Vector2<double> &v2);
template float VectorAngleR(const Vector3<float> &v1,
		const Vector3<float> &v2);
template double VectorAngleR(const Vector3<double> &v1,
		const Vector3<double> &v2);

template float VectorAngleD(const Vector2<float> &v1,
		const Vector2<float> &v2);
template double VectorAngleD(const Vector2<double> &v1,
		const Vector2<double> &v2);
template float VectorAngleD(const Vector3<float> &v1,
		const Vector3<float> &v2);
template double VectorAngleD(const Vector3<double> &v1,
		const Vector3<double> &v2);

template Vector2<float> Lerp(const Vector2<float> &v1,
		const Vector2<float> &v2, float t);
template Vector2<double>
Lerp(const Vector2<double> &v1, const Vector2<double> &v2, double t);
template Vector3<float> Lerp(const Vector3<float> &v1,
		const Vector3<float> &v2, float t);
template Vector3<double>
Lerp(const Vector3<double> &v1, const Vector3<double> &v2, double t);

template Vector2<float> ToUV(const Vector3<float> &vec, int nUIndex,
		int nVIndex);
template Vector2<double> ToUV(const Vector3<double> &vec,
		int nUIndex, int nVIndex);
template Vector3<float> To3D(const Vector2<float> &vec, int nUIndex,
		int nVIndex);
template Vector3<double> To3D(const Vector2<double> &vec,
		int nUIndex, int nVIndex);

template void
BarycentricCoords(const Vector3<float> &TriVtx1, const Vector3<float> &TriVtx2,
		const Vector3<float> &TriVtx3, const Vector3<float> &vVertex,
		float &fWeight1, float &fWeight2, float &fWeight3);
template void BarycentricCoords(const Vector3<double> &TriVtx1,
		const Vector3<double> &TriVtx2,
		const Vector3<double> &TriVtx3,
		const Vector3<double> &vVertex,
		double &fWeight1, double &fWeight2,
		double &fWeight3);
template void
BarycentricCoords(const Vector2<float> &TriVtx1, const Vector2<float> &TriVtx2,
		const Vector2<float> &TriVtx3, const Vector2<float> &vVertex,
		float &fWeight1, float &fWeight2, float &fWeight3);
template void BarycentricCoords(const Vector2<double> &TriVtx1,
		const Vector2<double> &TriVtx2,
		const Vector2<double> &TriVtx3,
		const Vector2<double> &vVertex,
		double &fWeight1, double &fWeight2,
		double &fWeight3);

template float Area(const Vector3<float> &TriVtx1,
		const Vector3<float> &TriVtx2,
		const Vector3<float> &TriVtx3);
template double Area(const Vector3<double> &TriVtx1,
		const Vector3<double> &TriVtx2,
		const Vector3<double> &TriVtx3);
template float Area(const Vector2<float> &TriVtx1,
		const Vector2<float> &TriVtx2,
		const Vector2<float> &TriVtx3);
template double Area(const Vector2<double> &TriVtx1,
		const Vector2<double> &TriVtx2,
		const Vector2<double> &TriVtx3);

template Vector3<float> Normal(const Vector3<float> &TriVtx1,
		const Vector3<float> &TriVtx2,
		const Vector3<float> &TriVtx3,
		float *pArea);
template Vector3<double> Normal(const Vector3<double> &TriVtx1,
		const Vector3<double> &TriVtx2,
		const Vector3<double> &TriVtx3,
		double *pArea);

template Vector3<float>
InterpNormal(const Vector3<float> &vTriVtx1, const Vector3<float> &vTriVtx2,
		const Vector3<float> &vTriVtx3, const Vector3<float> &vTriNorm1,
		const Vector3<float> &vTriNorm2, const Vector3<float> &vTriNorm3,
		const Vector3<float> &vPointInTri);
template Vector3<double>
InterpNormal(const Vector3<double> &vTriVtx1, const Vector3<double> &vTriVtx2,
		const Vector3<double> &vTriVtx3, const Vector3<double> &vTriNorm1,
		const Vector3<double> &vTriNorm2, const Vector3<double> &vTriNorm3,
		const Vector3<double> &vPointInTri);

template void
StretchMetric1(const Vector3<float> &q1, const Vector3<float> &q2,
		const Vector3<float> &q3, const Vector2<float> &p1,
		const Vector2<float> &p2, const Vector2<float> &p3, float &MaxSV,
		float &MinSV, float &L2Norm, float &LInfNorm);
template void
StretchMetric1(const Vector3<double> &q1, const Vector3<double> &q2,
		const Vector3<double> &q3, const Vector2<double> &p1,
		const Vector2<double> &p2, const Vector2<double> &p3,
		double &MaxSV, double &MinSV, double &L2Norm, double &LInfNorm);

template void
StretchMetric3(const Vector3<float> &q1, const Vector3<float> &q2,
		const Vector3<float> &q3, const Vector3<float> &p1,
		const Vector3<float> &p2, const Vector3<float> &p3, float &MaxSV,
		float &MinSV, float &L2Norm, float &LInfNorm);
template void
StretchMetric3(const Vector3<double> &q1, const Vector3<double> &q2,
		const Vector3<double> &q3, const Vector3<double> &p1,
		const Vector3<double> &p2, const Vector3<double> &p3,
		double &MaxSV, double &MinSV, double &L2Norm, double &LInfNorm);

template bool IsObtuse(const Vector2<float> &v1,
		const Vector2<float> &v2,
		const Vector2<float> &v3);
template bool IsObtuse(const Vector2<double> &v1,
		const Vector2<double> &v2,
		const Vector2<double> &v3);
template bool IsObtuse(const Vector3<float> &v1,
		const Vector3<float> &v2,
		const Vector3<float> &v3);
template bool IsObtuse(const Vector3<double> &v1,
		const Vector3<double> &v2,
		const Vector3<double> &v3);
} // namespace g3

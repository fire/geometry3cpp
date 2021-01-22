#pragma once

#include <g3types.h>

#include <VectorUtil.h>

/*
* Frame3 represents an orthonormal 3D coordinate frame, ie 3 orthogonal axes and an origin point.
* Although it is possible to pack all that information into a Matrix4, or a Matrix3+Vector3, it
* is not very convenient to work with. Making the axes available explicitly allows for much
* clearer code.
*
* Frame3 also provides functions to:
*    - manipulate the frame by (eg) setting one axis to point in a particular direction
*    - transform vectors and frames to/from frame-relative coordinates
*    - (TODO: ray-intersection, project to plane, ...)
*
* Internally, Frame3 stores the inverse of the rotation matrix that would take the XYZ axes
* to the Frame axes. This is so that the axes are the rows of the matrix, and hence in
* row-major storage the 3 floats of each axis are contiguous.
*
*/
namespace g3 {

template <class Real>
class Frame3 {
public:
	//! origin point of reference frame
	Vector3<Real> Origin;

	//! transpose of matrix that takes canonical axes to frame axes. Rows are frame axes.
	Matrix3<Real> InverseRotation;

public:
	//! create a reference frame at a specific origin
	Frame3(const Vector3<Real> &Origin = Vector3<Real>::Zero()) :
			Origin(Origin), InverseRotation(Matrix3<Real>::Identity()) {}

	//! create an orthogonal frame at an origin with a given z axis vector
	Frame3(const Vector3<Real> &vOrigin,
			const Vector3<Real> &vNormalizedAxis, int nAxis,
			bool bUseFastPerpVectors) :
			Origin(vOrigin), InverseRotation(Matrix3<Real>::Identity()) {
		if (bUseFastPerpVectors && nAxis == 2) {
			Vector3<Real> kTangent0, kTangent1;
			g3::ComputePerpVectors(vNormalizedAxis, kTangent0, kTangent1, true);
			InverseRotation.row(0) = kTangent0;
			InverseRotation.row(1) = kTangent1;
			InverseRotation.row(2) = vNormalizedAxis;
		} else {
			AlignAxis(nAxis, vNormalizedAxis);
		}
	}

	//! copy a reference frame
	Frame3(const Frame3 &copy) :
			Origin(copy.Origin), InverseRotation(copy.InverseRotation) {}

	~Frame3();

	//! get an axis of the frame
	Vector3<Real> Axis(unsigned int nAxis) const {
		return InverseRotation.row(nAxis);
	}

	//! access the axes of the frame
	Vector3<Real> X() const {
		return InverseRotation.row(0);
	}
	Vector3<Real> Y() const {
		return InverseRotation.row(1);
	}
	Vector3<Real> Z() const {
		return InverseRotation.row(2);
	}

	//! matrix that rotates canonical/unit XYZ axes to frame axes
	Matrix3<Real> GetRotation() const {
		return InverseRotation.transpose();
	}

	//! treat v as begin in frame-relative coords, transform to absolute/world coords
	void SetToWorldCoords(Vector3<Real> &v) const {
		v = Origin + InverseRotation.transpose() * v;
	}
	//! treat v as begin in frame-relative coords, transform to absolute/world coords
	Vector3<Real> ToWorldCoords(const Vector3<Real> &v) const {
		return Origin + InverseRotation.transpose() * v;
	}

	//! treat v as begin in absolute/world coords, transform to frame-relative coords
	void SetToAxisCoords(Vector3<Real> &v) const {
		v = InverseRotation * (v - Origin);
	}
	//! treat v as begin in absolute/world coords, transform to frame-relative coords
	Vector3<Real> ToAxisCoords(const Vector3<Real> &v) const {
		return InverseRotation * (v - Origin);
	}

	//! treat f as being relative to this frame, transform to absolute/world frame
	Frame3<Real> ToWorldCoords(const Frame3<Real> &f, bool bPreserveOrigin = false) const {
		Vector3<Real> x = InverseRotation.transpose() * f.X();
		Vector3<Real> y = InverseRotation.transpose() * f.Y();
		Vector3<Real> z = InverseRotation.transpose() * f.Z();
		Frame3<Real> l((bPreserveOrigin) ? f.Origin : ToWorldCoords(f.Origin));
		l.SetFrame(x, y, z);
		return l;
	}

	//! treat f as being an absolute/world frame, transform to be relative to this frame
	Frame3<Real> ToAxisCoords(const Frame3<Real> &f, bool bPreserveOrigin = false) const {
		Vector3<Real> x = InverseRotation * f.X();
		Vector3<Real> y = InverseRotation * f.Y();
		Vector3<Real> z = InverseRotation * f.Z();
		Frame3<Real> w((bPreserveOrigin) ? f.Origin : ToAxisCoords(f.Origin));
		w.SetFrame(x, y, z);
		return w;
	}

	//! translate the reference frame origin
	void Translate(const Vector3<Real> &vTranslate, bool bRelative) {
		if (bRelative)
			Origin += vTranslate;
		else
			Origin = vTranslate;
	}
	//! rotate the reference frame
	void Rotate(const Matrix3<Real> &mRotation, bool bReNormalize) {
		// [RMS] does not seem to work as expected...??
		Vector3<Real> rowX(mRotation * InverseRotation.row(0).transpose());
		Vector3<Real> rowY(mRotation * InverseRotation.row(1).transpose());
		Vector3<Real> rowZ(mRotation * InverseRotation.row(2).transpose());

		if (bReNormalize) {
			Vector3<Real> vCrossY = rowZ.cross(rowX);
			vCrossY.normalize();
			rowY = (vCrossY.dot(rowY) < 0) ? -vCrossY : vCrossY;
			Vector3<Real> vCrossX = rowZ.cross(rowY);
			vCrossX.normalize();
			rowX = (vCrossX.dot(rowX) < 0) ? -vCrossX : vCrossX;
			Vector3<Real> vCrossZ = rowX.cross(rowY);
			vCrossZ.normalize();
			rowZ = (vCrossZ.dot(rowZ) < 0) ? -vCrossZ : vCrossZ;
		}

		InverseRotation.row(0) = rowX;
		InverseRotation.row(1) = rowY;
		InverseRotation.row(2) = rowZ;
	}

	//! make axes perpendicular. can "preserve" one axis by setting nPreserveAxis (0=X,1=Y,2=Z)

	void ReNormalize(int nPreserveAxis) {
		Vector3<Real> rowX(InverseRotation.row(0));
		Vector3<Real> rowY(InverseRotation.row(1));
		Vector3<Real> rowZ(InverseRotation.row(2));

		switch (nPreserveAxis) {
			default:
			case -1: {
				Vector3<Real> vCrossY = rowZ.cross(rowX);
				vCrossY.normalize();
				rowY = (vCrossY.dot(rowY) < 0) ? -vCrossY : vCrossY;
				Vector3<Real> vCrossX = rowZ.cross(rowY);
				vCrossX.normalize();
				rowX = (vCrossX.dot(rowX) < 0) ? -vCrossX : vCrossX;
				Vector3<Real> vCrossZ = rowX.cross(rowY);
				vCrossZ.normalize();
				rowZ = (vCrossZ.dot(rowZ) < 0) ? -vCrossZ : vCrossZ;
			} break;

			case 0: {
				rowX.normalize();
				Vector3<Real> vCrossY(rowX.cross(rowZ));
				vCrossY.normalize();
				rowY = (vCrossY.dot(rowY) < 0) ? -vCrossY : vCrossY;
				Vector3<Real> vCrossZ(rowX.cross(rowY));
				vCrossZ.normalize();
				rowZ = (vCrossZ.dot(rowZ) < 0) ? -vCrossZ : vCrossZ;
			} break;

			case 1: {
				rowY.normalize();
				Vector3<Real> vCrossX(rowY.cross(rowZ));
				vCrossX.normalize();
				rowX = (vCrossX.dot(rowX) < 0) ? -vCrossX : vCrossX;
				Vector3<Real> vCrossZ(rowX.cross(rowY));
				vCrossZ.normalize();
				rowZ = (vCrossZ.dot(rowZ) < 0) ? -vCrossZ : vCrossZ;
			} break;

			case 2: {
				rowZ.normalize();
				Vector3<Real> vCrossX(rowY.cross(rowZ));
				vCrossX.normalize();
				rowX = (vCrossX.dot(rowX) < 0) ? -vCrossX : vCrossX;
				Vector3<Real> vCrossY(rowX.cross(rowZ));
				vCrossY.normalize();
				rowY = (vCrossY.dot(rowY) < 0) ? -vCrossY : vCrossY;
			} break;
		}

		InverseRotation.row(0) = rowX;
		InverseRotation.row(1) = rowY;
		InverseRotation.row(2) = rowZ;
	}

	//! intersect ray with plane defined by one of the axis vectors
	bool IntersectRay(const Vector3<Real> &vRayOrigin,
			const Vector3<Real> &vRayDirection,
			Vector3<Real> &vRayHit, int nPlaneNormalAxis) {
		Vector3<Real> N = Axis(nPlaneNormalAxis);
		Real d = -(Origin.dot(N));
		Real fDenom = vRayDirection.dot(N);
		if (fDenom < Math<Real>::ZERO_TOLERANCE)
			return false;

		Real t = -(vRayOrigin.dot(N) + d) / fDenom;
		vRayHit = vRayOrigin + t * vRayDirection;
		return true;
	}
	Vector3<Real> IntersectRay(const Vector3<Real> &vRayOrigin,
			const Vector3<Real> &vRayDirection,
			int nPlaneNormalAxis) {
		Vector3<Real> vHit = Vector3<Real>::Zero();
		bool bOK = IntersectRay(vRayOrigin, vRayDirection, vHit, nPlaneNormalAxis);
		// [TODO] return invalid vector here instead of ZERO
		return (bOK) ? vHit : Vector3<Real>::Zero();
	}

	//! create frame at an origin with a given x/y/z vectors
	Frame3(const Vector3<Real> &vOrigin, const Vector3<Real> &vXAxis,
			const Vector3<Real> &vYAxis, const Vector3<Real> &vZAxis) :
			Origin(vOrigin), InverseRotation(Matrix3<Real>::Identity()) {
		SetFrame(vXAxis, vYAxis, vZAxis);
	}

	//! create frame with a given rotation
	Frame3(const Matrix3<Real> &mRotation) :
			Origin(Vector3<Real>::Zero()), InverseRotation(mRotation.transpose()) {}

	//! create frame at an origin with a given rotation
	Frame3(const Vector3<Real> &vOrigin,
			const Matrix3<Real> &mRotation) :
			Origin(vOrigin), InverseRotation(mRotation.transpose()) {}

	//! allow external setting of matrix...
	void SetFrame(const Matrix3<Real> &matRotate) {
		InverseRotation = matRotate.transpose();
	}

	//! set the reference frame axes
	void SetFrame(const Vector3<Real> &vAxisX,
			const Vector3<Real> &vAxisY,
			const Vector3<Real> &vAxisZ) {
		InverseRotation.row(0) = vAxisX.normalized();
		InverseRotation.row(1) = vAxisY.normalized();
		InverseRotation.row(2) = vAxisZ.normalized();
	}

	//! rotate selected axis of this frame into toAxis
	void AlignAxis(int nAxis, const Vector3<Real> &toAxis,
			bool bNormalized) {
		Matrix3<Real> matAlign;
		g3::ComputeAlignAxisMatrix(
				Axis(nAxis), (bNormalized) ? toAxis : toAxis.normalized(), matAlign);
		Rotate(matAlign, false);
	}

	//! compute matrix that rotates this frame into toFrame
	void ComputeAlignmentMatrix(const Frame3<Real> &toFrame,
			Matrix3<Real> &matRotate) {
		// align vCurFrame.Z() with vDestFrame.Z()
		Matrix3<Real> matAlignZ;
		ComputeAlignAxisMatrix(this->Z(), toFrame.Z(), matAlignZ);
		Frame3<Real> vCopy(*this);
		vCopy.Rotate(matAlignZ);

		// compute rotation angle around vDestFrame.Z()
		Vector3<Real> vX1(toFrame.X());
		Vector3<Real> vX2(vCopy.X());
		Matrix3<Real> matAlignX;
		ComputeAlignAxisMatrix(vX2, vX1, matAlignX);

		matRotate = matAlignX * matAlignZ;

		vCopy = Frame3<Real>(*this);
		vCopy.Rotate(matRotate);
		Real fDotX = vCopy.X().dot(toFrame.X());
		Real fDotY = vCopy.Y().dot(toFrame.Y());
		Real fDotZ = vCopy.Z().dot(toFrame.Z());

		// check if Y is flipped - if it is, flip Y...
		// if ( vCopy.Y().dot( toFrame.Y() ) < 0 ) {
		//	Matrix3<Real> matFlip( 1,0,0, 0,-1,0, 0,0,1 );
		//	matRotate = matRotate * matFlip;
		//}
		//// check if Z is flipped - if it is, flip Z...
		if (vCopy.Z().dot(toFrame.Z()) < 0) {
			Matrix3<Real> matFlip;
			matFlip.row(0) = Vector3<Real>(1, 0, 0);
			matFlip.row(1) = Vector3<Real>(0, 1, 0);
			matFlip.row(2) = Vector3<Real>(0, 0, -1);
			matRotate = matRotate * matFlip;
		}

		// [RMS] does this do anything??
		vCopy = Frame3<Real>(*this);
		vCopy.Rotate(matRotate);
		fDotX = vCopy.X().dot(toFrame.X());
		fDotY = vCopy.Y().dot(toFrame.Y());
		fDotZ = vCopy.Z().dot(toFrame.Z());
	}
};
typedef Frame3<float> Frame3f;
typedef Frame3<double> Frame3d;
} // namespace g3
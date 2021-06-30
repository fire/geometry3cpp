#pragma once

/*
 * [RMS] this file defines and/or pre-declares many core types for the geometry3 library
 *   Note that order-of-include matters, because of nested includes, so some of the .h files
 *   are only included further down...
 */

#include <Eigen/Core>
#include <Eigen/Geometry> // required for MatrixBase.cross()  !!

#include <limits>
#include <memory>

#include <g3Config.h>

#include "DCurve3.h"
#include "src/mesh/DMesh3.h"

std::shared_ptr<DMesh> DMesh3Ptr;
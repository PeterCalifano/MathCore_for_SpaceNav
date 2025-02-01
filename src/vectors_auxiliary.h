/**
 * @file vectors_auxiliary.h
 * @author PeterC (petercalifano.gs@gmail.com)
 * @brief
 * @note Copy pasted from math_utils_for_CV.h in VO_directionOfMotion project
 * @version 0.1
 * @date 2025-01-13
 */
#pragma once

#include <Eigen/Dense>
#include <global_types.h>
#include <iostream>
#include <vector>

/**
 * @brief Function converting a Vector3 or any Eigen equivalent representation into a skew-symmetric matrix. Type is Derived::Scalar.
 *
 * @tparam Derived
 * @param VectCross
 * @return Matrix33<typename Derived::Scalar>
 */

template <typename T>
using Matrix33 = Eigen::Matrix<T, 3, 3>;


template <typename Derived>
Matrix33<typename Derived::Scalar> SkewSymmetricCross(const Eigen::MatrixBase<Derived> &VectCross)
{
    static_assert(Derived::RowsAtCompileTime == 3 && Derived::ColsAtCompileTime == 1,
                  "SkewSymmetricCross expects a 3x1 vector.");

    return (Matrix33<typename Derived::Scalar>() << 0, -VectCross(2), VectCross(1),
            VectCross(2), 0, -VectCross(0),
            -VectCross(1), VectCross(0), 0)
        .finished();
};
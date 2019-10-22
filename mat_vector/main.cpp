#define CATCH_CONFIG_MAIN
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "catch.h"

#define loop(a) for(int i = 0; i < a; i++)
#define loop_(a) for(int j = 0; j < a; j++)

const double e = 0.000000000001;

using namespace mat_vec;

TEST_CASE("Vector::Vector(size_t size, double value)") {
    Vector vector(5, 5);
    REQUIRE(vector.size() == 5);
    loop(5) {
        REQUIRE(vector[i] - 5.0 < e);
    }
}

TEST_CASE("Vector::Vector(const Vector &src)") {
    Vector vector_src(5, 5);
    Vector vector(vector_src);
    REQUIRE(vector.size() == 5);
    loop(5) {
        REQUIRE(vector[i] - 5.0 < e);
    }
}

TEST_CASE("Vector &Vector::operator=(const Vector &rhs)") {
    Vector vector_src(5, 5);
    Vector vector = vector_src;
    REQUIRE(vector.size() == 5);
    loop(5) {
        REQUIRE(vector[i] - 5.0 < e);
    }
}

TEST_CASE("size_t Vector::size() const") {
    Vector vector(7);
    REQUIRE(vector.size() == 7);
}

TEST_CASE("double Vector::operator[](size_t n) const") {
    Vector vector(5, 5);
    REQUIRE(vector[1] - 5.0 < e);
}

TEST_CASE("double &Vector::operator[](size_t n)") {
    Vector vector(5, 5);
    vector[1] = 7;
    REQUIRE(vector[1] - 7 < e);
}

TEST_CASE("double Vector::norm() const") {
    Vector vector(4, 5);
    REQUIRE(vector.norm() - 10 < e);
}

TEST_CASE("Vector Vector::normalized() const") {
    Vector vector(4, 5);
    Vector vector_src(vector);
    vector = vector_src.normalized();
    loop(4) {
        REQUIRE(vector[i] - 0.5 < e);
    }
}

TEST_CASE("void Vector::normalize()") {
    Vector vector(4, 5);
    vector.normalize();
    loop(4) {
        REQUIRE(vector[i] - 0.5 < e);
    }
}

TEST_CASE("Vector Vector::operator+(const Vector &rhs) const") {
    Vector vector(4, 5);
    Vector _vector(4, 7);
    Vector __vector(_vector);
    __vector = _vector + vector;
    loop(4) {
        REQUIRE(__vector[i] - 12 < e);
    }
}

TEST_CASE("Vector &Vector::operator+=(const Vector &rhs)") {
    Vector vector(4, 5);
    Vector _vector(4, 7);
    _vector += vector;
    loop(4) {
        REQUIRE(_vector[i] - 12 < e);
    }
}

TEST_CASE("Vector Vector::operator-(const Vector &rhs) const") {
    Vector vector(4, 5);
    Vector _vector(4, 7);
    Vector __vector(_vector);
    __vector = vector - _vector;
    loop(4) {
        REQUIRE(__vector[i] - 2 < e);
    }
}
TEST_CASE("Vector &Vector::operator-=(const Vector &rhs)") {
    Vector vector(4, 5);
    Vector _vector(4, 7);
    vector -= _vector;
    loop(4) {
        REQUIRE(vector[i] - 2 < e);
    }
}

TEST_CASE("Vector Vector::operator^(const Vector &rhs) const") {
    Vector vector(4, 5);
    Vector _vector(4, 2);
    Vector __vector(_vector);
    __vector = vector ^ _vector;
    loop(4) {
        REQUIRE(__vector[i] - 10 < e);
    }
}

TEST_CASE("Vector &Vector::operator^=(const Vector &rhs)") {
    Vector vector(4, 5);
    Vector _vector(4, 2);
    vector ^= _vector;
    loop(4) {
        REQUIRE(_vector[i] - 10 < e);
    }
}

TEST_CASE("double Vector::operator*(const Vector &rhs) const") {
    Vector vector(4, 5);
    Vector _vector(4, 2);
    double result = vector * _vector;
    REQUIRE(result - 40 < e);
}

TEST_CASE("Vector Vector::operator*(double k) const") {
    Vector vector(4, 5);
    Vector _vector(4, 2);
    vector = _vector * 5;
    loop(4) {
        REQUIRE(vector[i] - 10 < e);
    }
}

TEST_CASE("Vector &Vector::operator*=(double k)") {
    Vector vector(4, 2);
    vector *= 5;
    loop(4) {
        REQUIRE(vector[i] - 10 < e);
    }
}

TEST_CASE("Vector mat_vec::operator*(double k, const Vector &v)") {
    Vector vector(4, 5);
    Vector _vector(4, 2);
    vector = 5 * _vector;
    loop(4) {
        REQUIRE(vector[i] - 10 < e);
    }
}

TEST_CASE("Vector Vector::operator/(double k) const") {
    Vector vector(4, 6);
    Vector _vector(4, 2);
    _vector = vector / 2;
    loop(4) {
        REQUIRE(_vector[i] - 3 < e);
    }
}

TEST_CASE("Vector &Vector::operator/=(double k)") {
    Vector vector(4, 6);
    vector /= 2;
    loop(4) {
        REQUIRE(vector[i] - 3 < e);
    }
}

TEST_CASE("Vector Vector::operator*(const Matrix &mat) const") {
    Matrix matrix(2, 3, 0);
    Vector vector(2);
    loop(2) {
        loop_(3) {
            matrix[i * 3 + j] = i * 3 + j + 1;
        }
    }
    loop(2) {
        vector[i] = i + 1;
    }
    Vector _vector = vector * matrix;
    REQUIRE(_vector.size() == 3);
    REQUIRE(_vector[0] - 9 < e);
    REQUIRE(_vector[1] - 12 < e);
    REQUIRE(_vector[2] - 15 < e);
}

TEST_CASE("Vector &Vector::operator*=(const Matrix &mat)") {
    Matrix matrix(2, 3, 0);
    Vector vector(2);
    loop(2) {
        loop_(3) {
            matrix[i * 3 + j] = i * 3 + j + 1;
        }
    }
    loop(2) {
        vector[i] = i + 1;
    }
    vector *= matrix;
    REQUIRE(vector.size() == 3);
    REQUIRE(vector[0] - 9 < e);
    REQUIRE(vector[1] - 12 < e);
    REQUIRE(vector[2] - 15 < e);
}

TEST_CASE("bool Vector::operator==(const Vector &rhs) const") {
    Vector vector(4, 5);
    Vector _vector(4, 3);
    Vector rhs(4, 5);
    REQUIRE(rhs == vector);
    REQUIRE(!(_vector == vector));
}

TEST_CASE("bool Vector::operator!=(const Vector &rhs) const") {
    Vector vector(4, 5);
    Vector _vector(4, 5);
    Vector rhs(4, 4);
    REQUIRE(rhs != vector);
    REQUIRE(!(_vector != vector));
}

//
TEST_CASE("Matrix::Matrix(size_t size, double value)") {
    Matrix matrix(3, 5);
    REQUIRE(matrix.shape().first == 3);
    REQUIRE(matrix.shape().second == 3);
    loop(3) {
        loop_(3) {
            REQUIRE(matrix.get(i, j) - 5 < e);
        }
    }
}

TEST_CASE("Matrix Matrix::eye(size_t size)") {
    Matrix matrix(3);
    matrix = matrix.eye(3);
    loop(3) {
        loop_(3) {
            if (i == j) {
                REQUIRE(matrix.get(i, j) - 1 < e);
            } else {
                REQUIRE(matrix.get(i, j) - 0 < e);
            }
        }
    }
}

TEST_CASE("Matrix::Matrix(size_t rows, size_t cols, double value)") {
    Matrix matrix(2, 4, 3);
    loop(2) {
        loop_(4) {
            REQUIRE(matrix.get(i, j) - 3 < e);
        }
    }
}

TEST_CASE("Matrix::Matrix(const Matrix &src)") {
    Matrix matrix(5, 3);
    Matrix _matrix(matrix);
    loop(5) {
        loop_(5) {
            REQUIRE(_matrix.get(i, j) - 3 < e);
        }
    }
}

TEST_CASE("Matrix &Matrix::operator=(const Matrix &rhs)") {
    Matrix matrix(5, 3);
    Matrix _matrix = matrix;
    loop(5) {
        loop_(5) {
            REQUIRE(_matrix.get(i, j) - 3 < e);
        }
    }
}

TEST_CASE("void Matrix::reshape(size_t rows, size_t cols)\nstd::pair<size_t, size_t> Matrix::shape() const") {
    Matrix matrix(2, 4, 3);
    REQUIRE(matrix.shape().first == 2);
    REQUIRE(matrix.shape().second == 4);
    loop(matrix.shape().first) {
        loop_(matrix.shape().second) {
            matrix[i * matrix.shape().second + j] = i * matrix.shape().second + j + 1;
        }
    }
    matrix.reshape(4, 2);
    REQUIRE(matrix.shape().first == 4);
    REQUIRE(matrix.shape().second == 2);
    loop(matrix.shape().first) {
        loop_(matrix.shape().second) {
            REQUIRE(matrix.get(i, j) - (i * matrix.shape().second + j + 1) < e);
        }
    }
}

TEST_CASE("Matrix Matrix::operator+(const Matrix &rhs) const") {
    Matrix matrix(4, 5);
    Matrix _matrix(4, 7);
    Matrix __matrix(_matrix);
    __matrix = _matrix + matrix;
    loop(4) {
        loop_(4) {
            REQUIRE(__matrix.get(i, j) - 12 < e);
        }
    }
}

TEST_CASE("Matrix &Matrix::operator+=(const Matrix &rhs)") {
    Matrix matrix(4, 5);
    Matrix _matrix(4, 7);
    _matrix += matrix;
    loop(4) {
        loop_(4) {
            REQUIRE(_matrix.get(i, j) - 12 < e);
        }
    }
}

TEST_CASE("Matrix Matrix::operator-(const Matrix &rhs) const") {
    Matrix matrix(4, 5);
    Matrix _matrix(4, 7);
    Matrix __matrix(_matrix);
    __matrix = _matrix - matrix;
    loop(4) {
        loop_(4) {
            REQUIRE(__matrix.get(i, j) - 2 < e);
        }
    }
}

TEST_CASE("Matrix &Matrix::operator-=(const Matrix &rhs)") {
    Matrix matrix(4, 5);
    Matrix _matrix(4, 7);
    _matrix -= matrix;
    loop(4) {
        loop_(4) {
            REQUIRE(_matrix.get(i, j) - 2 < e);
        }
    }
}

TEST_CASE("Matrix Matrix::operator*(const Matrix &rhs) const") {
    Matrix matrix(2, 4, 0);
    matrix[0] = 1; matrix[1] = 2; matrix[2] = 3; matrix[3] = 4;
    matrix[4] = 5; matrix[5] = 6; matrix[6] = 7; matrix[7] = 8;
    Matrix _matrix(4, 3, 0);
    _matrix[0] = 1; _matrix[1] = 2; _matrix[2] = 3;
    _matrix[3] = 4; _matrix[4] = 5; _matrix[5] = 6;
    _matrix[6] = 7; _matrix[7] = 8; _matrix[8] = 9;
    _matrix[9] = 10; _matrix[10] = 11; _matrix[11] = 12;
    Matrix __matrix = _matrix * matrix;
    REQUIRE(__matrix.get(0, 0) - 70 < e);
    REQUIRE(__matrix.get(0, 1) - 80 < e);
    REQUIRE(__matrix.get(0, 2) - 90 < e);
    REQUIRE(__matrix.get(1, 0) - 158 < e);
    REQUIRE(__matrix.get(1, 1) - 184 < e);
    REQUIRE(__matrix.get(1, 2) - 210 < e);
}

TEST_CASE("Matrix &Matrix::operator*=(const Matrix &rhs)") {
    Matrix matrix(2, 4, 0);
    matrix[0] = 1; matrix[1] = 2; matrix[2] = 3; matrix[3] = 4;
    matrix[4] = 5; matrix[5] = 6; matrix[6] = 7; matrix[7] = 8;
    Matrix _matrix(4, 3, 0);
    _matrix[0] = 1; _matrix[1] = 2; _matrix[2] = 3;
    _matrix[3] = 4; _matrix[4] = 5; _matrix[5] = 6;
    _matrix[6] = 7; _matrix[7] = 8; _matrix[8] = 9;
    _matrix[9] = 10; _matrix[10] = 11; _matrix[11] = 12;
    _matrix *= matrix;
    REQUIRE(_matrix.get(0, 0) - 70 < e);
    REQUIRE(_matrix.get(0, 1) - 80 < e);
    REQUIRE(_matrix.get(0, 2) - 90 < e);
    REQUIRE(_matrix.get(1, 0) - 158 < e);
    REQUIRE(_matrix.get(1, 1) - 184 < e);
    REQUIRE(_matrix.get(1, 2) - 210 < e);
}

TEST_CASE("Matrix Matrix::operator*(double k) const") {
    Matrix matrix(4, 6);
    Matrix _matrix(4);
    _matrix = matrix * 2;
    loop(4) {
        loop_(4) {
            REQUIRE(_matrix.get(i, j) - 12 < e);
        }
    }
}

TEST_CASE("Matrix &Matrix::operator*=(double k)") {
    Matrix matrix(4, 6);
    matrix *= 2;
    loop(4) {
        loop_(4) {
            REQUIRE(matrix.get(i, j) - 12 < e);
        }
    }
}

TEST_CASE("Matrix Matrix::operator/(double k) const") {
    Matrix matrix(4, 6);
    Matrix _matrix(4);
    _matrix = matrix / 2;
    loop(4) {
        loop_(4) {
            REQUIRE(_matrix.get(i, j) - 3 < e);
        }
    }
}

TEST_CASE("Matrix &Matrix::operator/=(double k)") {
    Matrix matrix(4, 6);
    matrix /= 2;
    loop(4) {
        loop_(4) {
            REQUIRE(matrix.get(i, j) - 3 < e);
        }
    }
}

TEST_CASE("Matrix Matrix::transposed() const") {
    Matrix matrix(3);
    matrix[0] = 1; matrix[1] = 2; matrix[2] = 4;
    matrix[3] = 5; matrix[4] = 7; matrix[5] = 8;
    matrix[6] = 10; matrix[7] = 11; matrix[8] = 13;
    Matrix _matrix = matrix.transposed();
    REQUIRE(_matrix.get(0, 0) - 1 < e);
    REQUIRE(_matrix.get(0, 1) - 5 < e);
    REQUIRE(_matrix.get(0, 2) - 10 < e);
    REQUIRE(_matrix.get(1, 0) - 2 < e);
    REQUIRE(_matrix.get(1, 1) - 7 < e);
    REQUIRE(_matrix.get(1, 2) - 11 < e);
    REQUIRE(_matrix.get(2, 0) - 4 < e);
    REQUIRE(_matrix.get(2, 1) - 8 < e);
    REQUIRE(_matrix.get(2, 2) - 13 < e);
}

TEST_CASE("void Matrix::transpose()") {
    Matrix matrix(3);
    matrix[0] = 1; matrix[1] = 2; matrix[2] = 4;
    matrix[3] = 5; matrix[4] = 7; matrix[5] = 8;
    matrix[6] = 10; matrix[7] = 11; matrix[8] = 13;
    matrix.transpose();
    REQUIRE(matrix.get(0, 0) - 1 < e);
    REQUIRE(matrix.get(0, 1) - 5 < e);
    REQUIRE(matrix.get(0, 2) - 10 < e);
    REQUIRE(matrix.get(1, 0) - 2 < e);
    REQUIRE(matrix.get(1, 1) - 7 < e);
    REQUIRE(matrix.get(1, 2) - 11 < e);
    REQUIRE(matrix.get(2, 0) - 4 < e);
    REQUIRE(matrix.get(2, 1) - 8 < e);
    REQUIRE(matrix.get(2, 2) - 13 < e);
}

TEST_CASE("double Matrix::det() const") {
    Matrix matrix(3);
    matrix[0] = 1; matrix[1] = 2; matrix[2] = 4;
    matrix[3] = 5; matrix[4] = 7; matrix[5] = 8;
    matrix[6] = 10; matrix[7] = 11; matrix[8] = 13;
    double det = matrix.det();
    REQUIRE(det - (-27) < e);
}

TEST_CASE("Matrix Matrix::inv() const") {
    double eps = 0.000001;
    Matrix matrix(3);
    matrix[0] = 1; matrix[1] = 2; matrix[2] = 4;
    matrix[3] = 5; matrix[4] = 7; matrix[5] = 8;
    matrix[6] = 10; matrix[7] = 11; matrix[8] = 13;
    matrix = matrix.inv();
    REQUIRE(matrix.get(0, 0) - (-0.111111) < eps);
    REQUIRE(matrix.get(0, 1) - (-0.666667) < eps);
    REQUIRE(matrix.get(0, 2) - (0.444444) < eps);
    REQUIRE(matrix.get(1, 0) - (-0.555556) < eps);
    REQUIRE(matrix.get(1, 1) - (1) < eps);
    REQUIRE(matrix.get(1, 2) - (-0.444444) < eps);
    REQUIRE(matrix.get(2, 0) - (0.555556) < eps);
    REQUIRE(matrix.get(2, 1) - (-0.333333) < eps);
    REQUIRE(matrix.get(2, 2) - (0.111111) < eps);
}

TEST_CASE("Vector Matrix::operator*(const Vector &vec) const") {
    Matrix matrix(3, 2, 0);
    Vector vector(2);
    loop(3) {
        loop_(2) {
            matrix[i * 2 + j] = i * 2 + j + 1;
        }
    }
    loop(2) {
        vector[i] = i + 1;
    }
    Vector _vector = matrix * vector;
    REQUIRE(_vector.size() == 3);
    REQUIRE(_vector[0] - 5 < e);
    REQUIRE(_vector[1] - 11 < e);
    REQUIRE(_vector[2] - 17 < e);
}

TEST_CASE("bool Matrix::operator==(const Matrix &rhs) const") {
    Matrix matrix(5, 5);
    Matrix _matrix(5, 5);
    Matrix __matrix(5, 3);
    REQUIRE(matrix == _matrix);
    REQUIRE(!(matrix == __matrix));
}

TEST_CASE("bool Matrix::operator!=(const Matrix &rhs) const") {
    Matrix matrix(5, 5);
    Matrix _matrix(5, 3);
    Matrix __matrix(5, 3);
    REQUIRE(matrix != _matrix);
    REQUIRE(!(__matrix != _matrix));
}


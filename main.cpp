// Andrey Alexeev
// a.alexeev@innopolis.university
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

class Matrix {
public:
    int n, m;
    vector<vector<double>> values;

    Matrix(int nn, int mm, const vector<vector<double>> &arr) {
        n = nn;
        m = mm;
        values = arr;
    }

    Matrix *transpose() {
        vector<vector<double>> temp(m, vector<double>(n));
        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) {
                temp[j][i] = values[i][j];
            }
        }
        return new Matrix(m, n, temp);
    }

    int determinant() {
        double d = 1;
        for (int i = 0; i < n; i++) {
            d *= (values[i][i]);
        }
        return ceil(d * 100.0) / 100.0;
    }

    bool operator=(const Matrix &m2) {
        values = m2.values;
        n = m2.n;
        m = m2.m;
        return true;
    }

    friend istream &operator>>(istream &in, Matrix &mat) {
        in >> mat.n;
        mat.m = mat.n;
        mat.values.resize(mat.n, vector<double>(mat.m));
        for (int i = 0; i < mat.n; i++) {
            for (int j = 0; j < mat.m; j++) {
                in >> mat.values[i][j];
            }
        }
        return in;
    }

    friend ostream &operator<<(ostream &out, Matrix &mat) {
        for (int i = 0; i < mat.n; i++) {
            for (int j = 0; j < mat.m; j++) {
                if (j == mat.m - 1) {
                    if (mat.values[i][j] < 1e-10 && mat.values[i][j] > -1e-10) {
                        out << fixed << setprecision(4) << 0.00;
                    } else {
                        out << fixed << setprecision(4) << mat.values[i][j];
                    }
                } else {
                    if (mat.values[i][j] < 1e-10 && mat.values[i][j] > -1e-10) {
                        out << fixed << setprecision(4) << 0.00 << " ";
                    } else {
                        out << fixed << setprecision(4) << mat.values[i][j] << " ";
                    }
                }
            }
            out << endl;
        }
        return out;
    }
};

Matrix *operator+(Matrix m1, Matrix m2) {
    if (m1.n == m2.n && m1.m == m2.m) {
        vector<vector<double>> newVector(m1.n, vector<double>(m1.m));
        for (int i = 0; i < m1.n; i++) {
            for (int j = 0; j < m1.m; j++) {
                newVector[i][j] = m1.values[i][j] + m2.values[i][j];
            }
        }
        return new Matrix(m1.n, m1.m, newVector);
    } else {
        vector<vector<double>> emptyVector(0, vector<double>(0));
        cout << "Error: the dimensional problem occurred" << endl;
        return new Matrix(0, 0, emptyVector);
    }
}

Matrix *operator-(Matrix m1, Matrix m2) {
    if (m1.n == m2.n && m1.m == m2.m) {
        vector<vector<double>> newVector(m1.n, vector<double>(m1.m));
        for (int i = 0; i < m1.n; i++) {
            for (int j = 0; j < m1.m; j++) {
                newVector[i][j] = m1.values[i][j] - m2.values[i][j];
            }
        }
        return new Matrix(m1.n, m1.m, newVector);
    } else {
        vector<vector<double>> emptyVector(0, vector<double>(0));
        cout << "Error: the dimensional problem occurred" << endl;
        return new Matrix(0, 0, emptyVector);
    }
}

Matrix *operator*(Matrix m1, Matrix m2) {
    if (m1.m == m2.n) {
        int newN = m1.n;
        int newM = m2.m;
        vector<vector<double>> newVector(newN, vector<double>(newM));
        for (int i = 0; i < m1.n; i++) {
            for (int j = 0; j < m2.m; j++) {
                for (int k = 0; k < m1.m; k++) {
                    newVector[i][j] += m1.values[i][k] * m2.values[k][j];
                }
            }
        }
        return new Matrix(newN, newM, newVector);
    } else {
        vector<vector<double>> emptyVector(0, vector<double>(0));
        cout << "Error: the dimensional problem occurred" << endl;
        return new Matrix(0, 0, emptyVector);
    }
}


class Identity : public Matrix {
public:
    Identity(int nn, int mm, const vector<vector<double>> &arr) : Matrix(nn, mm, arr) {
        for (int i = 0; i < nn; i++) {
            vector<double> temp;
            for (int j = 0; j < mm; j++) {
                if (i == j) {
                    temp.push_back(1);
                } else {
                    temp.push_back(0);
                }
            }
            values.push_back(temp);
        }
    }
};

class Elimination : public Matrix {
public:
    Elimination(int nn, int mm, const vector<vector<double>> &arr, Matrix *host, Matrix *identity, int row, int col)
            : Matrix(nn, mm, arr) {
        row--;
        col--;
        values = identity->values;
        double pivot = host->values[row][col] / host->values[col][col];
        for (int i = 0; i < nn; i++) {
            values[row][i] -= values[col][i] * pivot;
        }
    }
};

class Permutation : public Matrix {
public:
    Permutation(int nn, int mm, const vector<vector<double>> &arr, Matrix *identity, int row, int col)
            : Matrix(nn, mm, arr) {
        row--;
        col--;
        values = identity->values;
        vector<double> temp;
        temp = values[row];
        values[row] = values[col];
        values[col] = temp;
    }
};



int main() {
    FILE* pipe = _popen("C:\\\\gnuplot\\\\bin\\\\gnuplot -persist", "w");

    int n;
    cin >> n;
    vector<vector<double>> b;
    double saveA[n];
    double saveB[n];
    vector<vector<double>> temp;
    for (int i = 0; i < n; i++) {
        vector<double> f;
        temp.push_back(f);
        temp[i].push_back(1);
        double Atemp, Btemp;
        cin >> Atemp >> Btemp;
        saveA[i] = Atemp;
        saveB[i] = Btemp;
        temp[i].push_back(Atemp);
        vector<double> bb;
        b.push_back(bb);
        b[i].push_back(Btemp);
    }
    int degree;
    cin >> degree;
    for (int i = 1; i < degree; i++) {
        for (int j = 0; j < n; j++) {
            temp[j].push_back(pow(temp[j][1], i + 1));
        }
    }
    Matrix A(n, degree + 1, temp);
    Matrix B(n, 1, b);
    cout << "A:" << endl;
    cout << A;
    cout << "A_T*A:" << endl;
    Matrix A_T = *A.transpose();
    Matrix A_TA = *(A_T * A);
    cout << A_TA;
    cout << "(A_T*A)^-1:" << endl;
    Identity ii(A_TA.n, A_TA.m, vector<vector<double>>());
    Matrix *A_TAREADY;
    A_TAREADY = &ii;


    for (int i = 0; i < A_TA.n; i++) {
        double nowPivot = A_TA.values[i][i];
        int index = -1;
        for (int j = i + 1; j < A_TA.n; j++) {
            if (abs(nowPivot) < abs(A_TA.values[j][i])) {
                index = j;
                nowPivot = A_TA.values[j][i];
            }
        }
        if (index != -1) {
            Identity eI(A_TA.n, A_TA.n, vector<vector<double>>());
            Permutation p(A_TA.n, A_TA.n, vector<vector<double>>(), &eI, i + 1, index + 1);
            A_TA = *(p * A_TA);
            A_TAREADY = p * *A_TAREADY;
        }
        for (int j = i + 1; j < A_TA.n; j++) {
            if (A_TA.values[j][i] != 0) {
                Identity eI(A_TA.n, A_TA.n, vector<vector<double>>());
                Elimination p(A_TA.n, A_TA.n, vector<vector<double>>(), &A_TA, &eI, j + 1, i + 1);
                A_TA = *(p * A_TA);
                A_TAREADY = p * *A_TAREADY;
            }
        }
    }
    for (int i = A_TA.n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (A_TA.values[j][i] != 0) {
                Identity eI(A_TA.n, A_TA.n, vector<vector<double>>());
                Elimination p(A_TA.n, A_TA.n, vector<vector<double>>(), &A_TA, &eI, j + 1, i + 1);
                A_TA = *(p * A_TA);
                A_TAREADY = p * *A_TAREADY;
            }
        }
    }

    for (int i = 0; i < A_TA.n; i++) {
        double delitel = A_TA.values[i][i];
        if (delitel != 0) {
            for (int j = 0; j < A_TA.m; j++) {
                A_TA.values[i][j] /= delitel;
                A_TAREADY->values[i][j] /= delitel;
            }
        }
    }

    cout << *A_TAREADY;

    cout << "A_T*b:" << endl;
    Matrix A_TB = *(A_T * B);
    cout << A_TB;

    cout << "x~:" << endl;
    Matrix ANS = *(*A_TAREADY * A_TB);
    cout << ANS;

    fprintf(pipe, "plot [-7 : 10] [-10 : 10] f(x)=%lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0, f(x) title 'f(x)', '-' using 1:2 with points title 'Point'\n", ANS.values[0][0], ANS.values[0][1], ANS.values[0][1], ANS.values[0][2]);
    for (int i = 0; i < n; i++) {
        fprintf(pipe, "%f\t%f\n", saveA[i], saveB[i]);
    }
    fprintf(pipe, "e\n");

    fflush(pipe);
#ifdef WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif

    return 0;
}
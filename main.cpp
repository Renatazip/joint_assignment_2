#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>

#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif


using namespace std;
int cnt = 0;

template<class T>
class IdentityMatrix;

template<class T>
class EliminationMatrix;

template<class T>
class PermutationMatrix;


template<class T>
class BaseMatrix {
public:
    T data[100][100];
    int n, m;

    BaseMatrix(int n, int m) : n(n), m(m) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                data[i][j] = 0;
            }
        }
    }

    friend istream &operator>>(istream &stream, BaseMatrix<T> &inst) {
        stream >> inst.n >> inst.m;
        for (int i = 0; i < inst.n; ++i) {
            for (int j = 0; j < inst.m; ++j) {
                stream >> inst.data[i][j];
            }
        }
        return stream;
    }

    friend ostream &operator<<(ostream &stream, const BaseMatrix<T> &inst) {
        for (int i = 0; i < inst.n; ++i) {
            for (int j = 0; j < inst.m; ++j) {
                printf("%.4f ", inst.data[i][j]);
                //stream << inst.data[i][j] << " ";
            }
            stream << "\n";
        }
        return stream;
    }

    BaseMatrix<T> operator+(BaseMatrix<T> &b) {
        if ((n != b.n) || (m != b.m)) {
            throw invalid_argument("not equal dimensions");
        } else {
            BaseMatrix<T> new_matrix(n, m);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    new_matrix.data[i][j] = data[i][j] + b.data[i][j];
                }
            }

            return new_matrix;
        }
    }

    BaseMatrix<T> operator-(BaseMatrix<T> &b) {
        if ((n != b.n) || (m != b.m)) {
            throw invalid_argument("not equal dimensions");
        } else {
            BaseMatrix<T> new_matrix(n, m);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    new_matrix.data[i][j] = data[i][j] - b.data[i][j];
                }
            }
            return new_matrix;
        }
    }

    void swapTwoRows(int i, int j) {
        for (int k = 0; k < m; ++k) {
            T tmp = data[i][k];
            data[i][k] = data[j][k];
            data[j][k] = tmp;
        }
    }

    BaseMatrix<T> operator*(BaseMatrix<T> &b) {
        if (m != b.n) {
            throw invalid_argument("not multiplicand");
        } else {
            BaseMatrix<T> new_matrix(n, b.m);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < b.m; ++j) {
                    for (int k = 0; k < m; ++k) {
                        new_matrix.data[i][j] += data[i][k] * b.data[k][j];
                    }
                }
            }
            return new_matrix;
        }
    }

    BaseMatrix<T> concatinate(BaseMatrix<T> &b) {
        if (n != b.n) {
            throw invalid_argument("can't");
        } else {
            BaseMatrix<T> new_matrix(n, b.m + m);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < b.m + m; ++j) {
                    if (j < m) {
                        new_matrix.data[i][j] = data[i][j];
                    } else {
                        new_matrix.data[i][j] = b.data[i][j];
                    }
                }
            }
            return new_matrix;
        }

    }

    BaseMatrix<T> Transpose() {
        BaseMatrix<T> new_matrix(m, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                new_matrix.data[j][i] = data[i][j];
            }

        }
        return new_matrix;
    }

    BaseMatrix<T> front_eliminations(bool with_steps) {
        BaseMatrix<T> copy = (*this);
        if (with_steps) {
            cout << "Direct way:\n";
        }

        for (int i = 0; i < copy.n; i++) {
            int mx = i;
            for (int j = i + 1; j < copy.m; j++) {
                if (abs(copy.data[j][i]) > abs(copy.data[mx][i])) {
                    mx = j;
                }
            }

            if (mx != i) {
                copy.swapTwoRows(i, mx);
                cnt += 1;
                if (with_steps) {
                    cout << "step #" << cnt << ": permutation\n";
                    cout << copy;
                }
            }


            for (int j = i + 1; j < copy.n; j++) {
                if (copy.data[j][i] == 0) {
                    continue;
                }
                double term = copy.data[j][i] / copy.data[i][i];
                for (int k = 0; k < copy.m; ++k) {
                    copy.data[j][k] = copy.data[j][k] - copy.data[i][k] * term;
                }
                cnt += 1;
                if (with_steps) {
                    cout << "step #" << cnt << ": elimination\n";
                    cout << copy;
                }
            }
        }
        return copy;
    }

    BaseMatrix<T> back_eliminations(bool with_steps) {
        BaseMatrix<T> copy = (*this);
        if (with_steps) {
            cout << "Way back:\n";
        }
        for (int i = copy.n - 1; i >= 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                if (copy.data[j][i] == 0) {
                    continue;
                }
                double term = copy.data[j][i] / copy.data[i][i];
                for (int k = 0; k < copy.m; ++k) {
                    copy.data[j][k] = copy.data[j][k] - copy.data[i][k] * term;
                }
                cnt += 1;
                if (with_steps) {
                    cout << "step #" << cnt << ": elimination\n";
                    cout << copy;
                }
            }
        }
        return copy;
    }

    BaseMatrix<T> diagonal_normalization(bool with_steps) {
        BaseMatrix<T> copy = (*this);

        for (int i = 0; i < copy.n; ++i) {
            for (int j = copy.m; j >= 0; --j) {
                copy.data[i][j] = copy.data[i][j] / copy.data[i][i];
            }

        }
        if (with_steps) {
            cout << "Diagonal normalization:\n";
            cout << copy;
        }
        return copy;
    }

    BaseMatrix<T> inverse(bool with_steps = false) {
        BaseMatrix<T> augmentMatrix(n, m * 2);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m * 2; ++j) {
                if (j < m) {
                    augmentMatrix.data[i][j] = data[i][j];
                } else {
                    if (i == j - m) {
                        augmentMatrix.data[i][j] = 1;
                    } else {
                        augmentMatrix.data[i][j] = 0;
                    }
                }
            }
        }
        if (with_steps) {
            cout << "step #0: Augmented Matrix\n";
            cout << augmentMatrix;
        }

        BaseMatrix<T> tmp = augmentMatrix.front_eliminations(with_steps);
        tmp = tmp.back_eliminations(with_steps);
        tmp = tmp.diagonal_normalization(with_steps);
        BaseMatrix<T> inverse_res = BaseMatrix<T>(tmp.n, tmp.m / 2);
        if (with_steps) {
            cout << "result:\n";


        }
        for (int i = 0; i < tmp.n; ++i) {
            for (int j = 0; j < tmp.m / 2; ++j) {
                inverse_res.data[i][j] = tmp.data[i][tmp.m / 2 + j];
                if (with_steps) {
                    printf("%.2f ", tmp.data[i][tmp.m / 2 + j]);
                }
            }
            if (with_steps) {
                cout << "\n";
            }
        }
        return inverse_res;
    }

};

template<class T>
class SquareMatrix : public BaseMatrix<T> {
public:
    SquareMatrix(int n) : BaseMatrix<T>(n, n) {};

    friend istream &operator>>(istream &stream, SquareMatrix<T> &inst) {
        stream >> inst.n;
        inst.m = inst.n;
        for (int i = 0; i < inst.n; ++i) {
            for (int j = 0; j < inst.m; ++j) {
                stream >> inst.data[i][j];
            }
        }
        return stream;
    }


    double det() {
        BaseMatrix<T> after_forward_eliminations = BaseMatrix<T>::front_eliminations();
        double d = 1.0;
        for (int i = 0; i < after_forward_eliminations.n; ++i) {
            d *= after_forward_eliminations.data[i][i];
        }
        printf("result:\n%.2f", d);
        return d;
    }
};

template<class T>
class IdentityMatrix : public SquareMatrix<T> {
public:
    IdentityMatrix(int n) : SquareMatrix<T>(n) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    this->data[i][j] = 1;
                } else {
                    this->data[i][j] = 0;
                }
            }

        }
    };
};


template<class T>
class EliminationMatrix : public IdentityMatrix<T> {
public:
    EliminationMatrix(BaseMatrix<T> base) : IdentityMatrix<T>(base.n) {
        for (int i = 1; i < base.n; i++) {
            for (int j = 0; j < base.n; j++) {
                this->data[i][j] = -1 * base.data[i][j] / base.data[j][j];
            }
        }
    }
};


template<class T>
class PermutationMatrix : public IdentityMatrix<T> {
public:
    PermutationMatrix(BaseMatrix<T> base) : IdentityMatrix<T>(base.n) {
        BaseMatrix<T> copy = base;

        for (int i = 0; i < copy.n; i++) {
            int mx = i;
            for (int j = i + 1; j < copy.m; j++) {
                if (abs(this->data[j][i]) > abs(this->data[mx][i])) {
                    mx = j;
                }
            }

            if (mx != i) {
                this->swapTwoRows(i, mx);
            }
        }
    }
};


int main() {
    FILE *pipe = popen(R"(C:\gnuplot\bin\gnuplot -persist)", "w");

    int m;
    cin >> m;
    vector<double> t(m), b(m);
    BaseMatrix<double> t_matrix(m, 1), b_matrix(m, 1);
    for (int i = 0; i < m; ++i) {
        cin >> t[i] >> b[i];
        b_matrix.data[i][0] = b[i];
        t_matrix.data[i][0] = t[i];
    }

    int n = 3;
    BaseMatrix<double> A(m, n + 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            A.data[i][j] = pow(t[i], j);
        }
    }

    BaseMatrix A_t = A.Transpose();
    BaseMatrix A_t_mul_a = A_t * A;
    BaseMatrix A_t_mul_a_inv = A_t_mul_a.inverse();
    BaseMatrix A_t_mul_b = A_t * b_matrix;
    BaseMatrix result = A_t_mul_a_inv * A_t_mul_b;

    fprintf(pipe, "plot [-2 : 8] [0 : 11] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' with points\n",
            result.data[3][0], result.data[2][0], result.data[1][0], result.data[0][0]);

    for (int i = 0; i < m; ++i) {
        fprintf(pipe, "%lf %lf\n", t_matrix.data[i][0], b_matrix.data[i][0]);
    }

    fprintf(pipe, "e\n");
    pclose(pipe);


    return 0;
}

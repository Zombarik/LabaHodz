#include <iostream>
#include <chrono>

#define N 4
using namespace std;

//=================================Bega Strassen algorithm==============================
void Print(double A[N][N], int n);
void ClassicMultiplication(double A[N][N], double B[N][N], double C[N][N], int n);
void Add(double A[N][N], double B[N][N], double C[N][N], int n);
void Subtract(double A[N][N], double B[N][N], double C[N][N], int n);
void Strassen(double A[N][N], double B[N][N], double C[N][N], int n);
//=================================Gordiichuk Gauss–Jordan elimination==============================
void inversion(double** A, int n);
//=================================Motrych LU inverse==============================
void subMatrix(double a[N][N], double temp[N][N], int p, int q, int n);
double Determinant(double a[N][N], int n);
double DeterminantLU(double l[N][N], double u[N][N], int n);
void LU(double a[N][N], double l[N][N], double u[N][N], int n);
void LX(double l[N][N], double x[N][N], int n);
void UA_(double u[N][N], double a_[N][N], double x[N][N], int n);
void LUinverse(double a[N][N], double l[N][N], double u[N][N], double a_[N][N], double x[N][N], int n);
//===============================Medvedchuk Newton==============================

//=======================================================================
void clearCin(); 

int main()
{
    double A[4][4] = {
        {3,4,5,1},
        {3,7,3,5},
        {2,-4,1,5},
        {7,2,3.5,8.2}
    };
    double B[4][4] = {
        {-1.2,3.3,4,1},
        {4,6,8,9},
        {3,6,7,-0.8},
        {1,4,3,6}
    };
    double L[4][4], U[4][4], A_[4][4], X[4][4];
    double C[4][4];
    while (true)
    {
        int input = 0;
        std::cout << "Choose(T**,R):\n[1] - Bega Strassen algorithm \n[2] - Gordiichuk GaussJordan elimination\n[3] - Motrych LU inverse\n[4] - Medvedchuk Newton\n" << std::endl;
        std::cin >> input;
        switch (input)
        {
        case 1:
        {
            system("cls");
            /*double A[4][4] = {
        {3,4,5,1},
        {3,7,3,5},
        {2,-4,1,5},
        {7,2,3.5,8.2}
            };
            double B[4][4] = {
                {-1.2,3.3,4,1},
                {4,6,8,9},
                {3,6,7,-0.8},
                {1,4,3,6}
            };
            double C[4][4];*/
            auto start = std::chrono::high_resolution_clock::now();
            Strassen(A, B, C, 4);
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            Print(C, 4);
            std::cout << "\nTook " << duration.count() << " microseconds";
            
            break;
        }
        case 2:
        {
            system("cls");
            int n;

            std::cout << "Enter N: ";
            std::cin >> n;

            double** matrix = new double* [n];

            for (int i = 0; i < n; i++)
                matrix[i] = new double[n];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    std::cout << "Enter matrix[" << i << "][" << j << "] = ";
                    std::cin >> matrix[i][j];
                }
            
                inversion(matrix, n);

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                        std::cout << matrix[i][j] << "  ";

                    std::cout << std::endl;
                }
            
            for (int i = 0; i < n; i++)
                delete[] matrix[i];

            delete[] matrix;

            //std::cin.get();
            break;
        }
        case 3:
        {
            system("cls");
            /*double A[4][4] = {
       {1,1,0,4},
       {2,1,3,2},
       {3,1,1,1},
       {3,4,1,6}
            };
            double L[4][4], U[4][4], A_[4][4], X[4][4];*/
            LUinverse(C, L, U, A_, X, 4);

            break;
        }
        case 4:
        {
            system("cls");
            
            break;
        }
        default:
        {
            clearCin();
        }
        }
    }
}
//=================================Bega Strassen algorithm==============================
void Print(double A[N][N], int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << "\n";
    }
}
void ClassicMultiplication(double A[N][N], double B[N][N], double C[N][N], int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            C[i][j] = 0;
            for (int k = 0; k < n; k++)
            {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}
void Add(double A[N][N], double B[N][N], double C[N][N], int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}
void Subtract(double A[N][N], double B[N][N], double C[N][N], int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}
void Strassen(double A[N][N], double B[N][N], double C[N][N], int n)
{
    double A11[N][N], A12[N][N], A21[N][N], A22[N][N],
        B11[N][N], B12[N][N], B21[N][N], B22[N][N],
        C11[N][N], C12[N][N], C21[N][N], C22[N][N],
        P1[N][N], P2[N][N], P3[N][N], P4[N][N], P5[N][N], P6[N][N], P7[N][N];

    if (n <= 2)
    {
        ClassicMultiplication(A, B, C, n);
        return;
    }
    for (int i = 0; i < n / 2; i++)
    {
        for (int j = 0; j < n / 2; j++)
        {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + n / 2];
            A21[i][j] = A[i + n / 2][j];
            A22[i][j] = A[i + n / 2][j + n / 2];
            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + n / 2];
            B21[i][j] = B[i + n / 2][j];
            B22[i][j] = B[i + n / 2][j + n / 2];
        }
    }

    double temp1[N][N], temp2[N][N];
    /*P1 = (A11 + A22)(B11 + B22)
    P2 = (A21 + A22)B11
    P3 = A11(B12 - B22)
    P4 = A22(B21 - B11)
    P5 = (A11 + A12)B22
    P6 = (A21 - A11)(B11 + B12)
    P7 = (A12 - A22)(B21 + B22)*/

    // P1 = (A11 + A22)(B11 + B22)
    Add(A11, A22, temp1, n / 2);
    Add(B11, B22, temp2, n / 2);
    Strassen(temp1, temp2, P1, n / 2);

    // P2 = (A21 + A22)B11
    Add(A21, A22, temp1, n / 2);
    Strassen(temp1, B11, P2, n / 2);

    // P3 = A11(B12 - B22)
    Subtract(B12, B22, temp2, n / 2);
    Strassen(A11, temp2, P3, n / 2);

    // P4 = A22(B21 - B11)
    Subtract(B21, B11, temp2, n / 2);
    Strassen(A22, temp2, P4, n / 2);

    // P5 = (A11 + A12)B22
    Add(A11, A12, temp1, n / 2);
    Strassen(temp1, B22, P5, n / 2);

    // P6 = (A21 - A11)(B11 + B12)
    Subtract(A21, A11, temp1, n / 2);
    Add(B11, B12, temp2, n / 2);
    Strassen(temp1, temp2, P6, n / 2);

    // P7 = (A12 - A22)(B21 + B22)
    Subtract(A12, A22, temp1, n / 2);
    Add(B21, B22, temp2, n / 2);
    Strassen(temp1, temp2, P7, n / 2);

    /*
    C11  = P1 + P4 - P5 + P7
    C12 = P3 + P5
    C21 = P2 + P4
    C22 = P1 - P2 + P3 +P6*/

    // C11  = P1 + P4 - P5 + P7
    Add(P1, P4, temp1, n / 2);
    Subtract(P7, P5, temp2, n / 2);
    Add(temp1, temp2, C11, n / 2);

    // C12 = P3 + P5
    Add(P3, P5, C12, n / 2);

    // C21 = P2 + P4
    Add(P2, P4, C21, n / 2);

    // C22 = P1 - P2 + P3 +P6
    Subtract(P1, P2, temp1, n / 2);
    Add(P3, P6, temp2, n / 2);
    Add(temp1, temp2, C22, n / 2);

    for (int i = 0; i < n / 2; i++)
    {
        for (int j = 0; j < n / 2; j++)
        {
            C[i][j] = C11[i][j];
            C[i][j + n / 2] = C12[i][j];
            C[i + n / 2][j] = C21[i][j];
            C[i + n / 2][j + n / 2] = C22[i][j];
        }
    }
}
//===============================Gordiichuk Gauss–Jordan elimination==============================
void inversion(double** A, int n)
{
    double temp;

    double** E = new double* [n];

    for (int i = 0; i < n; i++)
        E[i] = new double[n];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    for (int k = 0; k < n; k++)
    {
        temp = A[k][k];

        for (int j = 0; j < n; j++)
        {
            A[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < n; i++)
        {
            temp = A[i][k];

            for (int j = 0; j < n; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = n - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = A[i][k];

            for (int j = 0; j < n; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = E[i][j];

    for (int i = 0; i < n; i++)
        delete[] E[i];

    delete[] E;
}
//===============================Motrych LU inverse==============================
void subMatrix(double a[N][N], double temp[N][N], int p, int q, int n) {
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j] = a[row][col];
                j++;
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
double Determinant(double a[N][N], int n) {
    double determinant = 0;
    if (n == 1) {
        return a[0][0];
    }
    if (n == 2) {
        return (a[0][0] * a[1][1]) - (a[0][1] * a[1][0]);
    }
    double temp[N][N];
    int sign = 1;
    for (int i = 0; i < n; i++) {
        subMatrix(a, temp, 0, i, n);
        determinant += sign * a[0][i] * Determinant(temp, n - 1);
        sign = -sign;
    }
    return determinant;
}
double DeterminantLU(double l[N][N], double u[N][N], int n) {
    double d = 1;
   
    for (int i = 0; i < n; i++)
        d = d * l[i][i]*u[i][i];

    return d;
}
/*
L*U=A
l11  0   0   0    u11 u12 u13 u14   a11 a12 a13 a14
l21 l22  0   0  *  0  u22 u23 u24 = a21 a22 a23 a24
l31 l32 l33  0     0   0  u33 u34   a31 a32 a33 a34
l41 l42 l43 l44    0   0   0  u44   a41 a42 a43 a44

l11=a11                                                                              u11=1  u12=a12/l11  u13=a13/l11              u14=a14/l11
l21=a21  l22=a22-l21*u12                                                                    u22=1        u23=a23/l22-l21*u13/l22  u24=a24/l22-l21*u14/l22
l31=a31  l32=a32-l31*u12  l33=a33-l31*u13-l32*u23                                                        u33=1                    u34=a34/l33-l31*u14/l33-l32*u24
l41=a41  l42=a42-l41*u12  l43=a43-l41*u13-l42*u23  l44=a44-l41*u14-l42*u24-l43*u34                                                u44=1

*/
void LU(double a[N][N], double l[N][N], double u[N][N], int n) {
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (j < i)
                l[j][i] = 0;
            else {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++) {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < n; j++) {
            if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++) {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
}
/*
L*U=A
L*U*A^-1=I   l11  0   0     x11  0   0
X=U*A^-1   L=l21 l22  0   X=x21 x22  0
L*X=I        l31 l32 l33    x31 x32 x33

x11=1/l11                           x12=0                       x13=0               x14=0
x21=(l21*x11)/-l22                  x22=1/l22                   x23=0               x24=0
x31=(l31*x11+l32*x21)/-l33          x32=(l32*x22)/-l33          x33=1/l33           x34=0
x41=(l41*x11+l42*x21+l43*x31)/-l44  x42=(l42*x22+l43*x32)/-l44  x43=(l43*x33)/-l44  x44=1/l44
....

x[ðÿäîê][ñòîâï÷èê]
*/
void LX(double l[N][N], double x[N][N], int n) {
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (j > i) {
                x[j][i] = 0;
                for (k = 0; k < j; k++)
                    x[j][i] = x[j][i] + x[k][i] * l[j][k];
                x[j][i] = x[j][i] / (-l[j][j]);
            }
            if (j < i)
                x[j][i] = 0;
            else
                x[i][i] = 1 / l[i][i];
        }

    }


}
/*
u11 u12 u13 u14   y11 y12 y13 y14   x11 x12 x13 x14
 0  u22 u23 u24 * y21 y22 y23 y24 = x21 x22 x23 x24
 0   0  u33 u34   y31 y32 y33 y34   x31 x32 x33 x34
 0   0   0  u44   y41 y42 y43 y44   x41 x42 x43 x44
 y11=(x11-u12*y21-u13*y31-u14*y41)/u11
 y21=(x21-u23*y31-u24*y41)/u22
 y31=(x31-u34*y41)/u33
 y41=x41/u44
*/
void UA_(double u[N][N], double a_[N][N], double x[N][N], int n) {
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--) {
            a_[j][i] = x[j][i];
            for (k = j + 1; k < n; k++)
                a_[j][i] = a_[j][i] - u[j][k] * a_[k][i];
            a_[j][i] = a_[j][i] / u[j][j];
        }
    }
}
void LUinverse(double a[N][N], double l[N][N], double u[N][N], double a_[N][N], double x[N][N], int n) {
    LU(a, l, u, n);
    if(DeterminantLU( l, u, n) ==0)
        std::cout << "Inverse matrix does not exist" << std::endl;
    else {
        auto start = std::chrono::high_resolution_clock::now();
        LU(a, l, u, n);
        LX(l, x, n);
        UA_(u, a_, x, n);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        Print(a_, n);
        std::cout << "\nTook " << duration.count() << " microseconds";
    }
}
//===============================Medvedchuk==============================

//========================================================
void clearCin()
{
    std::cout << "Incorrect input" << std::endl;
    std::cin.clear();
    while (std::cin.get() != '\n')
    {
        std::cin.ignore();
    }
}
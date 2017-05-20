import java.util.Arrays;

/**
 * Created by denis__larin on 19.05.17.
 */
public class Main {
    public static void main(String[] args) {
        /*
        * уравнение:
        * y''-e^xy'-2xy = x^3
        * краевые условия:
        * 2y(0.1)-2.5y'(0.1) = 0
        * 3y(1.3)-3.4(1.3) = 5;
        *
        *
        *
        * p(x) = e^x
        * q(x) = 2x
        * f(x) = x^3
        * N = 5
        * A = 0
        * B = 5
        * alpha0 = 2
        * alpha1 = -2.5
        * beta0 = 3
        * beta1 = -3.4
        * a = 0.1
        * b = 1.3
        * */
        double alpha0 = 2, alpha1 = -2.5, Ac = 0, Bc = 5, beta0 = 3, beta1 = -3.4;
        double a0 = 0.1, b0 = 1.3;
        int n = 4;
        double[][] A = new double[n+1][n+1];
        double [] B = new double[n+1]; //вектор B[n+1]
        double [] X = new double[n+1]; //вектор X[n+1]

        //сетка
        double h = (b0-a0)/n;
        //ось X
        for (int i = 0; i <= n; i++) {
            X[i] = a0+i*h;
        }
        //собираем матрицы A и B
        for (int i = 0; i <= n-2; i++) {
            A[i][i] = h*h*q(X[i]) - h*p(X[i])+1;
            A[i][i+1] = h*p(X[i]) -2;
            A[i][i+2] = 1;
            B[i] = h*h*f(X[i]);
        }
        A[n-1][0] = alpha0*h - alpha1;
        A[n-1][1] = alpha1;
        A[n][n-1] = -beta1;
        A[n][n] = beta0*h+beta1;
        B[n-1] = h*Ac;
        B[n] = h*Bc;

        //вывод A
        for (int i = 0; i <=n; i++) {
            for (int j = 0; j <=n ; j++) {
                System.out.println("A["+i+"]["+j+"] = " + A[i][j]+ " ");
            }
            System.out.println();
        }
        for (int i = 0; i <= n; i++) {
            System.out.println("B["+i+"] = " + B[i]);
        }
        System.out.println("\n");

        //решаем уравнение A*X1 = B;
        double[] X1 = new double[n+1]; //вектор X1[n+1]
        gauss(A,B,X1,n+1);//решение методом гаусса

        //вывод ответа
        for (int i = 0; i <=n; i++) {
            System.out.println("X["+i+"] = " + X1[i]);
        }
    }

    private static void gauss(double[][] A, double[] B, double[] X, int n) {
        int m = n+1;
        //новая матрица размером n n+1
        double[][] C = new double[n][n+1];
        //копируем все из A и B в С
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                C[i][j] = A[i][j];
            }
            C[i][n] = B[i];
        }
        //вывод матрицы С
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n+1; j++) {
                System.out.println("C["+i+"]["+j+"] = " + C[i][j]+ " ");
            }
            System.out.println();
        }
        //прямой ход метода гаусса
        for (int k = 0; k < n-1; k++) {
            for (int i = k + 1; i < n; i++) {
                for (int j = m - 1; j >= k; j--) {
                    C[i][j] = C[i][j] - C[i][k] * C[k][j] / C[k][k];
                }
            }
        }
        //обратный ход
        X[n-1] = A[n-1][m-2]/A[n-1][m-2];
        for (int i = n-2; i >=0 ; i--) {
            double s = 0;
            for (int j = i+1; j <m-1 ; j++) {
                s = s+C[i][j]*X[j];
            }
            X[i] = (C[i][m-1] - s)/C[i][i];
        }
    }

    private static double f(double x) {
        return x*x*x;
    }

    private static double p(double x) {
        return Math.pow(Math.E,x);
    }

    private static double q(double x) {
        return 2*x;
    }
}

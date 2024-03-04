#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <numeric>

using namespace std;

double LEFT_BORDER, RIGHT_BORDER, h, sigma, eps = 0.000000001;
int N, n0, m = 2, k;
vector<double> X, Y, polinom, gradient, V, tmp;
bool first;

double P(double point, vector<double>& array) {
    double result = 0;
    for (int i = 0; i <= k; i++) {
        result = array[k - i] + result * point;
    }
    return result;
}

void get_gradient() {
    for (int j = 0; j < k + 1; j++) {
        gradient[j] = 0;
        for (int i = 0; i < N + 1; i++) {
            gradient[j] += 2 * pow(X[i], j) * (P(X[i], polinom) - Y[i]);
        }
    } 
}

double G(vector<double>& array) {
    double result = 0;
    for (int i = 0; i < N + 1; i++) {
        result += pow((Y[i] - P(X[i], array)), 2);
        
    }
    return result;
}

void get_V() {
    double tmpVal;
    if (first) {
        get_gradient();
        V = gradient;
        for (size_t j = 0; j < V.size(); j++) {
            V[j] *= -1;
        }
    } else {
        tmpVal = inner_product(gradient.begin(), gradient.end(), gradient.begin(), 0.0);
        get_gradient();
        for (size_t j = 0; j < V.size(); j++) {
            V[j] = -gradient[j] + (inner_product(gradient.begin(), gradient.end(), gradient.begin(), 0.0) / tmpVal) * V[j];
        }
    }
}

double lambda() {
    double m1, m2, m3;
    vector<double> polinom_minus_V(polinom.size());
    for (int i = 0; i < polinom.size(); i++) {
        polinom_minus_V[i] = polinom[i] - V[i];
    }
    vector<double> polinom_plus_V(polinom.size());
    for (int i = 0; i < polinom.size(); i++) {
        polinom_plus_V[i] = polinom[i] + V[i];
    }
    m1 = G(polinom_minus_V);
    m2 = G(polinom);
    m3 = G(polinom_plus_V);
    return (m1 - m3) / (2 * (m1 - 2 * m2 + m3));
}

double get_sigma() {
    double q, get_sigma_val = 1000;
    for (int i = 0; i < tmp.size(); i++) {
        tmp[i] -= polinom[i];
    }
    for (int i = 0; i < k + 1; i++) {
        if (polinom[i] != 0) {
            q = abs(tmp[i] / polinom[i]);
            if (q < get_sigma_val) {
                get_sigma_val = q;
            }
        }
    }
    return get_sigma_val;
}
void plot() {
    FILE *gp = popen("gnuplot","w"); 
	if(!gp) {
        printf("Error opening pipe to GNU plot\n");
        return;
    }
    fprintf(gp, "set terminal png\n");  
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set output 'problem.png'\n");
    fprintf(gp, "plot 'points.txt' index 0 with linespoints pt 7 lt rgb 'red' title 'exact','points.txt' index 1 with linespoints pt 7 lt rgb 'yellow' title 'polynom1','points.txt' index 1 with linespoints pt 7 lt rgb 'green' title 'polynom2'\n"); 
    fclose(gp);
    remove("points.txt"); 
}

int main() {
    cout << "_____________________________Please enter data________________________________\n\n";
    cout << "c: ";
    cin >> LEFT_BORDER;
    cout << "d: ";
    cin >> RIGHT_BORDER;
    cout << "N: ";
    cin >> N;
    cout << "no: ";
    cin >> n0;
    cout << endl;

    h = (RIGHT_BORDER - LEFT_BORDER) / N;
    X.resize(N + 1);
    for (int i = 0; i < N + 1; i++) {
        X[i] = LEFT_BORDER + i * h;
    }

    Y.resize(N + 1);
    cout << "Select function:\n";
    cout << "1. y = 2.14*x^3 - 0.6*x^2 + 5.342*x + 4\n";
    cout << "2. y = -2.14*(x^2) + 0.44*x - 4\n";
    cout << "3. y = sin(x) * exp(-x)\n";
    cout << "Your decision is: ";
    int choice;
    cin >> choice;
    cout << endl;

    switch (choice) {
        case 1:
            // Function 1
            for (int i = 0; i < N + 1; ++i) {
                Y[i] = 3.14 * pow(X[i], 3) - 1.7 * pow(X[i], 2) + 6.145 * X[i] + 1;
            }
            break;
        case 2:
            // Function 2
            for (int i = 0; i < N + 1; ++i) {
                Y[i] = -3.15 * pow(X[i], 2) + 1.21 * X[i] - 1;
            }
            break;
        case 3:
            // Function 3
            for (int i = 0; i < N + 1; ++i) {
                Y[i] = sin(X[i]) * exp(-X[i]);
            }
            break;
        default:
            cout << "Invalid choice!" << endl;
            return 1;
    }

    polinom.resize(n0 + m);
    gradient.resize(n0 + m);
    V.resize(n0 + m);
    tmp.resize(n0 + m);

    for (int i = 0; i < n0 + m; i++) {
        polinom[i] = gradient[i] = V[i] = tmp[i] = 0;
    }
    std::ofstream ofs("points.txt");
    for (k = n0; k < n0 + m; k++) {
        int i = 0;
        first = true;
        sigma = 1000;

        cout << "n: " << k << endl;

        while (sigma >= eps) {
            tmp = polinom;
            i++;
            get_V();
            double Lambda = lambda();
            
            for (int j = 0; j < n0 + m; j++) {
                polinom[j] += Lambda * V[j];
            }
            if (first) {
                first = false;
            } else {
                sigma = get_sigma();
            }

            cout << "iter number " << i << ": |grad(F)|=" << sqrt(inner_product(gradient.begin(), gradient.end(), gradient.begin(), 0.0)) << endl;
        }
        cout << endl;
        cout << "iter X(i)       Y(i)       P(i)       d" << endl;

        for (int i = 0; i < N + 1; i++) {
            cout << i << "    " << X[i] << "    " << Y[i] << "    " << P(X[i], polinom) << "    " << abs(Y[i] - P(X[i], polinom)) << endl;
        }

        cout << endl;
        cout << "a0..an" << endl;
        for (int i = 0; i <= k; i++) {
            cout << "a" << i << " = " << polinom[i] << endl;
        }
        cout << endl;
        for (int i = 0; i < N + 1; i++) {
            ofs << X[i] << ' ' << P(X[i], polinom) << std::endl;
        }
        ofs << std::endl << std::endl;
    }
    for (int i = 0; i < N + 1; i++) {
        ofs << X[i] << ' ' << Y[i] << std::endl;
    }
    ofs.close();
    plot();
    return 0;
}

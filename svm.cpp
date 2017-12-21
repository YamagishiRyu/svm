#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include "quadprog++.hh"

using namespace std;

int degree = 2;
double sigma = 10;

double inner_product(int way, std::pair<double, double> x0, std::pair<double, double> x1);
string input_filename();
int input_kernel_type();

int main (int argc, char *const argv[]) {

  // hear file name
  string filename = input_filename();

  // hear kernel type
  int kernel = input_kernel_type();

  // file open
  ifstream ifs(filename);
  if(ifs.fail()){
    cerr << "File do not exist: \'" << filename << "\'." << endl;
    exit(0);
  }

  // file scan
  string str;
  double a = 0, b = 0, c = 0;
  double maxX = 0, minX = 0;
  double maxY = 0, minY = 0;
  vector< pair<double, double> > x;
  vector<double> y;

  while(getline(ifs,str)){
    a = b = c = 0;
    sscanf(str.data(), "%lf %lf %lf", &a, &b, &c);
    // max min 
    if(maxX < a){maxX = a;}
    if(minX > a){minX = a;}
    if(maxY < b){maxY = b;}
    if(minY > b){minY = b;}

    x.push_back(make_pair(a, b));
    if(c != 1 && c != -1){
      cerr << "Data type is not appropriate." << endl;
      exit(0);
    }
    y.push_back(c);
  }

  cout << "size: " << x.size() << ", range: (" << minX << ", " << minY << ")~(" << maxX << ", " << maxY << ")." << endl;

  // quadratic programming
  double G[MATRIX_DIM][MATRIX_DIM], g0[MATRIX_DIM], 
          CE[MATRIX_DIM][MATRIX_DIM], ce0[MATRIX_DIM], 
          CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM], 
          alpha[MATRIX_DIM];
  int n, m, p;
  char ch;
  
  n = x.size();
  double delta = 0;
  if(kernel != 1){
    delta = 1.0e-9;
  }else{
    delta = 1.0e-7;
  }
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i == j){
        G[i][j] = pow(y[i], 2.0) * inner_product(kernel, x[i], x[i]) + delta;
      }else{
        double result = y[i] * y[j] * inner_product(kernel, x[i], x[j]);
        G[i][j] = result;
        G[j][i] = result;
      }
    }
    g0[i] = 0 - 1.0;
  }

  p = 1;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      CE[i][j] = y[i];
      ce0[j] = 0;
    }
  }

  m = n;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      if(i == j){
        CI[i][j] = 1.0;
      }else{
        CI[i][j] = 0.0;
      }
    }
    ci0[i] = 0;
  }

  // solve quadprog programming
  cout << "f: " << solve_quadprog(G, g0, n, CE, ce0, p, CI, ci0, m, alpha) << endl;

  // double checking
  double sum = 0;
  for(int i = 0; i < n; i++){
    sum += alpha[i] * y[i];
  }
  cout << "sum: " << sum << endl;

  // weight
  pair<double, double> weight = make_pair(0, 0);
  for(int i = 0; i < n; i++){
    weight.first += alpha[i] * y[i] * x[i].first;
    weight.second += alpha[i] * y[i] * x[i].second;
  }
  cout << "w: " << weight.first << ", " << weight.second << endl;

  // theta
  double maxA = alpha[0];
  int maxAindex = 0;
  for(int i = 0; i < n; i++){
    cout << "alpha[" << i << "] : " << alpha[i] << endl;
    if(alpha[i] > maxA){
      maxA = alpha[i];
      maxAindex = i;
    }
  }
  cout << "most index : " << maxAindex << endl;
  double theta = 0;
  for(int i = 0; i < n; i++){
    theta += alpha[i] * y[i] * inner_product(kernel, x[i], x[maxAindex]);
  }
  theta -= y[maxAindex];
  cout << "θ: " << theta << endl;

  // gnuplot plot
  FILE *gp;

  gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set multiplot\n");
  fprintf(gp, "set xrange [%f:%f]\n", minX, maxX);
  fprintf(gp, "set yrange [%f:%f]\n", minY, maxY);
  
  fprintf(gp, "set lmargin at screen 0.1\n");
  fprintf(gp, "set rmargin at screen 0.9\n");
  fprintf(gp, "set bmargin at screen 0.1\n");
  fprintf(gp, "set tmargin at screen 0.9\n");
  fprintf(gp, "set zrange [-10:10]\n");
  fprintf(gp, "unset xtics\n");
  fprintf(gp, "unset ytics\n");
  fprintf(gp, "set isosample 200, 200\n");
  fprintf(gp, "set ticslevel 1 / 4 \n");
  fprintf(gp, "set contour base\n");
  fprintf(gp, "set cntrparam levels discrete 0.0\n");
  fprintf(gp, "set nosurface\n");
  fprintf(gp, "set zeroaxis ls -1\n");
  fprintf(gp, "set view 0, 0\n");

  // margin plot
  if(kernel == 0){
    // linear
    fprintf(gp, "f(x, y) = 0");
    for(int i = 0; i < n; i++){
      fprintf(gp, " + %f * %f * (%f * x + %f * y)", alpha[i], y[i], x[i].first, x[i].second);
    }
    fprintf(gp, " - %f\n", theta);
    fprintf(gp, "splot f(x, y) with lines lt 2 lw 4 lc 5 title 'f'\n");
  }else if(kernel == 1){
    // polynomial
    fprintf(gp, "f(x, y) = 0");
    for(int i = 0; i < n; i++){
      fprintf(gp, " + %f * %f * (1 + %f * x + %f * y) ** %d", alpha[i], y[i], x[i].first, x[i].second, degree);
    }
    fprintf(gp, " - %f\n", theta);
    fprintf(gp, "splot f(x, y) with lines lt 2 lw 4 lc 5 title 'f'\n");
  }else{
    // Gauss
    /*
    fprintf(gp, "f(x, y) = 0");
    for(int i = n; i < n; i++){
      fprintf(gp, " + %f * %f * exp(0 - ((%f - x) ** 2 + (%f - y) ** 2) / 2 / (%f ** 2))", alpha[i], y[i], x[i].first, x[i].second, sigma);
    }
    fprintf(gp, " - %f\n", theta);
    */
    double dx = (maxX - minX) / 100;
    double dy = (maxY - minY) / 100;
    double epsilon = 0.4;
    double nowX, nowY;
    double sign;
    for(int i = 0; i < 100; i++){
      for(int j = 0; j < 100; j++){
        nowX = minX + dx * i;
        nowY = minY + dy * j;
        sign = 0;
        for(int k = 0; k < n; k++){
          sign += alpha[k] * y[k] * inner_product(kernel, x[k], make_pair(nowX, nowY));
        }
        sign -= theta;
        if(sign > - epsilon && sign < epsilon){
          fprintf(gp, "set object %d rect from %f, %f to %f, %f back linewidth 0 fillcolor rgb \"black\" fill solid 0.5\n", i * 100 + j, nowX, nowY, nowX + dx, nowY + dy);
        }
      }
    }
  }

  fprintf(gp, "set xtics\n");
  fprintf(gp, "set ytics\n");
  for(int i = 0; i < x.size(); i++){
    if(y[i] == 1){
      fprintf(gp, "plot '-' with points pt 9 lc 10 title '1' at 0.90, 0.93\n");
    }else{
      fprintf(gp, "plot '-' with points pt 1 lc 20 title '-1' at 0.90, 0.90\n");
    }

    fprintf(gp, "%f\t%f\n", x[i].first, x[i].second);
    fprintf(gp, "e\n");
  }

  fprintf(gp, "unset multiplot\n");
  pclose(gp);
  return 0;
}

string input_filename(){
  string filename;
  cout << "What is data file name ?" << endl;
  cin >> filename;
  return filename;
}

int input_kernel_type(){
  int kernel;
  cout << "Choose kernel type(0: none, 1: polynomial, 2: Gauss)." << endl;
  cin >> kernel;
  if(kernel < 0 || kernel > 2){
    cerr << "kernel type is not appropriate." << endl;
    exit(0);
  }
  return kernel;
}

double inner_product(int way, std::pair<double, double> x0, std::pair<double, double> x1){
  double ip = 0;
  if(way == 0){
    // linear
    ip = x0.first * x1.first + x0.second * x1.second;
  }else if(way == 1){
    // polynomial
    ip = pow(1 + x0.first * x1.first + x0.second * x1.second, degree);
  }else{
    // Gauss
    ip = exp(0 - (pow(x0.first - x1.first, 2) + pow(x0.second - x1.second, 2)) / 2 / pow(sigma, 2));
  }

  return ip;
}

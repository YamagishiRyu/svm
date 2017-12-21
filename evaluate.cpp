#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include "quadprog++.hh"

using namespace std;

// function definition
double inner_product(int way, double degree, double sigma, std::pair<double, double> x0, std::pair<double, double> x1);
string input_filename();
int input_kernel_type();
int input_divide_num();
int import_data(std::ifstream *ifs, std::vector< std::pair<double, double> > *x, std::vector<double> *y);
void solve(int kernel, double degree, double sigma, vector< pair<double, double> > *x, vector<double> *y, double alpha[MATRIX_DIM]);
double calculate_theta(int kernel, double degree, double sigma, vector< pair<double, double> > *x, vector<double> *y, double alpha[MATRIX_DIM]);
double test(int kernel, double degree, double sigma, vector< pair<double, double> > *train_x, vector<double> *train_y, vector< pair<double, double> > *test_x, vector<double> *test_y);
double calculate_score(int kernel, double degree, double sigma, int  divide_num, vector< pair<double, double> > *x, vector<double> *y);

int main (int argc, char *const argv[]) {
  
  string filename = input_filename();
  // file open
  ifstream ifs(filename);
  if(ifs.fail()){
    cerr << "File do not exist: \'" << filename << "\'." << endl;
    exit(0);
  }

  int kernel = input_kernel_type();
  int divide_num = input_divide_num();

  // import data
  vector< pair<double, double> > x;
  vector<double> y;
  int data_size = import_data(&ifs, &x, &y);
  cout << "data size: " << x.size() << endl;

  // parameter
  double degree = 0;
  double sigma = 0;
  double score;

  if(kernel == 0){
    score = calculate_score(kernel, degree, sigma, divide_num, &x, &y);
    cout << "score: " << score << endl;
  }else if(kernel == 1){
    for(int i = 1; i < 11; i++){
      score = calculate_score(kernel, i, sigma, divide_num, &x, &y);
      cout << "degree: " << i << ", score: " << score << endl;
    }
  }else{
    for(int i = 1; i < 10; i++){
      score = calculate_score(kernel, degree, i, divide_num, &x, &y);
      cout << "sigma: " << i << ", score: " << score << endl;
    }
  }

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

int input_divide_num(){
  int divide_num = 0;
  cout << "How many groups sample data will be divided into ?" << endl;
  cin >> divide_num;
  if(divide_num <= 0){
    cerr << "Division number should be natural number." << endl;
    exit(0);
  }
  return divide_num;
}

int import_data(std::ifstream *ifs, std::vector< std::pair<double, double> > *x, std::vector<double> *y){
  
  string str;
  double a = 0, b = 0, c = 0;
  
  while(getline(*ifs,str)){
    a = b = c = 0;
    sscanf(str.data(), "%lf %lf %lf", &a, &b, &c);

    x->push_back(make_pair(a, b));
    if(c != 1 && c != -1){
      cerr << "Data type is not appropriate." << endl;
      exit(0);
    }
    y->push_back(c);
  }

  return x->size();
}

void solve(int kernel, double degree, double sigma, vector< pair<double, double> > *x, vector<double> *y, double alpha[MATRIX_DIM]){
  // quadratic programming
  double G[MATRIX_DIM][MATRIX_DIM], g0[MATRIX_DIM], 
          CE[MATRIX_DIM][MATRIX_DIM], ce0[MATRIX_DIM], 
          CI[MATRIX_DIM][MATRIX_DIM], ci0[MATRIX_DIM];
  int n, m, p;
  
  n = x->size();
  double delta = 0;
  if(kernel != 1){
    delta = 1.0e-9;
  }else{
    delta = 1.0e-7;
  }
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i == j){
        G[i][j] = pow(y->at(i), 2.0) * inner_product(kernel, degree, sigma, x->at(i), x->at(i)) + delta;
      }else{
        double result = y->at(i) * y->at(j) * inner_product(kernel, degree, sigma, x->at(i), x->at(j));
        G[i][j] = result;
        G[j][i] = result;
      }
    }
    g0[i] = 0 - 1.0;
  }

  p = 1;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      CE[i][j] = y->at(i);
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
  double f = solve_quadprog(G, g0, n, CE, ce0, p, CI, ci0, m, alpha);
}

double inner_product(int way, double degree, double sigma, std::pair<double, double> x0, std::pair<double, double> x1){
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

double calculate_theta(int kernel, double degree, double sigma, vector< pair<double, double> > *x, vector<double> *y, double alpha[MATRIX_DIM]){
  double maxA = alpha[0];
  int maxAindex = 0;
  int n = x->size();
  for(int i = 0; i < n; i++){
    if(alpha[i] > maxA){
      maxA = alpha[i];
      maxAindex = i;
    }
  }
  double theta = 0;
  for(int i = 0; i < n; i++){
    theta += alpha[i] * y->at(i) * inner_product(kernel, degree, sigma, x->at(i), x->at(maxAindex));
  }
  theta -= y->at(maxAindex);
  return theta;
}

double calculate_score(int kernel, double degree, double sigma, int  divide_num, vector< pair<double, double> > *x, vector<double> *y){
  double score = 0; 
  vector< pair<double, double> > train_x;
  vector<double> train_y;
  vector< pair<double, double> > test_x;
  vector<double> test_y;
  
  for(int i = 0; i < divide_num; i++){
    // clear
    train_x.clear();
    train_y.clear();
    test_x.clear();
    test_y.clear();

    // divide
    for(int j = 0; j < x->size(); j++){
      if(j % divide_num == i){
        test_x.push_back(x->at(j));
        test_y.push_back(y->at(j));
      }else{
        train_x.push_back(x->at(j));
        train_y.push_back(y->at(j));
      }
    }

    // test
    score += test(kernel, degree, sigma, &train_x, &train_y, &test_x, &test_y);
  }
  return score / divide_num;
}

double test(int kernel, double degree, double sigma, vector< pair<double, double> > *train_x, vector<double> *train_y, vector< pair<double, double> > *test_x, vector<double> *test_y){
  
  // solve quadprog
  double alpha[MATRIX_DIM];
  solve(kernel, degree, sigma, train_x, train_y, alpha);

  // theta 
  double theta = calculate_theta(kernel, degree, sigma, train_x, train_y, alpha);
  
  int success_num = 0;
  double sign = 0;
  for(int i = 0; i < test_y->size(); i++){
    // calculate sign
    sign = 0;
    for(int j = 0; j < train_x->size(); j++){
      sign += alpha[j] * train_y->at(j) * inner_product(kernel, degree, sigma, train_x->at(j), make_pair(test_x->at(i).first, test_x->at(i).second));
    }
    sign -= theta;

    // count the number of success expample
    if(test_y->at(i) == 1 && sign > 0){
      success_num++;
    }else if(test_y->at(i) == -1 && sign < 0){
      success_num++;
    }
  }
  double result = success_num / (double)test_y->size();
  return result;
}



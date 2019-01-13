#include <cmath>
#include <vector>
#include "mesh.hpp"
#include "boundary_value_condition.hpp"
using std::vector;

#ifndef CRANK_NICOLSON_SOLVER
#define CRANK_NICOLSON_SOLVER

const double INF = 1e18;
const double EPS = 1e-12;

void gauss(vector<vector<double>>& augmented_matrix);



class crank_nicolson_solver {
    public:
    double theta;
    double dt;
    double dx;

    int T;

    mesh initial_mesh;

    vector<mesh> answer;

    boundary_value_condition* bv;

    double (*sigma)(double,double);
    double (*r)(double,double);

    crank_nicolson_solver(double theta,double (*sigma)(double,double),
    double (*r)(double,double),boundary_value_condition* bv,int T,mesh initial_mesh);    

    void solve();

};

#endif
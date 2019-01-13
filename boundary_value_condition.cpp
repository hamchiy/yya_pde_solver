#include "boundary_value_condition.hpp"


boundary_value_condition::boundary_value_condition(int N,int T) : N(N),T(T) {}

void boundary_value_condition::set_values(const vector<double>& v0, const vector<double>& v1) {
    this->v1 = v1;
    this->v0 = v0;
}

dirichlet_boundary_condition::dirichlet_boundary_condition(int N,int T) : boundary_value_condition(N,T) {};

void dirichlet_boundary_condition::get_eqn_at_time(double t,vector<double> & X_0,vector<double>& X_N) {
    X_0[0] = 1;
    X_0[N] = v0[t];

    X_N[N-1] = 1;
    X_N[N] = v1[t];
}



neumann_boundary_condition::neumann_boundary_condition(int N,int T,double dx) : boundary_value_condition(N,T), dx(dx) {};

void neumann_boundary_condition::get_eqn_at_time(double t,vector<double> & X_0,vector<double>& X_N) {
    X_0[1] = 1;
    X_0[0] = -1;
    X_0[N] = v0[t]*dx;

    X_N[N-1] = 1;
    X_N[N-2] = -1;
    X_N[N] = v1[t]*dx;

}

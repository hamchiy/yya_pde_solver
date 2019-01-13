#include "crank_nicolson_solver.hpp"
using namespace std;


void gauss(vector<vector<double>>& augmented_matrix) {
    int n = augmented_matrix.size();
    for(int i=0; i<n; ++i) {
        // find row from i with largest value at i
        int tr = i;
        for(int j=i+1; j<n; ++j) {
            if(fabs(augmented_matrix[j][i]) > fabs(augmented_matrix[tr][i])) tr = j;
        }
        // swap with row i
        for(int j=0; j<=n; ++j) swap(augmented_matrix[i][j], augmented_matrix[tr][j]);
       
        // eliminate i from all rows below
        for(int j=i+1; j<n; ++j) {
            for(int k=n; k>=i; --k) {
                augmented_matrix[j][k] -= augmented_matrix[i][k] * augmented_matrix[j][i] / augmented_matrix[i][i]; // be careful here
            }
        }
    }
    // restore
    for(int i=n-1; i>=0; --i) {
        double t = augmented_matrix[i][n];
        for(int j=i+1; j<n; ++j) t -= augmented_matrix[i][j] * augmented_matrix[j][n];
        t /= augmented_matrix[i][i];
        augmented_matrix[i][n] = t;
    }
    // value of the i'th variable is in augmented_matrix[i][n]
}




crank_nicolson_solver::crank_nicolson_solver(double theta,double (*sigma)(double,double),
double (*r)(double,double),boundary_value_condition* bv,int T,mesh initial_mesh) {
    this->theta = theta;
    this->sigma = sigma;
    this->r = r;
    this->bv = bv;
    this->T = T;
    this->initial_mesh = initial_mesh;
    dx = initial_mesh.range/(initial_mesh.arr[0].size()-1);
    dt = 1;

    answer.resize(T+1);
}    

void crank_nicolson_solver::solve() {
    answer[T].init(initial_mesh.N,initial_mesh.x_0,initial_mesh.range);
    answer[T].init(initial_mesh.arr[1]);
    int n = initial_mesh.N;

    vector<double> L(n);
    
    //storing L(i,T) at L[i] for i=1 to (n-2)
    for(int i=1;i<n-1;i++) {
        double f_curr = initial_mesh.get_val_at_index(i);
        double f_next = initial_mesh.get_val_at_index(i+1);
        double f_prev = initial_mesh.get_val_at_index(i-1);
        double x = initial_mesh.get_xval(i);
        double s = sigma(x,T);
        double r_val = r(x,T);
        L[i] = -0.5*s*s*(f_next-2*f_curr+f_prev)/(dx*dx) + 
        0.5*(s*s-r_val)*(f_next-f_prev)/(2*dx) + r_val*f_curr;
    }
    

    vector<vector<double>> augmented_matrix(n);
    for(int i=0;i<n;i++) {
        augmented_matrix[i].resize(n+1);
    }

    for(int t=T-1;t>=0;t--) {
        answer[t].init(initial_mesh.N,initial_mesh.x_0,initial_mesh.range);

        //clearing value of augmented matrix
        for(int i=0;i<n;i++) {
            augmented_matrix[i] = vector<double>(n+1,0);
        }
        
        bv->get_eqn_at_time(t,augmented_matrix[0],augmented_matrix[n-1]);
        for(int i=1;i<n-1;i++) {
            //filling the matrix
            double x = answer[t].get_xval(i);
            double s = sigma(x,t);
            double r_val = r(x,t);


            augmented_matrix[i][i+1] = -0.5*s*s*dt/(dx*dx)*theta + 
            0.5*(s*s-r_val)*dt/(2*dx)*theta;

            augmented_matrix[i][i] = -0.5*s*s*(-2)*dt/(dx*dx)*theta + r_val*theta*dt + 1;

            augmented_matrix[i][i-1] = -0.5*s*s*dt/(dx*dx)*theta - 
            0.5*(s*s-r_val)*dt/(2*dx)*theta;
            
            augmented_matrix[i][n] = answer[t+1].get_val_at_index(i) - (1-theta)*L[i]*dt;
        }

        gauss(augmented_matrix);

        for(int  i=0;i<n;i++) {
            answer[t].arr[1][i] = augmented_matrix[i][n];
        }

        //setting value of L(i,t) in L[i] for all i=1..(N-1)
        for(int i=1;i<n-1;i++) {
            double f_curr = answer[t].get_val_at_index(i);
            double f_next = answer[t].get_val_at_index(i+1);
            double f_prev = answer[t].get_val_at_index(i-1);
            double x = answer[t].get_xval(i);
            double s = sigma(x,T);
            double r_val = r(x,T);
            L[i] = -0.5*s*s*(f_next-2*f_curr+f_prev)/(dx*dx) + 
            0.5*(s*s-r_val)*(f_next-f_prev)/(2*dx) + r_val*f_curr;
        }
    }
}




#include "mesh.hpp"

mesh::mesh() {}
mesh::mesh(const mesh& m) {
    N = m.N;
    x_0 = m.x_0;
    range = m.range;
    arr[0] = m.arr[0];
    arr[1] = m.arr[1];
}

void mesh::init(int N,double x_0,double range) {
    this->N = N;
    this->x_0 = x_0;
    this->range = range;
    
    arr[0].resize(N); //x value
    arr[1].resize(N); //the value of f(x) at a certain time
    for(int i=0;i<N;i++) {
        arr[0][i] = x_0 - range/2 + range/(N-1)*i;
    }
}

void mesh::init(const vector<double>& v) {
    for(int i=0;i<N;i++) {
        arr[1][i] = v[i];
    }
}

double mesh::get_val_at_index(int i) {
    return arr[1][i];
}

double mesh::get_xval(int i) {
    return arr[0][i];
}

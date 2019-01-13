#include <vector>
using std::vector;


#ifndef MESH_HPP
#define MESH_HPP

class mesh {
public:
    int N;
    double x_0;
    double range;
    
    vector<double> arr[2];

    mesh();
    mesh(const mesh& m);
    void init(int N,double x_0,double range);
    void init(const vector<double>& v);
    double get_val_at_index(int i);
    double get_xval(int i);
};

#endif
#include <iostream>
#include "crank_nicolson_solver.hpp"
#include "closed_form.hpp"

using namespace std;

double sigma(double x,double t) {
    return 1;
}

double r(double x,double t) {
    return 1;
}

int main() {

	double spot;
	double strike;
	double volatility;
	double maturity;
	bool istrue = true;

	//Inital inputs for closed form 
	std::cout << "Underlying Spot Price:" << std::endl;
	std::cin >> spot;
	std::cout << "Strike Price:" << std::endl;
	std::cin >> strike;
	std::cout << "Volatility:" << std::endl;
	std::cin >> volatility;
	std::cout << "Maturity (in days):" << std::endl;
	std::cin >> maturity;
	std::cout << " " << std::endl;

	std::cout << "BS Price is:" << dauphine::bs_price(spot, strike, volatility, maturity, istrue) << std::endl;

    int N,T;
	cout << "Inputs for Crank Nicolson solver" << endl;
    cout << "Enter the number of Points in mesh : " << endl;
    cin >> N;

    cout << "Enter number of Days : " << endl;
    cin >> T;

    double x_0;
    cout << "Enter Central Value x_0 : " << endl;
    cin >> x_0;

    double range;
    cout << "Enter Range : " << endl;
    cin >> range;

    mesh init_mesh;

    vector<double> v(N);
    cout << "Enter the values of initial mesh : " << endl;
    for(int i=0;i<N;i++) cin >> v[i];

    init_mesh.init(N,x_0,range);
    init_mesh.init(v);

    cout << "Select Boundary Option : \n1)Dirichlet Boundary Condition\n" << 
    "2)Neumann Boundary Condition" << endl;

    int opt;
    cin >> opt;

    boundary_value_condition* bv;
    if(opt==1) {
        bv = new dirichlet_boundary_condition(N,T);
        
    } else {
        bv = new neumann_boundary_condition(N,T,range);
    }

    vector<double> v0(T+1),v1(T+1);
    cout << "Enter boundary values for i=0 : " << endl;
    for(int i=0;i<=T;i++) cin >> v0[i];

    cout << "Enter boundary values for i=N-1 : " << endl;
    for(int i=0;i<=T;i++) cin >> v1[i];

    bv->set_values(v0,v1);

    crank_nicolson_solver s(0.5,sigma,r,bv,T,init_mesh);
    s.solve();

    for(int i=0;i<=T;i++) {
        for(int j=0;j<N;j++) {
            cout << s.answer[i].get_val_at_index(j) << " ";
        }
        cout << endl;
    }





    




    return 0;
}
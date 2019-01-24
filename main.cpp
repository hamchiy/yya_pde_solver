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

	//initial inputs (inputs Spot = 100 / Strike = 10 / Vol = 0.2 / Maturity = 1 )
	std::cout << "Underlying Spot Price:" << std::endl;
	std::cin >> spot;
	std::cout << "Strike Price:" << std::endl;
	std::cin >> strike;
	std::cout << "Volatility:" << std::endl;
	std::cin >> volatility;
	std::cout << "Maturity (in years):" << std::endl;
	std::cin >> maturity;
	std::cout << " " << std::endl;

	std::cout << "BS Price is:" << BS::bs_price(spot, strike, volatility, maturity, istrue) << std::endl;
    int N,T;
    cout << "Enter the number of Points in mesh : " << endl;
    cin >> N;

    cout << "Enter number of years : " << endl;
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
    
    cin >> v0[0];
    for(int i=1;i<=T;i++) v0[i] = v0[0];
    
    cout << "Enter boundary values for i=N-1 : " << endl;
    cin >> v1[0];
    for(int i=1;i<=T;i++) v1[i] = v1[0];

    bv->set_values(v0,v1);

    double theta;
    cout << "Enter the value of Theta" << endl;
    cin >> theta;

    crank_nicolson_solver s(theta,sigma,r,bv,T,init_mesh);
    s.solve();

    for(int i=0;i<=T;i++) {
        for(int j=0;j<N;j++) {
            cout << s.answer[i].get_val_at_index(j) << " ";
        }
        cout << endl;
    }
    return 0;
}
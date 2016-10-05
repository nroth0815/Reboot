#include <numeric>
#include <iostream>
#include <vector>

using namespace std;

double dp(vector<double> v1, vector<double> v2){
  //int v1[] = { 1, 2, 3 };
  //int v2[] = { 4, 6, 8 };
  //cout << "the dot product of (1,2,3) and (4,6,8) is ";
  //cout << inner_product(begin(v1), end(v1), begin(v2), 0.0) << endl;
  return inner_product(begin(v1), end(v1), begin(v2), 0.0);//inner_product(v1, v1+3, v2, 0.0);
}

double dp_exp(vector<double> v1, vector<double> v2){

	double prodct=0.;

	for (int i=0; i<3; i++){
		prodct+=(v1[i]*v2[i]);
	}

	return prodct;

}


int main() {

	double vv1[3] = { 1, 2, 3 };
    double vv2[3] = { 4, 6, 8 };

	vector<double> v1 (begin(vv1), end(vv1));
    vector<double> v2 (begin(vv2), end(vv2));
    double res1, res2;

	//for (int j=0; j<10000; j++){
		res1=dp(v1, v2);
		res2=dp_exp(v1, v2);
	//	v1[0]+=1.;
	//}


	cout<< res1 << " " << res2<< endl;

return 0;
}

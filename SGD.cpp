#include <iostream>
#include <cmath>
#include <random>
#include <armadillo>
#include <map>
#include <vector>
#include <set>
#include <fstream>
#include <string>
#include <map>
#include <float.h>
#include <time.h>
using namespace std;
using namespace arma;

// The regularization term
#define LAMBDA  (1e-4)

std::mt19937 generator;

int iteration_num = 0;



void Linear_grad(vec& beta, vec& X, double y, vec& grad){
	double z = dot(beta, X) - y;
	grad = X * (z);          
}
void Logistic_grad(vec& beta, vec& X, double y, vec& grad){
	double z = y * dot(beta, X);
    if (z < 0) 
		grad = -1.0 / ( 1.0 + exp(z) )  * y * X;
    else  grad = - exp(-z) / ( 1.0 + exp(-z) ) * y * X;
}
void SVM_grad(vec& beta, vec& X, double y, vec& grad){
	double z = y * dot(beta, X);
    grad.zeros();
    if (z < 1.0)
       grad = - y * X;
}

void splitString(string s, char delimeter, vector<string>& a){
	a.clear();
	int loc = 0;
	int loc1;
	string str;
	int i =0;
	while(s.find(delimeter, loc) != string::npos){
		loc1 = s.find(delimeter, loc);
		str = s.substr(loc, loc1 - loc);
		a.push_back(str);
		loc = loc1 + 1;
	}
	if(loc != s.size()){
		str= s.substr(loc);
		a.push_back(str);
	}
}

void splitStringintoInt(string s, char delimeter, vector<int>& a){
	a.clear();
	int loc = 0;
	int loc1;
	string str;
	int i =0;
	while(s.find(delimeter, loc) != string::npos){
		loc1 = s.find(delimeter, loc);
		str = s.substr(loc, loc1 - loc);
		a.push_back(atoi(str.c_str()));
		loc = loc1 + 1;
	}
	if(loc != s.size()){
		str= s.substr(loc);
		a.push_back(atoi(str.c_str()));
	}
}


void loadData(string fdomain, string floc, mat& NX, rowvec& YY_linear, rowvec& YY_logistic){
	map<int, int> Category; 
	map<int, pair<double, double> > Numeric;

	fstream f; f.open(fdomain);
	string line;
	int count = 0;
	while(getline(f, line, '\n')){
		vector<string> a;
		splitString(line, ' ', a);
		if(a.at(0) == "C"){
			Category.insert(make_pair(count, atoi(a[1].c_str()) -1));
		}
		else{
			double min = atof(a[1].c_str());
			double max = atof(a[2].c_str());
			Numeric.insert(make_pair(count, pair<double, double>(min, max)));
		}
		count++;
	}
	f.close();
	
	f.open(floc);
	int n = 0;
	int d_org;
	double avg = 0;
	
	while(getline(f, line, '\n')){
		vector<int> original;
		splitStringintoInt(line, ' ', original);
		d_org = original.size();
		int _d = 0;
		avg += original.back();

		for(int i = 0; i < original.size()-1; i ++){
			map<int, pair<double, double> >::iterator itr = Numeric.find(i);
			if(itr != Numeric.end()){
				NX(_d, n) =2.0*(original[i] - itr->second.first)/(itr->second.second - itr->second.first)-1;
				_d += 1;    
			}
			else{
				map<int, int>::iterator  itr = Category.find(i);
				if(original[i] != 0){
					NX(_d+(original[i]-1),n) = 1;
				}
				_d += itr->second;
			} 
		}
		NX(NX.n_rows-1, n) = 1;
		YY_linear(n) = original.back();
		n++;
	}
	f.close();
	avg /= NX.n_cols;
	count = 0;
	for(int i = 0; i < NX.n_cols; i ++){
		if(YY_linear(i) >= avg)	YY_logistic(i) = 1; 
		else YY_logistic(i) = -1;
		map<int, pair<double, double> >::iterator itr = Numeric.find(d_org-1);
		YY_linear(i) = 2.0*(YY_linear(i) - itr->second.first)/(itr->second.second - itr->second.first)-1;
	}
}


// Algorithm 1 in the ICDE paper: Duchi et al.¡¯s Solution for One-Dimensional Numeric Data
double Duchi(double x, double eps) {
	if (x > 1) x = 1;
	if (x < -1) x = -1;
	double ne = exp(eps);
	double a = (x * (ne - 1) + ne + 1) / 2.0;
	double p = a / (ne + 1);
	int bit = 0;
	if (randu() <= p) bit = 1;
	double ce = (ne + 1) / (ne - 1);
	return (2*bit - 1)*ce ;
}


//Algorithm 2 in the ICDE paper: Piecewise Mechanism for One-Dimensional Numeric Data. 
double PM(double x, double eps) {
	double result;
	if (x > 1) x = 1;
	if (x < -1) x = -1;

	double z = exp(eps / 2);
	double P1 = (x + 1) / (2 + 2 * z);
	double P2 = z / (z + 1);
	double P3 = (1 - x) / (2 + 2 * z);

	double C = (z + 1) / (z - 1);
	double g1 = (C + 1)*x / 2 - (C - 1) / 2;
	double g2 = (C + 1)*x / 2 + (C - 1) / 2;

	double p = randu();
	if (p < P1) {// 
		result = -C + randu() * (g1 - (-C));
	}
	else if (p < P1 + P2) {
		result = (g2 - g1)*randu() + g1;
	}
	else {
		result = (C - g2)*randu() + g2;
	}
	return result;
}


//Hybrid mechanism
double HM(double x, double eps) {
	double result;
	if (eps < 0.61) {
		result = Duchi(x, eps);
	}
	else {
		double z = exp(-eps / 2);

		if (randu() <= z) {
			result = Duchi(x, eps);
		}
		else {
			result = PM(x, eps);
		}
	}
	return result;
}

/*
Each iteration asks a group of users to submit their gradients
*/
double MGD(string task, string method, vec& beta, vector<int> &train, vector<int>& test, mat &X, rowvec &Y_linear, rowvec &Y_logistic, float eps) {
	beta.zeros();
	int d = X.n_rows; int n = train.size();

	int bs =  int(d * log((float)d) / pow(eps, 2)); //batchsize, the number of users in one group
	if (bs < 100) bs = 100;
	
	int folds = n / bs;     // number of iterations

	int k = floor(eps / 2.5);//number of attributes each user submits
	if (k < 1) k = 1;
	if (k > d) k = d;

	vector<int> sampleDomain;
	for (int i = 0; i < d; i++)
		sampleDomain.push_back(i);               

	for (int i = 0; i < folds; i = i + 1) {//folds - nftest denotes the number the first 80% batch.
		vec sgrad = vec(d);
		sgrad.zeros();

		vec true_grad = vec(d);
		true_grad.zeros();
		
		iteration_num++;

		

		for (int id = 0; id < bs; id++) {//each user in one group
			

			int idd = train[i*bs + id];//the tuple id in the orginal data
			vec gradone = vec(d);//gradient of one user

			//compute gradients
			if (task == "linear") { Linear_grad(beta, X.unsafe_col(idd), Y_linear[idd], gradone); }
			else if (task == "log") { Logistic_grad(beta, X.unsafe_col(idd), Y_logistic[idd], gradone); }
			else { SVM_grad(beta, X.unsafe_col(idd), Y_logistic[idd], gradone); }

			gradone += beta * LAMBDA;//regularization term
			vector<int> currentSampleDomain = sampleDomain;
			for (int ki = 0; ki < k; ki++) {//choose k attributes to publish
				uniform_int_distribution<int> distribution(0, currentSampleDomain.size() - 1);
				int idxid = distribution(generator);
				int idx = currentSampleDomain[idxid];
				currentSampleDomain.erase(currentSampleDomain.begin() + idxid);
				double pai;
				if (method == "PM") pai = PM(gradone(idx), eps / k);//publish each value with privacy budget eps/k
				if (method == "HM") pai = HM(gradone(idx), eps / k);
				sgrad(idx) = sgrad(idx) + (double)d / k*pai;
				true_grad(idx) = true_grad(idx) + (double)d / k*gradone(idx);
			}

		}

		vec trueGrad = (true_grad / bs);
		vec grad = (sgrad / bs);
		beta = beta - 1 / sqrt(i + 1.0)*grad;//learing rate 1/sqrt(i+1)

		fstream tg;
		tg.open("true_grad.txt", ios::app);
		tg << trueGrad[0] << endl;
		tg.close();

		fstream ng;
		ng.open("noisy_grad.txt", ios::app);
		ng << grad[0] << endl;
		ng.close();

		fstream iterf;
		iterf.open("iteration_num.txt", ios::app);
		iterf << iteration_num << endl;
		iterf.close();

		fstream f;
		f.open("betao.txt", ios::app);
		f << beta[0] << endl;
		f.close();

		fstream ft;
		ft.open("betat.txt", ios::app);
		ft << beta[1] << endl;
		ft.close();


		double loss = 0;
		for (int j = 0; j < bs; j++) {
			int index= train[i*bs + j];
			if (task == "linear") {
				double value = dot(beta, X.unsafe_col(index));
				if (value < -1) value = -1;
				if (value > 1) value = 1;
				loss += pow(Y_linear[index] - value, 2);       // MSE
			}
		}
		loss = loss / bs;
		fstream fl;
		fl.open("loss.txt", ios::app);
		fl << loss << endl;
		fl.close();
	}
	//compute error
	double err = 0;
	for (int i = 0; i < test.size(); i++) {
		int idd = test[i];
		if (task == "linear") {
			double value = dot(beta, X.unsafe_col(idd));
			if(value < -1) value = -1;
			if(value > 1) value = 1;
			err += pow(Y_linear[idd] - value, 2);       // MSE
		}
		else {
			double z = dot(beta, X.unsafe_col(idd));
			double vl = -1.0;
			if (z >= 0) vl = 1.0;
			if (vl * Y_logistic[idd] < 0) err = err + 1;
		}
	}
	return err / test.size();
}



int main(){
	string task = "linear";//what kind of learning tasks (log, linear, svm) 

	//dataset
	string dataset = "mexico2000";
    int N = 4000000;//the number of tuples in the dataset

	//int d = 94;//after encoding the non-binary categorical attribute to binary, the number of dimensions
	int d = 94;

	/*string dataset = "BR";
	int N = 4283036;
	int d = 90;*/

	string fdomain =  dataset + ".domain";//path of the file for the attibute information: xx.domain
	string floc = dataset + ".dat";//path of the file for tuples: xx.dat
	
	vector<double> epsvector;
	epsvector.push_back(4);
	epsvector.push_back(2);
	epsvector.push_back(1);
	epsvector.push_back(0.5);
	
	int rep = 1; //each method runs rep times
	
	int crossnumber = 5;//crossnumber-cross validation

	mat NX(d, N); NX.zeros();	//tuples, each column is associated with one tuple
	rowvec YY_linear(N); YY_linear.zeros(); //predicted value for linear regressioin
	rowvec YY_logistic(N); YY_logistic.zeros();//predicted value for logistic regressioin and SVM
	loadData(fdomain, floc, NX, YY_linear, YY_logistic);//read data and transform non-binary categorical attributes to binary attributes
	for (int k = 0; k < 10000; k++){
		fstream f;
		f.open("data.txt", ios::app);
		f << NX(0,k) << " " << NX(1,k) << " " <<NX(93,k)<< endl;
		f.close();
	}
		
	std::random_device rd;
	auto seed = rd();
	arma_rng::set_seed( seed );  // set the seed to a random value
	generator.seed ( rd() );
	srand((unsigned)time(NULL));

	vector<int> id;                          
	for (int i = 0; i < N; i++) id.push_back(i);   

	for (int j = 0; j < epsvector.size(); j++) {
		double eps = epsvector[j];
		fstream epsf;
		epsf.open("eps.txt", ios::app);
		epsf << epsvector[j]<< endl;
		epsf.close();

		double err_PM = 0; double err_HM = 0;
		for(int repet = 0; repet < rep; repet ++){
			vector<int> idx = id;
			std::random_shuffle ( idx.begin(),idx.end() );//shuffle our data set
			int onesize = floor((double)idx.size()/ crossnumber);// the size of data for testing

			double error_each_PM = 0; double error_each_HM = 0;
			for(int i = 0; i < crossnumber; i ++){//do crossnumber-cross validation
				
				//get train data and test data
				vector<int> trainidx; vector<int> testidx;
				for(int m = 0; m < crossnumber; m ++){
					if(m == i){
						testidx.insert(testidx.end(), idx.begin() + onesize * m, idx.begin() + onesize * (m + 1));
					}
					else
						trainidx.insert(trainidx.end(), idx.begin() + onesize * m, idx.begin() + onesize * (m + 1));
				}
				//do stochastic gradient descent
				vec beta(d); //parameter
				double result_PM = MGD(task, "PM", beta, trainidx, testidx, NX, YY_linear, YY_logistic, eps);//use "piecewise mechanism" to publish gradients
				double result_HM = MGD(task, "HM", beta, trainidx, testidx, NX, YY_linear, YY_logistic, eps);//use "Hybrid mechanism" to publish gradients
				err_PM += result_PM;
				err_HM += result_HM;
			}
		}//eps
		fstream f;
		f.open("SGD.txt", ios::app);
		f <<dataset<<" "<<epsvector[j] << " "<< err_PM / (crossnumber * rep)<<" "<< err_HM / (crossnumber * rep)<< endl;
		f.close();		
	}
	//system("pause");
	return 0;
}


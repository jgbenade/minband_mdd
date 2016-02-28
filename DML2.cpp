//============================================================================
// Name        : dml.cpp
// Author      : Nam Ho-Nguyen
// Version     :
// Copyright   :
// Description : Distance Metric Learning in C++
//============================================================================

#include <iostream>
#include <armadillo>
#include <Eigen/Dense>

using namespace std;
using namespace arma;
using namespace Eigen;


mat std_vec_to_mat(vector<vector<double>> Xstd) {
	int N = Xstd.size(); int d = Xstd[0].size();
	// Flatten your Xstd into Xflat
	vector<double> Xflat;
	for (auto vec : Xstd) {
		for (auto el : vec) {
			Xflat.push_back(el);
		}
	}
	// Create your Armadillo matrix X
	mat X2 = conv_to<mat>::from(Xflat);
	X2.reshape(d,N); mat X = X2.t();
	return X;
}


vector<vector<double>> mat_to_std_vec(mat &A) {
	vector<vector<double>> V(A.n_rows);
    for (size_t i = 0; i < A.n_rows; ++i) {
        V[i] = conv_to< vector<double> >::from(A.row(i));
    };
    return V;
}


vector<vector<int>> mat_to_std_vec_int(imat &A) {
	vector<vector<int>> V(A.n_rows);
    for (size_t i = 0; i < A.n_rows; ++i) {
        V[i] = conv_to< vector<int> >::from(A.row(i));
    };
    return V;
}


MatrixXd arma_to_eigen_mat(mat A) {
	int m = A.n_rows; int n = A.n_cols;
	vector<vector<double>> V = mat_to_std_vec(A);
	MatrixXd B(m,n);
	for(int i = 0; i < m; i++) {
		B.row(i) = Map<VectorXd>(V[i].data(), n);
	}
	return B;
}


mat eigen_to_arma_mat(MatrixXd A) {
	int m = A.rows(); int n = A.cols();
	vector<double> W(A.data(), A.data() + m*n);
	mat B = conv_to<mat>::from(W);
	B.reshape(m,n);
	return B;
}


void my_eig(mat& eigval, mat& eigvec, mat A) {
	int m = A.n_rows; int n = A.n_cols;
	vec eigvalvec;
	eig_sym(eigvalvec, eigvec, A);
	eigval = zeros<mat>(m,n);
	eigval.diag() = eigvalvec;


	/*MatrixXd B = arma_to_eigen_mat(A);

	SelfAdjointEigenSolver<MatrixXd> eigs(B);
	MatrixXd Dz; Dz.setZero(m,n);
	//diagonal matrix of eigenvalues
	Dz.diagonal() = eigs.eigenvalues();
	eigval = eigen_to_arma_mat(Dz);
	//matrix of eigenvectors
	MatrixXd E = eigs.eigenvectors();
	eigvec = eigen_to_arma_mat(E);*/
}


double objective(mat& X, vector<vector<int> >& D, mat& A) {
	int N = X.n_rows;
	int d = X.n_cols;
	double sum_norms = 0;

	int nd = D[0].size(); int i; int j;
	for(int k=0; k < nd; k=k+1) {
		i = D[0][k]; j = D[1][k];
		vec d_ij = (X.row(i) - X.row(j)).t();
		double norm_ij = dot(d_ij, A * d_ij);
		norm_ij = sqrt(norm_ij);
		sum_norms = sum_norms + norm_ij;
	}
	return sum_norms;
}


mat gradD(mat& X, vector<vector<int> >& D, mat& A) {
	int N = X.n_rows;
	int d = X.n_cols;

	double fudge = 0.0000001;
	mat G = zeros<mat>(d,d); double sum_norms = 0;
	int nd = D[0].size(); int i; int j;
	for(int k=0; k < nd; k=k+1) {
		i = D[0][k]; j = D[1][k];
		vec d_ij = (X.row(i) - X.row(j)).t();
		double norm_ij = dot(d_ij, A * d_ij);
		norm_ij = sqrt(norm_ij);
		sum_norms = sum_norms + norm_ij;
		mat M_ij = d_ij * d_ij.t();
		G = G + 0.5*M_ij/(norm_ij + fudge);
	}

	G = G/sum_norms;

	return G;
}


mat gradProj(vec& g1, mat& g2, int d) {

	g2.reshape(d*d,1); //cout << "g2: " << g2 << endl;
	g2 = normalise(g2); //cout << "g2n: " << g2 << endl;

	vec g = g1 - dot(g1,g2)*g2; //cout << "g: " << g << endl;
	g = normalise(g); //cout << "gn: " << g << endl;

	mat gm = conv_to<mat>::from(g);
	gm.reshape(d,d); //cout << "gmat: " << gm << endl;

	return gm;
}


mat iter_projection(mat X, vector<vector<int> > S, vector<vector<int> > D, mat A, vec w, double t, int maxiter) {

	int N = X.n_rows;
	int d = X.n_cols;
	double epsilon = 0.01;
	double threshold = 0.01;
	int maxcount = 100;

	vec w1 = normalise(w);
	double t1 = t/norm(w);

	int count = 1;
	double alpha = 0.1;

	vec grad1 = w;
	mat grad2 = gradD(X,D,A);
	mat M = gradProj(grad1, grad2, d);

	mat A_last = A;
	bool done = 0;
	bool converged = 0;
	while(done==0) {
		//projection onto domain
		bool satisfy = 0;
		int projection_iters = 0;
		while((projection_iters<maxiter)&&(satisfy==0)) {
			mat A0 = A;
			mat x0 = A0; x0.reshape(d*d,1);
			if(dot(w,x0) <= t) {
				A = A0;
			}
			else {
				vec xm = x0 + (t1 - dot(w1,x0))*w1;
				mat x = conv_to<mat>::from(xm);
				x.reshape(d,d);
				A = x;
			}

			A = 0.5*(A + A.t());
			mat D1; mat P;
			my_eig(D1, P, A);
			D1 = max(D1, zeros<mat>(d,d));
			A = P * D1 * P.t();
			mat Atest = A; Atest.reshape(d*d,1);
			double error = (dot(Atest,w) - t)/t;
			if(error > epsilon) {
				satisfy = 0;
			}
			else {
				satisfy = 1;
			}
			projection_iters += 1;
		}

		double obj_prev = objective(X, D, A_last);
		double obj = objective(X, D, A);
		if( ((obj > obj_prev) || (count==1)) && (satisfy==1) ) {
			alpha = alpha * 1.05; A_last = A;
			grad1 = w;
			grad2 = gradD(X, D, A);
			M = gradProj(grad1, grad2, d);
			A = A + alpha*M;
		}
		else {
			alpha = alpha/2;
			A = A_last + alpha*M;
		}

		mat Awow = 0.5*(A + A.t());
		mat Dwow; mat Pwow;
		my_eig(Dwow, Pwow, Awow);
		Dwow = max(Dwow, zeros<mat>(d,d));
		Awow = Pwow * Dwow * Pwow.t();

		double delta = norm(alpha*M,"fro")/norm(A_last,"fro");
		count = count + 1;
		if((count==maxcount)||(delta<threshold)) {
			done = 1;
			if(delta > threshold) {
				converged = 0;
			}
			else {
				converged = 1;
			}
		}
		else {
			done = 0;
		}
	}

    cout << "iter_projection converged: " << converged <<endl;

	A = A_last;

	return A_last;
}

mat opt(mat X, vector<vector<int> > S, vector<vector<int> > D, int maxiter) {
	//mat X = std_vec_to_mat(Xstd);

	int N = X.n_rows;
	int d = X.n_cols;

	mat A = eye(d,d)*0.1;

	if (D[1].size() == 0) {
		cout << "no dissimilar pairs return identity"<< endl;
		return A;
	}

	mat W = zeros<mat>(d,d);
	int ns = S[0].size(); int i; int j;
	for(int k=0; k < ns; k=k+1) {
		i = S[0][k]; j = S[1][k];
		vec d_ij = (X.row(i) - X.row(j)).t();
		W = W + d_ij*d_ij.t();
	}

	mat Wnew = W; mat Anew = A;
	Wnew.reshape(d*d,1); Anew.reshape(d*d,1);
	vec w = Wnew;
	vec avec = Anew/100;
	double t = dot(w, avec);

	A = iter_projection(X, S, D, A, w, t, maxiter);

	return A;
}


/*
int main() {
	int N = 10; int d = 4;
	mat A = randi<mat>(d,d, distr_param(0,20));

	A = 0.5*(A + A.t());
	mat E; mat P;
	my_eig(E, P, A);
	E = max(E, zeros<mat>(d,d));
	A = P * E * P.t();

	imat Sm = randi(N, N, distr_param(0,1));
	imat Dm = ones<imat>(N,N);
	Dm = Dm - Sm;
	Dm = (Dm + Dm.t())/2;
	Dm.diag() = zeros<ivec>(N);
	Sm = (Sm + Sm.t())/2;
	Sm.diag() = zeros<ivec>(N);
	int nd = accu(Dm)/2; int ns = accu(Sm)/2;
	imat DS = zeros<imat>(nd,2);
	int k=0;
	for(int i=0; i < N; i=i+1) {
		for(int j=0; j < i; j=j+1) {
			if(Dm(i,j)==1) {
				DS(k,0) = i;
				DS(k,1) = j;
				k = k + 1;
			}
		}
	}
	imat SS = zeros<imat>(ns,2);
	k=0;
	for(int i=0; i < N; i=i+1) {
		for(int j=0; j < i; j=j+1) {
			if(Sm(i,j)==1) {
				SS(k,0) = i;
				SS(k,1) = j;
				k = k + 1;
			}
		}
	}
	vector<vector<int> > S = mat_to_std_vec_int(SS);
	vector<vector<int> > D = mat_to_std_vec_int(DS);

	mat X = randu<mat>(N,d);
	vector<vector<double> > Xstd = mat_to_std_vec(X);


	mat F = opt(X, S, D, 1000);
	cout << "opt: " << endl << F << endl;

	mat eigvalF; mat eigvecF;
	my_eig(eigvalF, eigvecF, F);
	cout << "eigvalF: " << endl << eigvalF.diag()  << endl;
	cout << "eigvecF: " << endl << eigvecF << endl;


	return 0;
}*/







/*--testing--
double f = objective(X, D, A);
cout << "objective: " << f << endl;

mat F1 = gradD(X,D,A);
cout << "gradD: " << F1 << endl;

mat Anew = randi<mat>(d,d, distr_param(0,10));
Anew = 0.5*(Anew + Anew.t());
my_eig(E, P, Anew);
//cout << E << endl;
E = max(E, zeros<mat>(3,3));
//cout << E << endl;
Anew = P * E * P.t();
Anew.reshape(d*d,1);
vec w = Anew; //mat w2 = conv_to<mat>::from(w);
//w2.reshape(d,d);
//cout << w << w2 << endl;
mat F2 = gradProj(w, A, d);
cout << "gradProj: " << F2 << endl;

mat F3 = iter_projection(X, S, D, A, w, 1, 1000);
cout << "iter_projection: " << F3 << endl;
------*/

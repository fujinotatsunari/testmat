#include"matrix.hpp"
#include<iostream>
#include<vector>
#include <random>

int random(int low, int high)
{
	return low + rand() % (high - low + 1);
}

int main(void) {
	int dim = 10;
	vector<double> bvec;
	vector<double> xvec;
	bvec.resize(dim);
	xvec.resize(dim);
	cout << "bvec = [" << endl;
	for (int i = 0; i < dim; i++) {
		bvec[i] = (double)random(1, 9);
		cout << bvec[i] << endl;
	}
	cout << "]" << endl;


	Matrix mat(dim, dim);

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (i == j) {
				mat[i][j] = 1;
			}
			else {
				mat[i][j] = random(0, 3);
			}

			
		}
	}
	mat.print();

	LU_solve(mat, xvec, bvec);

	cout << "after solved" << endl;
	mat.print();
	cout << "xvec = [" << endl;
	for (int i = 0; i < dim; i++) {
		cout << xvec[i] << endl;
	}
	cout << "]" << endl;


	return 0;

}
#pragma once
#include"bicg.hpp"
#include"matrix.hpp"
#include<iostream>
#include<vector>
using namespace std;
//bicg-stab法class設計部分
//必要そうなもの
//0:bicg法の実装に必要そうなobject

//1:bicg法による行列解法のglobalな関数 (正味これだけでもええかも)
//LU分解はmatrix.hppに実装してます.わかんないところがあれば藤野まで

//以下classの宣言 設計部分


//以下行列solver実装部分
// void hoge_solver( Matrix& A, vector<double>& x, vector<double>& b);
//の形式で宣言　設計してください
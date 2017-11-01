#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/GF2.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>

NTL_CLIENT

using namespace std;

int main() {

	/*****************************************************/
	//Opening the file which contains the hash value.
	/*****************************************************/
	
	FILE *fptr;
	char filename[15];
	char ch;

	int i, j;

	vec_GF2 h[8];
	for(i=0; i<8; ++i) {
		h[i].SetLength(64);
//		random(h[i], 64);	//change this
	}	

	cout<<"Enter the filename to be opened : ";
	cin>>filename;

	ifstream infile; 
   	infile.open(filename);
	if (!infile)
	{
		printf("Cannot open file \n");
		exit(0);
	}

	i = 0;
	infile.get(ch);
	while (i < 512)
	{
		h[i/64][i%64] = ch - '0';
		++i;
		infile.get(ch);
	}

	infile.close();
	
	/*****************************************************/
	//End of reading file
	/*****************************************************/

	h[0][0] = h[0][0] + 1; // iota inverse

	vec_GF2 hh[9];
	for(i = 0; i<9; ++i)
		hh[i].SetLength(64);

	for (i = 0; i < 64; ++i) // chi inverse
	{
		hh[8][i] = 1;
		hh[7][i] = h[7][i];
		hh[6][i] = (h[6][i]+h[7][i]+1);
		hh[5][i] = (h[5][i]+(h[6][i]+h[7][i])*h[7][i]);
		
		hh[0][i] = (h[0][i]+(h[1][i]+1)*(h[2][i]+(h[3][i]+1)*h[4][i]));
		hh[1][i] = (h[1][i]+(h[2][i]+1)*(h[3][i]+(h[4][i]+1)*h[0][i]));
		hh[2][i] = (h[2][i]+(h[3][i]+1)*(h[4][i]+(h[0][i]+1)*h[1][i]));
		hh[3][i] = (h[3][i]+(h[4][i]+1)*(h[0][i]+(h[1][i]+1)*h[2][i]));
		hh[4][i] = (h[4][i]+(h[0][i]+1)*(h[1][i]+(h[2][i]+1)*h[3][i]));
	}

	h[0][0] = h[0][0] + 1;	//to get back the hash value

	vec_GF2 d[9], e[9];
	for(i = 0; i < 9; ++i) {
		d[i].SetLength(64);
		e[i].SetLength(64);
	}

	random(e[7], 64);
	random(e[8], 64);
	e[8][63] = 1;

	for(i = 0; i < 64; ++i) {
		d[0][i] = 0;		//d_0 = 0
		d[2][i] = 0; 		//d_2 = 0
		d[4][i] = 0; 		//d_4 = 0
		d[3][(i + 5) % 64] = hh[4][i] + hh[6][(i + 58) % 64];	//d_3(5) = hh_4(0) + hh_6(58)
		e[3][(i + 28) % 64] = hh[5][i] + hh[3][(i + 7) % 64];	//e_3(28) = hh_5(0) + hh_3(7)
		e[6][(i + 44) % 64] = hh[1][i] + 1;	//e_6(44) = hh_1(0) + 1
	}

	/*****************************************************/
	//Building a system of linear equation. Ax = b where x is a vector containing lane variables.
	/*****************************************************/

	vec_GF2 b, x, RC;
	b.SetLength(384);
	x.SetLength(384);
	RC.SetLength(64);
	RC[0] = 1;

	for(i = 0; i<64; ++i) {
		b[(0 * 64) + i] = hh[0][i] + RC[i] + d[3][(i + 55) % 64] + d[3][(i + 54) % 64] + e[6][(i + 63) % 64];

		b[(1 * 64) + i] = hh[2][i] + d[3][(i + 34) % 64] + e[6][(i + 43) % 64] + d[3][(i + 6) % 64] + e[3][(i + 42) % 64] + e[8][(i + 42) % 64];

		b[(2 * 64) + i] = hh[3][i] + e[7][(i + 21) % 64] + d[3][(i + 11) % 64];

		b[(3 * 64) + i] = hh[6][i] + d[3][(i + 48) % 64] + e[8][(i + 20) % 64] + d[3][(i + 47) % 64] + RC[(i + 19) % 64] + e[3][(i + 20) % 64];

		b[(4 * 64) + i] = hh[7][i] + d[3][(i + 58) % 64] + d[3][(i + 57) % 64] + e[6][(i + 2) % 64];
	
		b[(5 * 64) + i] = 1 + RC[(i +45) % 64] + d[3][(i + 9) % 64] + e[7][(i + 44) % 64];
		
	}

	mat_GF2 A;
	A.SetDims(384, 384);
	clear(A);


	for(i=0; i<64; ++i) {
		//First equation
		A[0 + i][(0 * 64) + ((i + 44) % 64)] = 1;	//d_1(44)
		A[0 + i][(0 * 64) + ((i + 43) % 64)] = 1;	//d_1(43)
		A[0 + i][(1 * 64) + ((i + 0) % 64)] = 1;	//e_0(0)
		A[0 + i][(2 * 64) + ((i + 63) % 64)] = 1;	//e_1(63)
		A[0 + i][(4 * 64) + ((i + 0) % 64)] = 1;	//e_4(0)

		//Second equation
		A[64 + i][(0 * 64) + ((i + 23) % 64)] = 1;	//d_1(23)
		A[64 + i][(0 * 64) + ((i + 43) % 64)] = 1;	//d_1(43)
		A[64 + i][(2 * 64) + ((i + 43) % 64)] = 1;	//e_1(43)

		//Third equation
		A[128 + i][(0 * 64) + ((i + 0) % 64)] = 1;	//d_1(0)
		A[128 + i][(3 * 64) + ((i + 21) % 64)] = 1;	//e_2(21)
		A[128 + i][(4 * 64) + ((i + 20) % 64)] = 1;	//e_4(20)

		//Fourth equation
		A[192 + i][(0 * 64) + ((i + 21) % 64)] = 1;	//d_1(21)
		A[192 + i][(0 * 64) + ((i + 20) % 64)] = 1;	//d_1(20)
		A[192 + i][(1 * 64) + ((i + 19) % 64)] = 1;	//e_0(19)
		A[192 + i][(5 * 64) + ((i + 19) % 64)] = 1;	//e_5(19)

		//Fifth equation
		A[256 + i][(0 * 64) + ((i + 4) % 64)] = 1;	//d_1(4)
		A[256 + i][(0 * 64) + ((i + 47) % 64)] = 1;	//d_1(47)
		A[256 + i][(0 * 64) + ((i + 46) % 64)] = 1;	//d_1(46)
		A[256 + i][(2 * 64) + ((i + 2) % 64)] = 1;	//e_1(2)
		A[256 + i][(4 * 64) + ((i + 3) % 64)] = 1;	//e_4(3)

		//Sixth equation
		A[320 + i][(0 * 64) + ((i + 46) % 64)] = 1;	//d_1(46)
		A[320 + i][(1 * 64) + ((i + 45) % 64)] = 1;	//e_0(45)
		A[320 + i][(3 * 64) + ((i + 44) % 64)] = 1;	//e_2(44)
		A[320 + i][(5 * 64) + ((i + 45) % 64)] = 1;	//e_5(45)
	}

	/*****************************************************/
	//End of building linear equations
	/*****************************************************/

	GF2 det;
	solve(det, A, x, b);	//solving the system of linear equation

	/*****************************************************/
	//Extracting the solution of the system of linear equation
	/*****************************************************/

	for(i = 0; i<64; ++i) {
		d[1][i] = x[(0 * 64) + i];
		e[0][i] = x[(1 * 64) + i];
		e[1][i] = x[(2 * 64) + i];
		e[2][i] = x[(3 * 64) + i];
		e[4][i] = x[(4 * 64) + i];
		e[5][i] = x[(5 * 64) + i];
	}

	/*****************************************************/
	//End
	/*****************************************************/

	d[5] = d[0];
	d[6] = d[1];
	d[7] = d[2];
	d[8] = d[3];

	cout<<"****************HASH VALUE**************\n";
	for(i=0; i<8; ++i) {
		for(j=0; j<64; ++j) {
			cout<<h[i][j];
		}	
//		cout<<endl;
	}
	cout<<"****************************************\n\n";

	cout<<"The preimage is stored in the file named \"preimage\".\n";
	ofstream outfile;
	outfile.open("preimage", ios::out);
	if(!outfile)
	{
		cout<<"Cannot open file \n";
		exit(1);
	}	

	cout<<"*****************MESSAGE D***************\n";

	for(i=0; i<9; ++i) {
		for(j=0; j<64; ++j) {
			cout<<d[i][j];
			outfile<<d[i][j];
		}	
//		cout<<endl;
	}
	cout<<"\n****************************************\n\n";

	cout<<"*****************MESSAGE E***************\n";
	for(i=0; i<9; ++i) {
		for(j=0; j<64; ++j) {
			cout<<e[i][j];
			outfile<<e[i][j];
		}	
//		cout<<endl;
	}
	cout<<"\n****************************************\n\n";

	outfile.close();
	
}
/*
* Author : Shou-Fu Lo
* Last Modified : 2017/07/31
*
* This library includes the LCS Algorithms as following:
* Our Algorithm_MLCS
* Rahman and Rahman Algorithm_MLCS
* Deorowicz and Danek Algorithm_MLCS
* Peng et al. Algorithm_MLCS
* Huang et al. Algorithm_MLCS
*
* Copyright (c) 2017, Shou-Fu Lo
* All rights reserved.
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Shou-Fu Lo may not be used to endorse or promote products derived from this
*       software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY Shou-Fu Lo "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
* FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
*  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <utility>
#include "time.h"
#include <queue>
#include <functional>
#include <list>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <windows.h>
using namespace std;

int MLCS_DP(vector<int>T, vector<int>A, vector<int>B, int Tlen, int Alen, int Blen);
int MLCS_Ours(vector<int> T, vector<int> A, vector<int>B, int Tlen, int Alen, int Blen, int size);
void Dominate(int *L1x, int *L1y, int *L2x, int *L2y, int *Rx, int *Ry, int Alenplusone, int Blenplusone);
int MLCS_Bitpar(vector<int> T, vector<int> A, vector<int> B, int Tlen, int Alen, int Blen, int size, int WordSize);
int MLCS_PengSparse(vector<int> T, vector<int> A, vector<int>B, int Tlen, int Alen, int Blen, int size);
int MLCS_SPDP_Boundedheap(vector<int>T, vector<int>A, vector<int>B, int Tlen, int Alen, int Blen, int size);
int binarysearch(vector< pair<int, int> > P, int i, int k);
class triple {
public:
	int x;
	int y;
	//int z;
	//triple *parent;
	//triple(int posA, int posB, int posT){ x = posA; y = posB; z = posT; parent = nullptr; }
	//triple(int posA, int posB, int posT){ x = posA; y = posB; z = posT; }
	triple(int posA, int posB) { x = posA; y = posB; }
};
int main()
{
	int r, alphaSize, Case;
	float m ;

    // inputLength
    cout<<"Tlen= {1000, 2000, 5000, 10000},"<<endl<< "Ratio = |A| / |T| = {0.1, 0.2, 0.3, 0.4, 0.5}"<<endl<<"Alphabet sizes {4, 64, 1000},"<<endl<<endl;
    cout<<"Input the length of string T: ";
    cin>>r;
    cout<<"Input the ratio: ";
    cin>>m;
    cout<<"Alphabet sizes : ";
    cin>>alphaSize;
    cout<<"The number of executions(1,2,...,100) : ";
    cin>>Case;

    int Tlen, Alen, Blen, i, j, k, p, similar, C;
	Tlen = r, Alen = r * m, Blen = r - Alen;

	vector<int> A(Alen + 2);
	vector<int> B(Blen + 2);
	vector<int> T(Tlen + 2);

	if(alphaSize == 4)
    {
        similar=80;
        p=2;
    }
    else if (alphaSize == 64)
    {
        similar=30;
        p=5;
    }
    else
    {
        similar=10;
        p=5;
    }

	for(int S = similar; S <= 100; S+=p)
    {
        string str, data;
        stringstream ss;
        ss << "table_" << Tlen << "_" << Alen << "_" << Blen << "_" << alphaSize << "_" << S << ".csv";
        string filename;
        ss >> filename;
        ifstream fptr(filename.c_str());
        if (!fptr)
        {
            cout << filename << " not exist" << endl;
            continue;
        }
        else
        {
            cout << filename << " is running" << endl;
        }

        C = Case;
        getline(fptr, str, '\n');
        while (C--)
        {
            getline(fptr, str, ',');
            getline(fptr, str, ',');
            istringstream templine(str);
            getline(templine, data, '|'); //"-1"
            k = 1;
            while (getline(templine, data, '|'))
            {
                T[k] = atoi(data.c_str());
                k++;
            }
            getline(fptr, str, ',');
            istringstream templine2(str);
            getline(templine2, data, '|'); //"-2"
            i = 1;
            while (getline(templine2, data, '|'))
            {
                A[i] = atoi(data.c_str());
                i++;
            }

            getline(fptr, str, '\n');
            istringstream templine3(str);
            getline(templine3, data, '|'); //"-3"
            j = 1;
            while (getline(templine3, data, '|'))
            {
                B[j] = atoi(data.c_str());
                j++;
            }
            ss.clear();

            cout << "----------MLCS---------" << endl;
            cout<<"Tlen = "<<Tlen<<", Alen = "<<Alen<<", Blen = "<<Blen<<", Alphabetsize = "<<alphaSize<<", Similar = " <<S<<"%"<< endl;
            cout << "MLCS_Ours               : " << MLCS_Ours(T, A, B, Tlen, Alen, Blen, alphaSize) << endl;
            cout << "MLCS_Bitpar             : " << MLCS_Bitpar(T, A, B, Tlen, Alen, Blen, alphaSize, 64) << endl;
            cout << "MLCS_DP                 : " << MLCS_DP(T, A, B, Tlen, Alen, Blen) << endl;
            cout << "MLCS_PengSparse         : " << MLCS_PengSparse(T, A, B, Tlen, Alen, Blen, alphaSize) << endl;
            cout << "MLCS_SPDP_Boundedheap   : " << MLCS_SPDP_Boundedheap(T, A, B, Tlen, Alen, Blen, alphaSize) << endl;



        }
	}
	return 0;
}

int MLCS_DP(vector<int>T, vector<int>A, vector<int>B, int Tlen, int Alen, int Blen)
{
	vector<vector<int> > AT(Alen + 1, vector<int>(Tlen + 1));
	vector<vector<int> > BT(Blen + 1, vector<int>(Tlen + 1));
	vector<vector<int> > ABT1(Alen + 1, vector<int>(Blen + 1));
	vector<vector<int> > ABT2(Alen + 1, vector<int>(Blen + 1));

	AT[0][0] = 0;
	for (int i = 1; i <= Alen; i++)
		AT[i][0] = 0;
	for (int i = 1; i <= Tlen; i++)
		AT[0][i] = 0;
	for (int i = 1; i <= Alen; i++)
		for (int j = 1; j <= Tlen; j++)
		{
			if (A[i] == T[j])
				AT[i][j] = AT[i - 1][j - 1] + 1;
			else
				AT[i][j] = max(AT[i][j - 1], AT[i - 1][j]);
		}
	BT[0][0] = 0;
	for (int i = 1; i <= Blen; i++)
		BT[i][0] = 0;
	for (int i = 1; i <= Tlen; i++)
		BT[0][i] = 0;
	for (int i = 1; i <= Blen; i++)
		for (int j = 1; j <= Tlen; j++)
		{
			if (B[i] == T[j])
				BT[i][j] = BT[i - 1][j - 1] + 1;
			else
				BT[i][j] = max(BT[i][j - 1], BT[i - 1][j]);
		}

	int layer = 1;
	while (1)
	{
		if (layer % 2 == 1) {
			for (int i = 1; i <= Alen; i++) {
				ABT1[i][0] = AT[i][layer];
				for (int j = 1; j <= Blen; j++) {
					ABT1[0][j] = BT[j][layer];
					ABT1[i][j] = max(max(ABT1[i - 1][j], ABT1[i][j - 1]), ABT2[i][j]);
					if (A[i] == T[layer])
					{
						ABT1[i][j] = max(ABT1[i][j], ABT2[i - 1][j] + 1);

					}
					if (B[j] == T[layer])
					{
						ABT1[i][j] = max(ABT1[i][j], ABT2[i][j - 1] + 1);
					}
				}
			}
			if (layer == Tlen)
				return ABT1[Alen][Blen];
			++layer;
		}
		else {
			for (int i = 1; i <= Alen; i++) {
				ABT2[i][0] = AT[i][layer];
				for (int j = 1; j <= Blen; j++) {
					ABT2[0][j] = BT[j][layer];
					ABT2[i][j] = max(max(ABT2[i - 1][j], ABT2[i][j - 1]), ABT1[i][j]);
					if (A[i] == T[layer])
					{
						ABT2[i][j] = max(ABT2[i][j], ABT1[i - 1][j] + 1);

					}
					if (B[j] == T[layer])
					{
						ABT2[i][j] = max(ABT2[i][j], ABT1[i][j - 1] + 1);
					}
				}
			}
			if (layer == Tlen)
				return ABT2[Alen][Blen];
			++layer;
		}
	}
}

int MLCS_Ours(vector<int> T, vector<int> A, vector<int>B, int Tlen, int Alen, int Blen, int size)
{
	vector<vector<int> > nextA(size + 1, vector<int>(Alen + 1, 0));
	vector<vector<int> > nextB(size + 1, vector<int>(Blen + 1, 0));
	int L = 0;

	int Alenplusone = Alen + 1;
	int Blenplusone = Blen + 1;
	int Lmaxsize = (Alen + Blen < Tlen ? Alen + Blen : Tlen); // Max LCS length

	int Dmaxsize = (Alen < Blen ? Blen : Alen); // Max size of a dominating set

	int **Dx = new int *[Lmaxsize + 2];
	int **Dy = new int *[Lmaxsize + 2];
	for (int i = 0; i < Lmaxsize + 2; ++i) {
		Dx[i] = new int[Dmaxsize + 2];
		Dy[i] = new int[Dmaxsize + 2];
	}
	int *WAx = new int[Dmaxsize + 2];
	int *WAy = new int[Dmaxsize + 2];
	int *WBx = new int[Dmaxsize + 2];
	int *WBy = new int[Dmaxsize + 2];
	int *tempDx = new int[Dmaxsize + 2];
	int *tempDy = new int[Dmaxsize + 2];
	// Dmaxsize=smaller of Alen and Blen, cannot set to Tlenvector<vector<int> > Dx(Lmaxsize + 2, vector<int>(Dmaxsize, -1));
	/*vector<vector<int> > Dx(Lmaxsize + 2, vector<int>(Dmaxsize + 2, -1));
	vector<vector<int> > Dy(Lmaxsize + 2, vector<int>(Dmaxsize + 2, -1));
	vector<int> WAx(Dmaxsize + 2, -1);
	vector<int> WAy(Dmaxsize + 2, -1);
	vector<int> WBx(Dmaxsize + 2, -1);
	vector<int> WBy(Dmaxsize + 2, -1);
	vector<int> tempDx(Dmaxsize + 2, -1);
	vector<int> tempDy(Dmaxsize + 2, -1);*/



	//construct array of nextA
	for (int i = 1; i <= size; i++) {
		nextA[i][Alen] = Alen + 1;
		nextB[i][Blen] = Blen + 1;
	}
	for (int j = Alen; j >= 1; j--) {
		for (int i = 1; i <= size; i++)
			nextA[i][j - 1] = nextA[i][j];
		nextA[A[j]][j - 1] = j;  //  need update only if (A[j] == i)
	}

	//construct array of nextB
	for (int j = Blen; j >= 1; j--) {
		for (int i = 1; i <= size; i++)
			nextB[i][j - 1] = nextB[i][j];
		nextB[B[j]][j - 1] = j;  //  need update only if (B[j] == i)
	}


	Dx[0][1] = 0;
	Dy[0][1] = 0;
	Dx[0][2] = Alenplusone;
	Dy[0][2] = 0;

	for (int s = 1; s <= Lmaxsize; ++s) {
		Dx[s][1] = Alenplusone;
		Dy[s][1] = 0;
	}

	int PosT = 0;
	for (int i = 1; i <= Tlen; ++i) {
		if (i > Tlen - L)
			break;
		for (int s = 1; s <= Tlen - i + 1; ++s) {
			PosT = i + s - 1;
			//(Dk−1,s−1) Extension of A, extension of B
			int j, k;
			int sMinusOne = s - 1;
			for (j = 1, k = 1; Dx[sMinusOne][j] != Alenplusone; ++j) {
				WAx[j] = nextA[T[PosT]][Dx[sMinusOne][j]];
				WAy[j] = Dy[sMinusOne][j];
				WBx[j] = Dx[sMinusOne][j];   // some <x, Blenplusone> are added
				WBy[j] = nextB[T[PosT]][Dy[sMinusOne][j]];
			}
			WAx[j] = Alenplusone;
			WBx[j] = Alenplusone;

			tempDx[1] = Alenplusone;

			// Dominate(Dk−1,s ∪ ExtA )
			if ((Dx[s][1] != Alenplusone) || (WAx[1] != Alenplusone)) {
				Dominate(Dx[s], Dy[s], WAx, WAy, tempDx, tempDy, Alenplusone, Blenplusone);
			}
			// Dominate(Dominate(Dk−1,s ∪ ExtA)∪ExtB)
			if ((tempDx[1] != Alenplusone) || (WBx[1] != Alenplusone)) {
				Dominate(tempDx, tempDy, WBx, WBy, Dx[s], Dy[s], Alenplusone, Blenplusone);
			}


			if ((Dx[s][1] == Alenplusone))
				break;

			if (s > L)
			{
				L = s;
			}

		}//end for j
	}//end for i

	 // delete the allocated array Dx[][], Dy[][]
	for (int i = 0; i < Lmaxsize + 2; ++i) {
		delete[] Dx[i];
		delete[] Dy[i];
	}
	delete[] Dx; 	delete[] Dy;
	delete[] WAx; 	delete[] WAy; 	delete[] WBx;	delete[] WBy;
	delete[] tempDx;	delete[] tempDy;
	return L;
}

void Dominate(int *L1x, int *L1y, int *L2x, int *L2y, int *Rx, int *Ry, int Alenplusone, int Blenplusone) {

	int i = 1, j = 1, tempx, tempy;
	int Rsizeminusone = 0;   // compare starting from [0]
							 //cout<<L2x[1]<<" "<<L2y[1]<<endl;

	Rx[0] = -1;
	Ry[0] = Blenplusone;   // [0] used for discard <x, Blenplusone>
	Rx[1] = Alenplusone;
	Ry[1] = Blenplusone;
	while ((L1x[i] != Alenplusone) || (L2x[j] != Alenplusone)) {
		//<is , js> ← the non-null smaller one of the ﬁrst 2-tuples of L1 and L2
		if (L1x[i] <= L2x[j]) { //<3 , 2> is smaller than <4 , 1>
								//	tempx = L1x[i];
								//	tempy = L1y[i];
			if (L1y[i] >= Ry[Rsizeminusone])  // && L2x[j] >= Rx[Rsizeminusone]
				;   // no operation, discard L2
			else if (L1x[i] > Rx[Rsizeminusone] && L1y[i] < Ry[Rsizeminusone]) {
				++Rsizeminusone;
				Rx[Rsizeminusone] = L1x[i];
				Ry[Rsizeminusone] = L1y[i];
			}
			else {
				Rx[Rsizeminusone] = L1x[i];
				Ry[Rsizeminusone] = L1y[i];
			} // else discard L1
			++i;
		}
		else {
			//	tempx = L2x[j];
			//	tempy = L2y[j];
			if (L2y[j] >= Ry[Rsizeminusone])  // && L2x[j] >= Rx[Rsizeminusone]
				;   // no operation, discard L2
			else if (L2x[j] > Rx[Rsizeminusone] && L2y[j] < Ry[Rsizeminusone]) {
				++Rsizeminusone;
				Rx[Rsizeminusone] = L2x[j];
				Ry[Rsizeminusone] = L2y[j];
			}
			//    else if (L2x[j] <= Rx[Rsizeminusone] && L2y[j] <= Ry[Rsizeminusone]) {
			else {
				Rx[Rsizeminusone] = L2x[j];
				Ry[Rsizeminusone] = L2y[j];
			}
			++j;
		}
	}
	++Rsizeminusone;
	Rx[Rsizeminusone] = Alenplusone;
	Ry[Rsizeminusone] = 0;

}

//BitPar-MergedLCS
int MLCS_Bitpar(vector<int> T, vector<int> A, vector<int> B, int Tlen, int Alen, int Blen, int size, int WordSize) {
	vector<vector<uint64_t> > W1(Blen + 1, vector<uint64_t>(ceil((Tlen + 1) / WordSize + 1)));
	vector<vector<uint64_t> > W2(Blen + 1, vector<uint64_t>(ceil((Tlen + 1) / WordSize + 1)));
	vector<vector<uint64_t> > Y(size, vector<uint64_t>(ceil((Tlen + 1) / WordSize + 1), 0));
	uint64_t carrybitW = 0;
	uint64_t MaxNUM64 = 0XFFFFFFFFFFFFFFFF;//2^64 - 1
	uint64_t modLastWord = 0X0000020000000000;//2^42
											  //uint64_t modLastWord = 0X0000000000000010;
	int WordNum = ceil((Tlen + 1) / WordSize + 1) - 1;

	//Constructing array of Y

	uint64_t shiftNum = 1;
	for (int i = 1; i <= min(Tlen, WordSize - 1); ++i) {
		shiftNum = shiftNum << 1;
		Y[T[i] - 1][0] |= shiftNum;
	}

	if (Tlen > 63)
		for (int i = WordSize; i <= Tlen; ++i) {
			shiftNum = shiftNum << 1;
			if (i % WordSize == 0)
				shiftNum = 1;
			Y[T[i] - 1][floor(i / WordSize)] |= shiftNum;
		}

	//Initialisation
	for (int i = 0; i <= WordNum; ++i) {
		W2[0][i] = MaxNUM64;


	}
	W2[0][WordNum] %= modLastWord;

	//Calculating boundaries
	uint64_t U = 0, tmpa = 0;
	for (int k = 1; k <= Blen; ++k) {
		carrybitW = 0;
		for (int i = 0; i <= WordNum - 1; ++i)
		{
			tmpa = W2[k - 1][i];
			U = tmpa & Y[B[k] - 1][i];
			W2[k][i] = (tmpa + U + carrybitW) | (W2[k - 1][i] - U);
			//detect overflow
			if (U > MaxNUM64 - tmpa)
				carrybitW = 1;
			else
				carrybitW = 0;

		}
		tmpa = W2[k - 1][WordNum];
		U = tmpa& Y[B[k] - 1][WordNum];
		W2[k][WordNum] = ((tmpa + U + carrybitW) % modLastWord) | (tmpa - U);
		//detect overflow
		if (U > MaxNUM64 - tmpa)
			carrybitW = 1;
		else
			carrybitW = 0;
	}

	int layer = 1;
	while (1) {
		if (layer % 2 == 1) {
			carrybitW = 0;
			for (int i = 0; i <= WordNum - 1; ++i)
			{
				tmpa = W2[0][i];
				U = tmpa & Y[A[layer] - 1][i];
				W1[0][i] = (tmpa + U + carrybitW) | (tmpa - U);
				//detect overflow
				if (U > MaxNUM64 - tmpa)
					carrybitW = 1;
				else
					carrybitW = 0;

			}
			tmpa = W2[0][WordNum];
			U = tmpa & Y[A[layer] - 1][WordNum];
			W1[0][WordNum] = ((tmpa + U + carrybitW) % modLastWord) | (tmpa - U);
			//detect overflow
			if (U > MaxNUM64 - tmpa)
				carrybitW = 1;
			else
				carrybitW = 0;

			//Main calculations
			uint64_t Up = 0, Upp = 0, Wp = 0, Wpp = 0, Ut = 0, Wt = 0, Vt = 0, carrybitWp = 0, carrybitWpp = 0;
			int tmpv = ceil(log2(WordSize)) - 1;
			for (int k = 1; k <= Blen; ++k) {
				carrybitWp = 0, carrybitWpp = 0;
				bool f = false;
				for (int i = 0; i <= WordNum - 1; ++i)
				{
					tmpa = W2[k][i];
					Up = tmpa & Y[A[layer] - 1][i];
					Wp = (tmpa + Up + carrybitWp) | (tmpa - Up);
					if (Up > MaxNUM64 - tmpa)
						carrybitWp = 1;
					else
						carrybitWp = 0;

					tmpa = W1[k - 1][i];
					Upp = tmpa & Y[B[k] - 1][i];
					Wpp = (tmpa + Upp + carrybitWpp) | (tmpa - Upp);
					if (Upp > MaxNUM64 - tmpa)
						carrybitWpp = 1;
					else
						carrybitWpp = 0;

					Ut = Wp | Wpp;
					Wt = Wp ^ Wpp;
					Vt = Wt;

					for (int ipp = 0; ipp <= tmpv; ++ipp) {
						Vt = Vt ^ (Vt << (1 << ipp));
					}
					if (f == true) Vt = ~Vt;
					if ((Vt & 0X8000000000000000) > 0) f = true;
					else f = false;
					W1[k][i] = ~(Wt & Vt) & Ut;
				}
				tmpa = W2[k][WordNum];
				Up = tmpa & Y[A[layer] - 1][WordNum];
				Wp = ((tmpa + Up + carrybitWp) % modLastWord) | (tmpa - Up);
				if (Up > MaxNUM64 - tmpa)
					carrybitWp = 1;
				else
					carrybitWp = 0;

				tmpa = W1[k - 1][WordNum];
				Upp = tmpa & Y[B[k] - 1][WordNum];
				Wpp = ((tmpa + Upp + carrybitWpp) % modLastWord) | (tmpa - Upp);
				if (Upp > MaxNUM64 - tmpa)
					carrybitWpp = 1;
				else
					carrybitWpp = 0;

				Ut = Wp | Wpp;
				Wt = Wp ^ Wpp;
				Vt = Wt;

				for (int ipp = 0; ipp <= tmpv; ++ipp) {
					Vt = Vt ^ (Vt << (1 << ipp));
				}
				if (f == true) Vt = ~Vt;
				if ((Vt & 0X8000000000000000) > 0) f = true;
				else f = false;
				W1[k][WordNum] = ~(Wt & Vt) & Ut;
			}

			if (layer == Alen)
			{
				int z = 0;
				uint64_t Vz = 0;
				for (int i = 0; i <= WordNum - 1; ++i) {
					Vz = ~W1[Blen][i];
					while (Vz != 0) {
						Vz = Vz & (Vz - 1);
						++z;
					}
				}
				Vz = (~W1[Blen][WordNum]) % modLastWord;
				while (Vz != 0) {
					Vz = Vz & (Vz - 1);
					++z;
				}
				return z;
			}
			++layer;
		}
		else {

			carrybitW = 0;
			for (int i = 0; i <= WordNum - 1; ++i)
			{
				tmpa = W1[0][i];
				U = tmpa & Y[A[layer] - 1][i];
				W2[0][i] = (tmpa + U + carrybitW) | (tmpa - U);
				//detect overflow
				if (U > MaxNUM64 - tmpa)
					carrybitW = 1;
				else
					carrybitW = 0;

			}
			tmpa = W1[0][WordNum];
			U = tmpa & Y[A[layer] - 1][WordNum];
			W2[0][WordNum] = ((tmpa + U + carrybitW) % modLastWord) | (tmpa - U);
			//detect overflow
			if (U > MaxNUM64 - tmpa)
				carrybitW = 1;
			else
				carrybitW = 0;

			//Main calculations
			uint64_t Up = 0, Upp = 0, Wp = 0, Wpp = 0, Ut = 0, Wt = 0, Vt = 0, carrybitWp = 0, carrybitWpp = 0;
			int tmpv = ceil(log2(WordSize)) - 1;
			for (int k = 1; k <= Blen; ++k) {
				carrybitWp = 0, carrybitWpp = 0;
				bool f = false;
				for (int i = 0; i <= WordNum - 1; ++i)
				{
					tmpa = W1[k][i];
					Up = tmpa & Y[A[layer] - 1][i];
					Wp = (tmpa + Up + carrybitWp) | (tmpa - Up);
					if (Up > MaxNUM64 - tmpa)
						carrybitWp = 1;
					else
						carrybitWp = 0;

					tmpa = W2[k - 1][i];
					Upp = tmpa & Y[B[k] - 1][i];
					Wpp = (tmpa + Upp + carrybitWpp) | (tmpa - Upp);
					if (Upp > MaxNUM64 - tmpa)
						carrybitWpp = 1;
					else
						carrybitWpp = 0;

					Ut = Wp | Wpp;
					Wt = Wp ^ Wpp;
					Vt = Wt;

					for (int ipp = 0; ipp <= tmpv; ++ipp) {
						Vt = Vt ^ (Vt << (1 << ipp));
					}
					if (f == true) Vt = ~Vt;
					if ((Vt & 0X8000000000000000) > 0) f = true;
					else f = false;
					W2[k][i] = ~(Wt & Vt) & Ut;
				}
				tmpa = W1[k][WordNum];
				Up = tmpa & Y[A[layer] - 1][WordNum];
				Wp = ((tmpa + Up + carrybitWp) % modLastWord) | (tmpa - Up);
				if (Up > MaxNUM64 - tmpa)
					carrybitWp = 1;
				else
					carrybitWp = 0;

				tmpa = W2[k - 1][WordNum];
				Upp = tmpa & Y[B[k] - 1][WordNum];
				Wpp = ((tmpa + Upp + carrybitWpp) % modLastWord) | (tmpa - Upp);
				if (Upp > MaxNUM64 - tmpa)
					carrybitWpp = 1;
				else
					carrybitWpp = 0;

				Ut = Wp | Wpp;
				Wt = Wp ^ Wpp;
				Vt = Wt;

				for (int ipp = 0; ipp <= tmpv; ++ipp) {
					Vt = Vt ^ (Vt << (1 << ipp));
				}
				if (f == true) Vt = ~Vt;
				if ((Vt & 0X8000000000000000) > 0) f = true;
				else f = false;
				W2[k][WordNum] = ~(Wt & Vt) & Ut;
			}
			if (layer == Alen)
			{
				int z = 0;
				uint64_t Vz = 0;
				for (int i = 0; i <= WordNum - 1; ++i) {
					Vz = ~W2[Blen][i];
					while (Vz != 0) {
						Vz = Vz & (Vz - 1);
						++z;
					}
				}
				Vz = (~W2[Blen][WordNum]) % modLastWord;
				while (Vz != 0) {
					Vz = Vz & (Vz - 1);
					++z;
				}
				return z;
			}
			++layer;
		}
	}
}

int MLCS_PengSparse(vector<int> T, vector<int> A, vector<int>B, int Tlen, int Alen, int Blen, int size)
{
	vector<vector<int> > nextA(size + 1, vector<int>(Alen + 1, 0));
	vector<vector<int> > nextB(size + 1, vector<int>(Blen + 1, 0));
	int len = 1;
	vector<vector<triple> > QL1(Tlen + 2, vector<triple>());
	vector<vector<triple> > QL2(Tlen + 2, vector<triple>());
	vector<vector<triple> > tempQL(Tlen + 2, vector<triple>());
	int qlsize = 0;
	int locationA = 0;
	int locationB = 0;
	//construct array of nextA
	for (int i = 1; i <= size; ++i) {
		locationA = Alen + 1;
		for (int j = Alen; j > 0; --j) {
			if (j == Alen)
				nextA[i][j] = locationA;
			if (A[j] == i) {
				locationA = j;
				nextA[i][j - 1] = locationA;
			}
			else {
				nextA[i][j - 1] = locationA;
			}
		}
	}

	//construct array of nextB
	for (int i = 1; i <= size; ++i) {
		locationB = Blen + 1;
		for (int j = Blen; j > 0; --j) {
			if (j == Blen)
				nextB[i][j] = locationB;
			if (B[j] == i) {
				locationB = j;
				nextB[i][j - 1] = locationB;
			}
			else {
				nextB[i][j - 1] = locationB;
			}
		}
	}


	int PosA = 0;
	int PosB = 0;
	int layer = 1;
	while (1) {
		if (layer % 2 == 1) {
			QL1[0].push_back(triple(0, 0));
			for (int l = 1; l <= len; l++)
			{
				for (std::vector<triple>::iterator mytriple = QL2[l - 1].begin(); mytriple != QL2[l - 1].end(); ++mytriple)
				{
					PosA = nextA[T[layer - 1]][mytriple->x];
					PosB = nextB[T[layer - 1]][mytriple->y];
					if (PosA <= Alen) {
						QL1[l].push_back(triple(PosA, mytriple->y));
					}
					if (PosB <= Blen) {
						QL1[l].push_back(triple(mytriple->x, PosB));
					}
				}
				if (QL2[l].size() > 0) {
					for (std::vector<triple>::iterator mytriple = QL2[l].begin(); mytriple != QL2[l].end(); ++mytriple)
					{
						QL1[l].push_back(triple(mytriple->x, mytriple->y));
					}
				}

				//2-D minima algorithm
				qlsize = QL1[l].size();
				if (qlsize > 1) {
					for (int i = 0; i < qlsize; ++i) {
						triple temptri = QL1[l][i];
						if (!tempQL[temptri.x].empty() && temptri.y < tempQL[temptri.x][0].y) {
							tempQL[temptri.x].erase(tempQL[temptri.x].begin());
							tempQL[temptri.x].push_back(temptri);
						}
						else if (tempQL[temptri.x].empty() == true) {
							tempQL[temptri.x].push_back(temptri);
						}
					}
					QL1[l].clear();
					int prevMaxY = Blen + 1;
					for (int i = 0; i <= Alen; ++i) {
						if (!tempQL[i].empty()) {
							if (tempQL[i][0].y < prevMaxY) {
								prevMaxY = tempQL[i][0].y;
								QL1[l].push_back(tempQL[i][0]);
							}
							tempQL[i].clear();
						}
					}
				}
			}//end for j
			if (QL1[len].size() > 0)
				len++;
			if (layer == Tlen + 1)
				break;
			++layer;
		}
		else {
			QL2[0].push_back(triple(0, 0));
			for (int l = 1; l <= len; l++)
			{
				for (std::vector<triple>::iterator mytriple = QL1[l - 1].begin(); mytriple != QL1[l - 1].end(); ++mytriple)
				{
					PosA = nextA[T[layer - 1]][mytriple->x];
					PosB = nextB[T[layer - 1]][mytriple->y];
					if (PosA <= Alen) {
						QL2[l].push_back(triple(PosA, mytriple->y));
					}
					if (PosB <= Blen) {
						QL2[l].push_back(triple(mytriple->x, PosB));
					}
				}
				if (QL1[l].size() > 0) {
					for (std::vector<triple>::iterator mytriple = QL1[l].begin(); mytriple != QL1[l].end(); ++mytriple)
					{
						QL2[l].push_back(triple(mytriple->x, mytriple->y));
					}
				}

				//2-D minima algorithm
				qlsize = QL2[l].size();
				if (qlsize > 1) {
					for (int i = 0; i < qlsize; ++i) {
						triple temptri = QL2[l][i];
						if (!tempQL[temptri.x].empty() && temptri.y < tempQL[temptri.x][0].y) {
							tempQL[temptri.x].erase(tempQL[temptri.x].begin());
							tempQL[temptri.x].push_back(temptri);
						}
						else if (tempQL[temptri.x].empty() == true) {
							tempQL[temptri.x].push_back(temptri);
						}
					}
					QL2[l].clear();
					int prevMaxY = Blen + 1;
					for (int i = 0; i <= Alen; ++i) {
						if (!tempQL[i].empty()) {
							if (tempQL[i][0].y < prevMaxY) {
								prevMaxY = tempQL[i][0].y;
								QL2[l].push_back(tempQL[i][0]);
							}
							tempQL[i].clear();
						}
					}
				}
			}//end for j
			if (QL2[len].size() > 0)
				len++;
			if (layer == Tlen + 1)
				break;
			++layer;
		}
	}

	return len - 1;
}

//sparse DP (build matchpair and boundedheap) [2014] Effective Sparse Dynamic Programming Algorithms for Merged and Block Merged LCS Problems
int MLCS_SPDP_Boundedheap(vector<int>T, vector<int>A, vector<int>B, int Tlen, int Alen, int Blen, int size)
{
	vector< pair<int, int> > mylistmatchAT, mylistmatchBT;
	vector< pair<int, int> > ::iterator it;
	vector<vector<int> >Boundedheap1(Alen + 1, vector<int>(Blen + 1));
	vector<vector<int> >Boundedheap2(Alen + 1, vector<int>(Blen + 1));

	int glMLCS = 0, temp;
	for (int i = 1; i <= Alen; i++) //AT match pair
	{
		for (int k = 1; k <= Tlen; k++)
		{
			if (A[i] == T[k])
				mylistmatchAT.push_back(pair<int, int>(i, k));
		}
	}
	for (int j = 1; j <= Blen; j++) //BT match pair
	{
		for (int k = 1; k <= Tlen; k++)
		{
			if (B[j] == T[k])
				mylistmatchBT.push_back(pair<int, int>(j, k));
		}
	}
	//hybrid
	if (mylistmatchAT.size() + mylistmatchBT.size() <= Tlen) {
		int layer = 1;
		while (1) {
			if (layer % 2 == 1) {

				for (int i = 0; i <= Alen; i++)
				{
					for (int j = 0; j <= Blen; j++) // H(k,i)=max(H(k-1,i),H(k,i-1))
					{
						if (i - 1 == -1)
							Boundedheap1[i][j] = Boundedheap2[i][j];
						else if (layer - 1 == 0)
							Boundedheap1[i][j] = Boundedheap1[i - 1][j];
						else
							Boundedheap1[i][j] = max(Boundedheap2[i][j], Boundedheap1[i - 1][j]);
					}
					for (it = mylistmatchBT.begin(); it != mylistmatchBT.end(); it++)//(*it).first = j (*it).second = k
					{
						if (layer == (*it).second)
						{
							for (int m = (*it).first - 1; m >= 0; m--)
							{
								Boundedheap1[i][(*it).first] = max(Boundedheap1[i][(*it).first], Boundedheap2[i][m] + 1);//boundedmax
								if (glMLCS < Boundedheap1[i][(*it).first])
									glMLCS = Boundedheap1[i][(*it).first];
							}
						}
					}

					if (binarysearch(mylistmatchAT, i, layer))
					{
						for (int j = 0; j <= Blen; j++)
						{
							for (int m = j; m >= 0; m--)
							{
								Boundedheap1[i][j] = max(Boundedheap1[i][j], Boundedheap2[i - 1][m] + 1);//boundedmax
								if (glMLCS < Boundedheap1[i][j])
									glMLCS = Boundedheap1[i][j];
							}
						}
					}
				}
				if (layer == Tlen)
					break;
				++layer;
			}
			else {
				for (int i = 0; i <= Alen; i++)
				{
					for (int j = 0; j <= Blen; j++) // H(k,i)=max(H(k-1,i),H(k,i-1))
					{
						if (i - 1 == -1)
							Boundedheap2[i][j] = Boundedheap1[i][j];
						else if (layer - 1 == 0)
							Boundedheap2[i][j] = Boundedheap2[i - 1][j];
						else
							Boundedheap2[i][j] = max(Boundedheap1[i][j], Boundedheap2[i - 1][j]);
					}
					for (it = mylistmatchBT.begin(); it != mylistmatchBT.end(); it++)//(*it).first = j (*it).second = k
					{
						if (layer == (*it).second)
						{
							for (int m = (*it).first - 1; m >= 0; m--)
							{
								Boundedheap2[i][(*it).first] = max(Boundedheap2[i][(*it).first], Boundedheap1[i][m] + 1);//boundedmax
								if (glMLCS < Boundedheap2[i][(*it).first])
									glMLCS = Boundedheap2[i][(*it).first];
							}
						}
					}

					if (binarysearch(mylistmatchAT, i, layer))
					{
						for (int j = 0; j <= Blen; j++)
						{
							for (int m = j; m >= 0; m--)
							{
								Boundedheap2[i][j] = max(Boundedheap2[i][j], Boundedheap1[i - 1][m] + 1);//boundedmax
								if (glMLCS < Boundedheap2[i][j])
									glMLCS = Boundedheap2[i][j];
							}
						}
					}
				}
				if (layer == Tlen)
					break;
				++layer;
			}

		}

		return glMLCS;
	}
	else
		return MLCS_PengSparse(T, A, B, Tlen, Alen, Blen, size);

}

int binarysearch(vector< pair<int, int> > P, int i, int k)
{
	int l = 0, r = P.size(), middle;

	while (l <= r)
	{
		middle = (l + r) / 2;
		if (P[middle].first == i)
		{
			if (P[middle].second == k)
				return 1;
			else if (P[middle].second > k)
				r = middle - 1;
			else
				l = middle + 1;
		}
		else if (P[middle].first > i)
			r = middle - 1;
		else
			l = middle + 1;
	}

	return 0;
}


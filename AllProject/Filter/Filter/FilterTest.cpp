// Filter.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "convolutions.h"
#include <string>
#include "Filter.h"
#include <stdlib.h>
#include <windows.h>
#include "AlgProcessTime.h"
using std::string;

int _tmain(int argc, _TCHAR* argv[])
{
	string strSrcFile = "D:\\Data\\can_tmr\\can_tmr.img";
	string strDstFile = "D:\\m3n3LaCFilter2.tif";

	/*float fCoefArray[9] = { 1.0/9,  1.0/9, 1.0/9,
							1.0/9, 1.0/9, 1.0/9,
							1.0/9, 1.0/9, 1.0/9, };*/
	/*float fCoefArray[25] = { 1, 1, 1, 1, 1,
							1, 1, 1, 1, 1,
							1, 1, 1, 1, 1, 
							1, 1, 1, 1, 1,
							1, 1, 1, 1, 1 };*/

	/*float fCoefArray[9] = { 1, 1, 1,
							 1, 1, 1, 
							 1, 1, 1 };*/
	float fCoefArray[9] = { 0, -1, 0,
							-1, 4, -1, 
							0, -1, 0 };
	float coef = 1;
	int m = 3;
	int n = 3;

	//CFilter::BoxLowPassMxNGeneric(strSrcFile, strDstFile, m, n, fCoefArray, coef);
	//LowPassMxNGeneric(strSrcFile, strDstFile, m, n, fCoefArray, coef);

	//LowPass(strSrcFile, strDstFile, 3, BORDER_REPLOCATE);
	//LowPassMxN(strSrcFile, strDstFile, 5, 5, fCoefArray);

	CAlgProcessTime::Alg_start();
	//CFilter::MedianLowPassMxNGeneric(strSrcFile, strDstFile, m, n);
	CFilter::LaplacianHighPassMxNGeneric(strSrcFile, strDstFile, m, n, fCoefArray, coef);

	CAlgProcessTime::Alg_end();

	printf("time = %lf s\n", CAlgProcessTime::GetAlgProcessTime());

	::system("pause");
	return 0;
}


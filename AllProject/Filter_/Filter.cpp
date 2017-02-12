// Filter.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "lowerpass.h"
#include <stdlib.h>

int _tmain(int argc, _TCHAR* argv[])
{
	
	const char *pSrcFile = "D:\\Data\\Filter\\qb_boulder_msi";
	const char *pDstFile = "D:\\Data\\Temp\\3.tif";
	int a = 1;

	TemplateOperation(pSrcFile, pDstFile, a);

	::system("pause");
	return 0;
}


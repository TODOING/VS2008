#include <stdlib.h>
#include "Eigen/Dense"
#include "StdRXD.h"

int main()
{
	const char *pszSrcFile = "D:\\Data\\rxd\\nvis_sub1_hsi.img";
	const char *pszRXDFile = "D:\\rxd.tif";

	CHighSpecRXD rxd(pszSrcFile);

	rxd.ExecuteStdRXD(pszRXDFile);


	::system("pause");
	return 0;
}
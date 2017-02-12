#include <stdlib.h>
#include "pca.h"

int main()
{
	//const char *src_file = "D:\\Data\\PCA\\qb_boulder_msi.img";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\qb_pca_my4_0_精确统计.tif";
	
	//const char *src_file = "D:\\Data\\PCA\\can_tmr.img";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr_my6_0.tif";
	//const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr6_0.sta";

	//const char *src_file = "D:\\Data\\PCA\\MCD43B3.A2016185.Albedo_BSA_Band_JL.tif";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\MCD43B3.tif";

	const char *src_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr_my6_0.tif";
	const char *inverse_pca_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr_my6_0_inverse_pca.tif";
	const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr6_0.sta";

	PCA pca(src_file);
	//pca.ExecutePCA(pca_file, statistic_file, 6);
	pca.ExecuteInversePCA(inverse_pca_file, statistic_file);

	system("pause");
	return 0;
}
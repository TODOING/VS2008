#include <stdlib.h>
#include "pca.h"
#include "AlgProcessTime.h"

int main()
{
	//const char *src_file = "D:\\Data\\PCA\\qb_boulder_msi.img";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\qb_pca_my4_0_精确统计.tif";
	
	//const char *src_file = "D:\\Data\\PCA\\can_tmr.img";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr_my6_12.tif";
	//const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr6_12.sta";

	const char *src_file = "D:\\Data\\PCA\\CASI_2012_09_08_145508_g.img";
	const char *pca_file = "D:\\Data\\PCA\\Temp\\my\\CASI_2012_09_08_145508_g_my24_0.tif";
	const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\CASI_2012_09_08_145508_g24_0.sta";

	//const char *src_file = "D:\\Data\\PCA\\MCD43B3.A2016185.Albedo_BSA_Band_JL.tif";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\MCD43B3.tif";

	/*const char *src_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr_my6_11.tif";
	const char *inverse_pca_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr_my6_21_inverse_pca.tif";
	const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr6_11.sta";*/

	PCA pca(src_file);
	CAlgProcessTime::Alg_start();

	pca.ExecutePCA(pca_file, statistic_file, 24, GDT_Float32);
	//pca.ExecuteInversePCA(inverse_pca_file, statistic_file, GDT_Int16);

	CAlgProcessTime::Alg_end();

	printf("time = %lf s\n", CAlgProcessTime::GetAlgProcessTime());

	system("pause");
	return 0;
}
#include <stdlib.h>
#include "pca.h"
#include "AlgProcessTime.h"
#include "mnf.h"

int main()
{
	//const char *src_file = "D:\\Data\\PCA\\CASI_2012_09_08_145508_g.img";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\my\\CASI_2012_09_08_145508_g_my6_int32_0.tif";
	//const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\CASI_2012_09_08_145508_g_my6_int32_0.sta";

	//const char *src_file = "D:\\Data\\PCA\\EO1H1320422004339110PZ_B.img";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\my\\EO1H1320422004339110PZ_B_my6_int32_0.tif";
	//const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\EO1H1320422004339110PZ_B_my6_int32_0.sta";

	//const char *src_file = "D:\\Data\\PCA\\CASI_2012_09_08_145508_g_BSQ.img";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\my\\CASI_2012_09_08_145508_g_my6_BSQ_2.tif";
	//const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\CASI_2012_09_08_145508_g_my6_BSQ_2.sta";

	//////////////////////////////////////////////////////////////////////////
	// 正变换
	const char *src_file = "D:\\Data\\PCA\\can_tmr.img";
	const char *pca_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr_my6_new_21.tif";
	const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr6_new_21.sta";

	//const char *src_file = "D:\\Data\\PCA\\EO1H1320422004339110PZ_B.img";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\my\\EO1H1320422004339110PZ_B_my163_Float32_2.tif";
	//const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\EO1H1320422004339110PZ_B_my163_Float32_2.sta";

	//const char *src_file = "D:\\Data\\PCA\\CASI_2012_09_08_145508_g_BSQ.img";
	//const char *pca_file = "D:\\Data\\PCA\\Temp\\my\\CASI_2012_09_08_145508_g_my6_BSQ_5.tif";
	//const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\CASI_2012_09_08_145508_g_my6_BSQ_5.sta";


	// 逆变换
	//const char *src_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr_my6_new_15.tif";
	//const char *inverse_pca_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr_my6_15_inverse_pca.tif";
	//const char *statistic_file = "D:\\Data\\PCA\\Temp\\my\\can_tmr6_new_15.sta";
	//////////////////////////////////////////////////////////////////////////

	/*PCA pca(src_file);
	

	pca.ExecutePCA(pca_file, statistic_file, 24, GDT_Float32);*/
	
	/*CAlgProcessTime::Alg_start();
	pca.ExecuteInversePCA(inverse_pca_file, statistic_file, GDT_Byte);
	CAlgProcessTime::Alg_end("执行InversePCA算法总耗时为");*/

	MNF mnf(src_file);
	mnf.ExecuteMNF(pca_file, statistic_file, 6, GDT_Float32);

	system("pause");
	return 0;
}
/*
 * This file is part of IIC-JKU QFR library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum/ for more information.
 */

#include "QuantumComputation.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/GoogleRandomCircuitSampling.hpp"
#include <ctime>

using namespace std;
using namespace chrono;

int main(int argc, char *argv[]) {

	clock_t start,end;
	double endtime;

	//读取real文件
	std::string fileName = "./circuits/test.real";
	if(argc==1)
    {
		std::cout<<"没有输入参数,使用默认电路"<<std::endl;
 	}
	else fileName = argv[1];

	// std::cout << std::endl << std::endl << "********************" << std::endl << "电路：" << fileName << std::endl;
	qc::QuantumComputation qc(fileName);
	//qc::QFT qc(16);

	//创建一个dd实例
	auto dd = make_unique<dd::Package>(); // create an instance of the DD package
	auto functionality = qc.buildFunctionality(dd);
	auto dd_size = dd->size(functionality);
	//qc.printMatrix(dd, functionality, std::cout); //输出矩阵
	//dd::export2Dot(functionality, "test-real.dot"); //输出QMDD表示文件
	// std::cout << "原始dd大小：" << dd->size(functionality) << "    "; //输出原始dd的大小

	qc::permutationMap varMap = qc.initialLayout; //复制电路的初始变量映射
	
	// qc.printPermutationMap(varMap); //打印初始映射

	// 调用sifting重排算法
	start=clock();
	auto reorderdd = dd->dynamicReorder(functionality, varMap, dd::DynamicReorderingStrategy::Sifting);
	// end=clock();
	//qc.printMatrix(dd, reorderdd, std::cout); //输出矩阵
	//dd::export2Dot(reorderdd, "reorderdd.dot"); //输出QMDD表示文件
	int sifting_size = dd->size(reorderdd);
	int one_sifting_size = sifting_size;
	// endtime=(double)(end-start)/CLOCKS_PER_SEC;
	auto one_sifting_time = endtime;

	//sifting 到收敛
/* 	int old_size=0;
	for(int i=0; i<3; ++i)
	{
		reorderdd = dd->dynamicReorder(reorderdd, varMap, dd::DynamicReorderingStrategy::Sifting);
		sifting_size = dd->size(reorderdd);
		if(sifting_size == old_size) {
			break;
		} else {
			old_size = sifting_size;
		}
	} */

	end=clock();
	endtime=(double)(end-start)/CLOCKS_PER_SEC;
	auto conv_sifting_time = endtime;
	int conv_sifting_size = sifting_size;
	cout << fileName << ' ' << dd_size << ' ' << one_sifting_size << ' ' << one_sifting_time << ' ' << conv_sifting_size << ' ' << conv_sifting_time << endl;
	return 0;
}


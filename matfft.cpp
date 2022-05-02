/*
 *zhangchang2317@mails.jlu.edu.cn
 *
 *
 *2022/4/29 
* 重力托普利兹矩阵fft加速
 *张敞 zhangchang2317@mails.jlu.edu.cn
*/
#define EIGEN_FFTW_DEFAULT
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
#include ".\Eigen\Dense"
#include <iostream>
#include ".\unsupported\Eigen\FFT"
#include ".\fftw3.h"

using namespace matlab::data;
using namespace std;
using namespace Eigen;
using matlab::mex::ArgumentList;
class MexFunction : public matlab::mex::Function {
// // // // // // // // // // // // // // // // // // // // // //     
    	//! Extracts the pointer to underlying data from the non-const iterator (`TypedIterator<T>`).
	/*! This function does not throw any exceptions. */
	template <typename T>
	inline T* toPointer(const matlab::data::TypedIterator<T>& it) MW_NOEXCEPT {
		static_assert(std::is_arithmetic<T>::value && !std::is_const<T>::value,
			"Template argument T must be a std::is_arithmetic and non-const type.");
		return it.operator->();
	}
	template <typename T>
	inline T* getPointer(matlab::data::TypedArray<T>& arr) MW_NOEXCEPT {
		static_assert(std::is_arithmetic<T>::value, "Template argument T must be a std::is_arithmetic type.");
		return toPointer(arr.begin());
	}
	template <typename T>
	inline const T* getPointer(const matlab::data::TypedArray<T>& arr) MW_NOEXCEPT {
		return getPointer(const_cast<matlab::data::TypedArray<T>&>(arr));
	}
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //     
public:

    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        ArrayFactory factory;
        Eigen::FFT<double> fft;
        int nt = inputs[0][0];
        int nr = inputs[0][1];
        int ns = inputs[0][2];
        TypedArray<double> Gz1_input = std::move(inputs[1]);
        auto Gz1_ptr = getPointer(Gz1_input);
        MatrixXd Gz1 = Map<MatrixXd>(Gz1_ptr,ns*nt*nr*3,ns*nt);
        
       
        TypedArray<double> AAA_input = std::move(inputs[2]);
        auto AAA_ptr = getPointer(AAA_input);
        MatrixXd AAA = Map<MatrixXd>(AAA_ptr,nt,ns);
        
        VectorXd p_zero(1,ns);
        p_zero = VectorXd::Zero(ns);
        MatrixXd Z(nt*nr*3,ns);
        Z = MatrixXd::Zero(nt*nr*3,ns);
        VectorXd D1(ns*nt);
        VectorXd d1(ns);
        VectorXd d2(ns);
        VectorXd gg_ext1(2*ns);
        VectorXd t_ext1(2*ns);
        VectorXd p_ext1(2*ns);
        VectorXcd tmpOut1(2*ns);
        VectorXcd tmpOut2(2*ns);
        VectorXcd tmpOut3(2*ns);
for( int zz = 0; zz < nt*nr*3; zz++){
        
        
    for(int ii= 0; ii<nt; ii++){
        D1=Gz1.row((zz)*ns);
        d1=D1.segment(ii*ns,ns);

        d2=AAA.row(ii);
        t_ext1.head(ns) = d1;
        t_ext1(ns) = 0;
        t_ext1.tail(ns-1) = d1.tail(ns-1).reverse();
        p_ext1.head(ns) = d2;  
        p_ext1.tail(ns) = p_zero;
       
        fft.fwd(tmpOut1, t_ext1);
        fft.fwd(tmpOut2, p_ext1);
        tmpOut3 = tmpOut1.cwiseProduct(tmpOut2);
 
        fft.inv(gg_ext1, tmpOut3);
//             cout<<gg_ext1.transpose().rows()<<gg_ext1.transpose().cols()<<endl;
        Z.row(zz) = Z.row(zz) + gg_ext1.head(ns).transpose();

        
    }

}
       std::vector<size_t> size_output(1,nt*nr*3);
//        size_output.insert(size_output.begin(),nt);
       size_output.insert(size_output.end(),ns);       
       double* Zptr = Z.data();
       outputs[0] = factory.createArray<double>(size_output,Zptr,Zptr + ns*nt*nr*3);
    }

    void arrayProduct(matlab::data::TypedArray<double>& matrix, double multiplier) {
//     cout <<"matrix:" <<  matrix << endl;
//     cout <<"&matrix:" << & matrix << endl;
    }


};
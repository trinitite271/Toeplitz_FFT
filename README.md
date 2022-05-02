# 通过托普利兹矩阵加速重力反演
目标是使用Eigen代替matlab循环起到加速作用，实际测试表明matlabFFT比FFTW快很多。  
期待Eigen早日支持MKL 后端。
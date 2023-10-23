# vmd
Variational Mode Decomposition (VMD) in C++

# Dependencies
- Eigen
- fftw


# Example
```cpp
int T = 1000;
float fs = 1. / T;
Eigen::VectorXf t = Eigen::VectorXf::LinSpaced(T, 0, 1);
Eigen::VectorXf input(T);

float f1 = 2.0f, f2 = 24.0f, f3 = 288.0f;
for (int i = 0; i < T; i++) {
  input(i) = cos(2 * M_PI * f1 * t(i)) 
       + cos(2 * M_PI * f2 * t(i)) 
       + cos(2 * M_PI * f3 * t(i));
}
input += 0.1 * Eigen::VectorXf::Random(T);

int IMFs = 3;                                   // number of modes                    
vmd result(input, IMFs);          

cout << result.omega;
Eigen::MatrixXf & u = result.u;	                // the collection of decomposed modes
Eigen::MatrixXf & u_hat = result.u_hat;	        // spectra of the modes
Eigen::VectorXcf & omega = result.omega;        // estimated mode center-frequencies
```
- see `vmd.h` for more info.

# note
- It works most of the time, but not every time. Please keep in mind that it may not give you the same result as the vmd in Matlab do.
- Please kindly let me know if you find anything wrong.
- This implementation would be overhauled later to obtain better accuracy and efficiency.


# ref
1. Function for decomposing a signal according to the Variational Mode Decomposition ([Dragomiretskiy and Zosso, 2014](https://doi.org/10.1109/TSP.2013.2288675)) method.
2. DodgeHo's work, see https://github.com/DodgeHo/VMD_cpp

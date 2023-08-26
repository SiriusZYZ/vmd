# vmd
Variational Mode Decomposition (VMD) in C++

# Dependencies
- Eigen
- fftw


# Example
```cpp
// full test will be updated later. you can try it yourself first
int IMFs = 5;
vmd result(signal, IMFs);
Eigen::MatrixXf & u = result.u;				//	the collection of decomposed modes
Eigen::MatrixXf & u_hat = result.u_hat;		//	spectra of the modes
Eigen::VectorXcf & omega = result.omega;	//	estimated mode center-frequencies
```
- see `vmd.h` for more info.

# note
- It works most of the time, but not every time. Please keep in mind that it may not give you the same result as the vmd in Matlab do.
- Please kindly let me know if you find anything wrong.

# ref
Function for decomposing a signal according to the Variational Mode Decomposition ([Dragomiretskiy and Zosso, 2014](https://doi.org/10.1109/TSP.2013.2288675)) method.
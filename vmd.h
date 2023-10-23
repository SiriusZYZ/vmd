#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#ifndef EIGEN_FFTW_DEFAULT
#define EIGEN_FFTW_DEFAULT
#endif // !EIGEN_FFTW_DEFAULT

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif


/// @brief Perform variational mode decomposition on the given signal.
/// @example
///	vmd result(signal, IMFs);
///	Eigen::MatrixXf &u = result.u;	//	the collection of decomposed modes
///	Eigen::MatrixXf &u_hat = result.u_hat;	//	spectra of the modes
///	Eigen::VectorXcf &omega = result.omega;	//	estimated mode center-frequencies
class vmd
{


public:
	///	the collection of decomposed modes
	Eigen::MatrixXf u;

	///	spectra of the modes
	Eigen::MatrixXcf u_hat;

	///	estimated mode center-frequencies
	Eigen::VectorXf omega;

private:

	int IterationNum = 0;
	float Udiff;

public:


	/// @brief Perform variational mode decomposition on a given signal.
	/// @param signal : The time domain signal(1D) to be decomposed.
	/// @param tol : tolerance of convergence criterion; (default as 1e-7)
	/// @param K : The number of modes to be recovered.
	/// @param alpha : The balancing parameter of the data - fidelity constraint (default as 1000).
	/// @param tau : Time - step of the dual ascent(pick 0 for noise - slack) (default as 0).
	/// @param DC : True if the first mode is put and kept at DC(0 - freq) (default as 0).
	/// @param init :	0 as all mode center-frequencies will start at 0; 
	///					1 as all mode center-frequencies start uniformly distributed; 
	///					2 as all mode center-frequencies initialized randomly; (currently unsupported)
	///					(default as 0)
	vmd(const Eigen::VectorXf& signal, int K, float tol = 1e-7, float alpha = 1000, float tau = 0, int DC = 0, int init = 0);

	/// @brief return how many iterations has been done.
	/// @return iteration number.
	int return_iteration_number();

	/// @brief return Udiff
	/// @return Udiff
	float retrun_diff();


private:

	/// @brief Perform fft on the given signal.
	/// @param x : The input 1D vector to be calculate.
	/// @return The fft of x, with the same length as x.
	inline Eigen::VectorXcf fft(const Eigen::VectorXf& x);

	/// @brief Shift zero-frequency component to center of spectrum
	/// @param X : The input 1D vector to be shifted.
	/// @return The shifted signal.
	inline Eigen::VectorXcf fftshift(const Eigen::VectorXcf& X);

	/// @brief Inverse fast Fourier transform
	/// @param X : The input 1D vector to be calculate.
	/// @return The ifft of X, with the same length as X.
	inline Eigen::VectorXcf ifft(const Eigen::VectorXcf& X);

	/// @brief Inverse zero-frequency shift
	/// @param X : The input 1D vector to be calculate.
	/// @return The shifted signal.
	inline Eigen::VectorXcf ifftshift(const Eigen::VectorXcf& X);


};

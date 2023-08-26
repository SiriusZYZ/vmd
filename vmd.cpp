#include "vmd.h"

Eigen::VectorXcf vmd::fft(const Eigen::VectorXf& x)
{
	Eigen::FFT<float> fft;

	Eigen::VectorXcf res;

	fft.fwd(res, x, x.size());

	return res;


}

Eigen::VectorXcf vmd::fftshift(const Eigen::VectorXcf& X)
{
	int n = X.size();
	int k = (n + 1) / 2;
	Eigen::VectorXcf w(n);

	w(Eigen::seq(0, n - k - 1)) = X(Eigen::seq(k, n - 1));
	w(Eigen::seq(n - k, n - 1)) = X(Eigen::seq(0, k - 1));

	return w;

}

Eigen::VectorXcf vmd::ifft(const Eigen::VectorXcf& X)
{

	Eigen::FFT<float> fft;
	Eigen::VectorXcf res;
	fft.inv(res, X, X.size());
	return res;

}

Eigen::VectorXcf vmd::ifftshift(const Eigen::VectorXcf& X)
{
	int n = X.size();
	int k = (n + 1) / 2;
	Eigen::VectorXcf w(n);
	w(Eigen::seq(0, k - 1)) = X(Eigen::seq(n - k, n - 1));
	w(Eigen::seq(k, n - 1)) = X(Eigen::seq(0, n - k - 1));
	return w;
}


vmd::vmd(const Eigen::VectorXf& signal, int K, float tol, float alpha, float tau, int DC, int init)
{

	/*--------Preparations--------*/

	// Period and sampling frequency of input signal
	int save_T = signal.size();

	bool even = (save_T % 2 == 0) ? true : false;
	Eigen::VectorXf E_signal;

	// replicate the last element if the length of signal is odd.
	if (even)
	{
		E_signal = signal;
	}
	else
	{
		E_signal = Eigen::VectorXf(save_T + 1);
		E_signal.segment(0, save_T) = signal;
		E_signal(save_T) = signal(save_T - 1);
	}



	save_T = E_signal.size();

	float fs = 1. / save_T;

	//  extend the signal by mirroring
	int T = save_T;
	Eigen::VectorXf f_mirror(2 * T);
	f_mirror.segment(0, T / 2) = E_signal.segment(0, T / 2).reverse();
	f_mirror.segment(T / 2, T) = E_signal.segment(0, T);
	f_mirror.segment(3. / 2 * T, T / 2) = E_signal.segment(T / 2, T / 2).reverse();

	//  Time Domain 0 to T (of mirrored signal)
	T = f_mirror.rows();
	Eigen::VectorXf t = Eigen::VectorXf::LinSpaced(T, 1, T).array() / T;

	//  Spectral Domain discretization
	Eigen::VectorXf freqs = t.array() - (0.5 + 1. / T);

	//  Maximum number of iterations.
	int N = 500;

	Eigen::VectorXf Alpha = Eigen::VectorXf::Ones(K).array() * alpha;

	//  Construct and center f_hat
	auto f_hat = fftshift(fft(f_mirror));
	auto f_hat_plus = f_hat;
	f_hat_plus.segment(0, T / 2).setZero();

	//  matrix keeping track of every iterant // could be discarded for mem

	std::vector<Eigen::MatrixXcf> u_hat_plus = std::vector<Eigen::MatrixXcf>(2, Eigen::MatrixXcf::Zero(freqs.rows(), K));
	Eigen::MatrixXcf omega_plus = Eigen::MatrixXcf::Zero(N, K);

	switch (init)
	{
	case(1):
		omega_plus.row(0) = Eigen::VectorXcf::LinSpaced(K, 0, K - 1).array() * (0.5 / K);
		break;
	case(2):
		// not done yet. here, if init == 2, Initialize the central frequencies as random numbers distributed uniformly in the interval [0,0.5] cycles/sample. 
	default:
		break;
	}

	if (DC)
	{
		omega_plus(0, 0) = std::complex<float>(0., 0.);
	}

	// Eigen::MatrixXcf lambda_hat = Eigen::MatrixXcf::Zero(N, freqs.rows());
	Eigen::VectorXcf lambda_hat = Eigen::VectorXcf::Zero(freqs.rows());

	float eps = 1.1921e-07;

	Udiff = tol + eps;

	int n = 0, k;

	Eigen::VectorXcf sum_uk = Eigen::VectorXcf::Zero(T);

	while (Udiff > tol && n < N - 1)
	{
		//   update first mode accumulator
		k = 0;
		sum_uk = u_hat_plus[0].col(K - 1) + sum_uk - u_hat_plus[0].col(0);

		//  update spectrum of first mode through Wiener filter of residuals

		Eigen::MatrixXcf Dividend_vec = f_hat_plus - sum_uk - (lambda_hat / 2.0);
		Eigen::MatrixXcf Divisor_vec = (1 + alpha * ((freqs.array() - omega_plus(n, k))).array().square());
		u_hat_plus[1].col(k) = Dividend_vec.cwiseQuotient(Divisor_vec);

		if (!DC)
		{
			Eigen::VectorXcf left = freqs(Eigen::seq(T / 2, T - 1));
			Eigen::VectorXcf right = u_hat_plus[1].col(k)(Eigen::seq(T / 2, T - 1)).array().abs().square();
			omega_plus(n + 1, k) = left.transpose() * right;
			omega_plus(n + 1, k) /= u_hat_plus[1].col(k)(Eigen::seq(T / 2, T - 1)).array().abs().square().sum();

		}

		//  update of any other mode
		for (int k = 1; k < K; k++)
		{
			// accumulator
			sum_uk = u_hat_plus[1].col(k - 1) + sum_uk - u_hat_plus[0].col(k);

			// mode spectrum
			Eigen::MatrixXcf Dividend_vec = f_hat_plus - sum_uk - (lambda_hat / 2.0);
			Eigen::MatrixXcf Divisor_vec = (1 + alpha * ((freqs.array() - omega_plus(n, k))).array().square());
			u_hat_plus[1].col(k) = Dividend_vec.cwiseQuotient(Divisor_vec);

			//  center frequencies
			Eigen::VectorXcf left = freqs(Eigen::seq(T / 2, T - 1));
			Eigen::VectorXcf right = u_hat_plus[1].col(k)(Eigen::seq(T / 2, T - 1)).array().abs().square();
			omega_plus(n + 1, k) = left.transpose() * right;
			omega_plus(n + 1, k) /= u_hat_plus[1].col(k)(Eigen::seq(T / 2, T - 1)).array().abs().square().sum();

			/*std::complex<float> Dividend{ 0,0 }, Divisor{ 0, 0 }, Addend{ 0, 0 };
			for (int i = 0; i < T - T / 2; i++) {
				Addend = abs(u_hat_plus[n](T / 2 + i, k)) * abs(u_hat_plus[n](T / 2 + i, k));
				Divisor += Addend;
				Dividend += freqs[T / 2 + i] * Addend;
			}
			omega_plus(n, k) = Dividend / Divisor;*/

		}

		// Dual ascent

		lambda_hat = lambda_hat + tau * (u_hat_plus[1].rowwise().sum() - f_hat_plus);

		// loop counter 
		n += 1;

		Udiff = eps;
		std::complex<float> acc{ eps, 0 };

		for (int i = 0; i < K; i++)
		{
			auto tmp = (u_hat_plus[1].col(i) - u_hat_plus[0].col(i));
			acc += 1. / T * tmp.adjoint() * tmp;

		}

		u_hat_plus[0] = u_hat_plus[1];

		Udiff = abs(acc);

		//std::cout << n << "  " << Udiff << std::endl;
	}


	/*-------postprocessing and cleanup--------*/
	N = (N > n) ? n : N;
	IterationNum = N;
	omega = omega_plus.row(N).real();

	u_hat = Eigen::MatrixXcf::Zero(T, K);

	//u_hat.block(T / 2, 0, T / 2, K) = u_hat_plus[N].block(T / 2, 0, T / 2, K);

	for (int i = 0; i < T / 2; i++)
	{
		u_hat.row(T / 2 + i) = u_hat_plus[1].row(T / 2 + i);
		u_hat.row(T / 2 - i) = u_hat_plus[1].row(T / 2 + i).adjoint().transpose();
	}


	u_hat.row(0) = u_hat.row(T - 1).transpose().adjoint();

	Eigen::MatrixXf u_tmp = Eigen::MatrixXf::Zero(K, T);
	for (int k = 0; k < K; k++)
	{
		u_tmp.row(k) = ifft(ifftshift(u_hat.col(k))).real();
	}
	u = Eigen::MatrixXf(K, signal.size());
	u = u_tmp.block(0, T / 4, K, signal.size());

	for (int k = 0; k < K; k++)
	{
		u_hat.col(k) = fftshift(fft(u_tmp.row(k))).transpose();
	}


}

int vmd::return_iteration_number()
{
	return IterationNum;
}

float vmd::retrun_diff()
{
	return Udiff;
}
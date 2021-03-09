 #ifndef PLANCK_XCOMPLEX_H
 #define PLANCK_XCOMPLEX_H
 
 #include <complex>
 
 #define xcomplex std::complex
 typedef xcomplex<float>  fcomplex;
 typedef xcomplex<double> dcomplex;
 
 template<typename T> inline xcomplex<T> times_i(xcomplex<T> inp)
   { return xcomplex<T>(-inp.imag(),inp.real()); }
 
 #endif

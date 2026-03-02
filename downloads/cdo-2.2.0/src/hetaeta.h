#ifndef _HETAETA_H
#define _HETAETA_H

template <typename T>
void hetaeta(bool ltq, int ngp, const int *imiss, int nlev1, const double *ah1, const double *bh1, const Varray<double> &fis1,
             const Varray<double> &ps1, const Varray<T> &t1, const Varray<T> &q1, int nlev2, const double *ah2, const double *bh2,
             const Varray<double> &fis2, Varray<double> &ps2, Varray<T> &t2, Varray<T> &q2, int nvars, const Varray2D<T> &vars1,
             Varray2D<T> &vars2, Varray<double> &scor, Varray<double> &pscor, Varray<double> &secor);

#endif /* _HETAETA_H */

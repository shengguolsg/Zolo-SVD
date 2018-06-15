# Zolo-SVD

This package is written by Shengguo Li, National University of Defense Technology, Changsha, China

It includes QDWH-PD, QDWH-SVD and Zolo-PD, Zolo-SVD algorithms.
It uses ELPA package to compute the eigendecompositions when computing SVD, and
we used "elpa-2016.05.003".

A structured QR factorization algorithm is also included, which are modified from Scalapack routines,
and they are mpdgeqrf.f, mpdorgqr.f, respectively.
They are usually faster than Scalapack routines.

If you have any questions, you can send email to
nudtlsg@gmail.com
shengguolsg@126.com

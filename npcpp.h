
extern"C" {
  void setreac(int L, int P, int Z, int A);

  void setsigt(int T);
  void setsigu(int U);
  void setsiga(int A);
  void kinemin(double * Emin);
  void kinen(double E, double * numin, double * numax);
  void kinenq(double E, double nu, double * Q2min, double * Q2max);
  void kineqn(double E, double Q2, double * numin, double * numax);
  void kinenqt(double E, double nu, double Q2, double * tmin, double * tmax);
  void kinex(double E, double * xmin, double * xmax);
  void kinexy(double E, double x, double * ymin, double * ymax);
  void kiney(double E, double * ymin, double * ymax);
  void kineyx(double E, double y, double * xmin, double * xmax);
  void kinexyt(double E, double x, double y, double * tmin, double tmax);
  void kinexynq(double E, double x, double y, double * nu, double * Q2);
  void kinenqxy(double E, double nu, double Q2, double * x, double * y);
  double signa(double E);
  double signas(double sN);
  double signam(double mu);
  double signatt();
  double signat(double t);
  double signanu(double E, double nu);
  double signanq(double Q2);
  double signaq2(double E, double Q2);
  double siganqn(double nu);
  double signax(double E, double x);
  double signaxx(double y);
  double signay(double E, double y);
  double signayy(double x);
  double signaxy(double E, double x, double y);
  //  double signaxyt(double E, double x, double y, double t);
  double signanqt(double E, double nu, double Q2, double t);
  void phicoh(double * PhiC);
  double phicohb(double b);
  void phiinc(double * PhiI);
}

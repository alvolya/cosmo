// Benchmark helper (not part of CoSMo): builds H exactly as DEXHVJTcsb does,
// then reports nonzeros, build time and the time of one stored-H matvec.
// Usage: csb_matvec_bench system.xml [n_matvec=20]
#include <cstdint>
typedef uint16_t spsint;
typedef unsigned long long uint_nbasis;
#define INFO yes
#include <av/argparser.cxx>
#include <av/statusclass.cxx>
#include <omp.h>
#include <xsm/DEXHVcsb.cxx>
int main(int argc, char** argv) {
  if (argc < 2) { std::cerr << "usage: csb_matvec_bench sys.xml [nmv]\n"; return 1; }
  int nmv = argc > 2 ? atoi(argv[2]) : 20;
  XMLNode xMainNode = XMLNode::openFileHelper(argv[1], xsm::cosmo);
  XMLNode xSystemNode = xml::SelectNode(xMainNode, xsm::system, "name");
  SM::SingleParticleStates sps; Make(sps, xMainNode.getChildNode(xsm::valencespace));
  SM::Many_Body_States st; SM::Get(st, xMainNode, xSystemNode);
  tav::STensor<spsint,double> ESP(2), VSPX(4);
  xReadSPInt(ESP, VSPX, xSystemNode, xMainNode, sps);
  SM::PPMatrix HH(st); HH.push_back(&ESP); HH.push_back(&VSPX);
  double t0 = omp_get_wtime();
  uint16_t br = std::sqrt(1.*st.n) + 1;
  tav::SparseMatrixCSB<double,uint32_t,uint16_t> FHH(st.n, br, br, true);
  SM::PPMatrix2CSBBlk(FHH, HH);
  double tb = omp_get_wtime() - t0;
  unsigned long long nnz = 0;
  for (unsigned i = 0; i < FHH.numBlockRow; ++i) for (unsigned j = 0; j < FHH.numBlockCol; ++j) nnz += FHH.Blocks[i][j].nnz;
  std::vector<double> x(st.n), y(st.n);
  for (uint_nbasis i = 0; i < st.n; ++i) x[i] = std::sin(0.37 * i + 1.0);
  FHH.TimesVector(y.data(), x.data());  // warm-up
  double t1 = omp_get_wtime();
  for (int k = 0; k < nmv; ++k) FHH.TimesVector(y.data(), x.data());
  double tm = (omp_get_wtime() - t1) / nmv;
  // stored = upper triangle incl. diagonal; full-matrix nonzeros = 2*stored - D (approx; diagonal stored once)
  double nnz_full = 2.0 * nnz - double(st.n);
  printf("BENCH D=%llu threads=%d nnz_stored=%llu nnz_full~%.4g nnz_per_row~%.1f build_s=%.3f matvec_s=%.5f ns_per_stored_nz=%.3f ns_per_full_nz=%.3f build_us_per_full_nz=%.4f\n",
         (unsigned long long)st.n, omp_get_max_threads(), nnz, nnz_full, nnz_full / st.n, tb, tm,
         1e9 * tm / nnz, 1e9 * tm / nnz_full, 1e6 * tb / nnz_full);
  return 0;
}

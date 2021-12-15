#include <boost/timer/progress_display.hpp>
#include <libint2.hpp>
#include <stdexcept>
#include <vector>

using std::cout;
using std::endl;
using std::ifstream;
using std::invalid_argument;
using std::ofstream;
using std::stod;
using std::string;
using std::vector;

int main(int argc, char *argv[]) {
  // Setup
  if (argc < 6) {
    throw invalid_argument("Arguments: xyzfile, basis name, gamma, "
                           "threshold, output");
  }
  string xyzfile = argv[1];
  string bname = argv[2];
  double gamma_val = stod(argv[3]);
  double threshold = stod(argv[4]);
  string out = argv[5];
  libint2::initialize();

  cout << "Parameters:" << endl;
  cout << "  "
       << "xyzfile: " << xyzfile << endl;
  cout << "  "
       << "basis: " << bname << endl;
  cout << "  "
       << "gamma: " << gamma_val << endl;
  cout << "  "
       << "threshold: " << threshold << endl;
  cout << "  "
       << "output: " << out << endl;

  // Geometry and Basis
  ifstream ifile(xyzfile);
  auto geom = libint2::read_dotxyz(ifile);
  libint2::BasisSet basis(bname, geom);

  // Compute the number of electrons
  int nel = 0;
  for (int i = 0; i < geom.size(); ++i) {
    nel += geom[i].atomic_number;
  }

  // Create the Engine
  libint2::Engine eng(libint2::Operator::yukawa, basis.max_nprim(),
                      basis.max_l());
  eng.set_params(gamma_val);

  // Results stored here
  vector<int> rows, cols;
  vector<double> vals;

  // Compute
  const auto &buffer = eng.results();
  auto shell_lookup = basis.shell2bf();
  boost::timer::progress_display prog(basis.size());
  for (int i = 0; i < basis.size(); ++i) {
    for (int j = 0; j < basis.size(); ++j) {
      if (i > j)
        continue;

      eng.compute(basis[i], basis[i], basis[j], basis[j]);
      auto data = buffer[0];

      // check screening
      if (data == nullptr)
        continue;

      // unpack
      int offset_i = shell_lookup[i];
      int offset_j = shell_lookup[j];
      int n = basis[i].size();
      int m = basis[j].size();
      for (int k = 0; k < n; ++k) {
        for (int l = 0; l < m; ++l) {
          if (data[k * n + l] > threshold) {
            vals.push_back(data[k * n + l]);
            rows.push_back(offset_i + k);
            cols.push_back(offset_j + l);
          }
        }
      }
    }
    ++prog;
  }

  // Write out statistics
  cout << "System Information:" << endl;
  cout << "  Sparsity: " << double(vals.size()) / (basis.nbf() * basis.nbf()) << endl;
  cout << "  NBF: " << basis.nbf() << endl;
  cout << "  Electrons: " << nel << endl;

  // Write out to file
  ofstream ofile(out);
  ofile << "%%MatrixMarket matrix coordinate real symmetric" << endl;
  ofile << "%" << endl;
  ofile << basis.nbf() << " " << basis.nbf() << " " << vals.size() << endl;
  for (int i = 0; i < vals.size(); ++i) {
    ofile << rows[i] + 1 << " " << cols[i] + 1 << " " << vals[i] << endl;
  }
  ofile.close();

  // Cleanup
  libint2::finalize();
  return 0;
}

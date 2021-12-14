#include <libint2.hpp>
#include <stdexcept>
#include <vector>

int main(int argc, char* argv[]) {
  // Setup
  if (argc < 6) {
    throw std::invalid_argument("Arguments: xyzfile, basis name");
  }
  std::string xyzfile = argv[1];
  std::string bname = argv[2];
  double gamma_val = std::stod(argv[3]);
  double threshold = std::stod(argv[4]);
  std::string out = argv[5];
  libint2::initialize();

  // Geometry and Basis
  std::ifstream ifile(xyzfile);
  auto geom = libint2::read_dotxyz(ifile);
  libint2::BasisSet basis(bname, geom);

  // Create the Engine
  libint2::Engine eng(libint2::Operator::yukawa, basis.max_nprim(), basis.max_l());
  eng.set_params(gamma_val);

  // Results stored here
  std::vector<int> rows, cols;
  std::vector<double> vals;

  // Compute
  const auto& buffer = eng.results();
  auto shell_lookup = basis.shell2bf();
  for (int i = 0; i < basis.size(); ++i) {
    for (int j = 0; j < basis.size(); ++j) {
      if (i > j) continue;

      eng.compute(basis[i], basis[i], basis[j], basis[j]);
      auto data = buffer[0];

      // check screening
      if (data == nullptr) continue;

      // unpack
      int offset_i = shell_lookup[i];
      int offset_j = shell_lookup[j];
      int n = basis[i].size();
      int m = basis[j].size();
      for (int k = 0; k < n; ++k) {
	for (int l = 0; l < m; ++l) {
	  if (data[k*n + l] > 1e-8) {
	    vals.push_back(data[k*n + l]);
	    rows.push_back(offset_i + k);
	    cols.push_back(offset_j + l);
	  }
	}
      }
    }
  }

  // Write out to file
  std::ofstream ofile(out);
  ofile << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
  ofile << "%" << std::endl;
  ofile << basis.nbf() << " " << basis.nbf() << " " << vals.size() << std::endl;
  for (int i = 0; i < vals.size(); ++i) {
    ofile << rows[i] + 1 << " " << cols[i] + 1 << " " << vals[i] << std::endl;
  }
  ofile.close();

  // Cleanup
  libint2::finalize();
  return 0;
}


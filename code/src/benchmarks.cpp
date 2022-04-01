#include <iostream>
#include <sstream>
#include <chrono>
#include <array>

namespace benchmarks {

  template <typename Test, int avg_over = 10, int n_runs = 10>
  struct Benchmark {
    void generate() {
      double start, end;
      for (int i = 0; i < avg_over; i++) {
        Test test;
        for (int j = 0; j < n_runs; j++) {
          start = time();
          test.run();
          end = time();
          this->times[i][j] = end - start;
        };
      };
    };

    std::string get_results() const {
      std::stringstream results;
      std::array<double, n_runs> sums;
      for (int i = 0; i < avg_over; i++) {
        for (int j = 0; j < n_runs; j++) {
          sums[j] += this->times[i][j];
        };
      };
      results << "showing results for " << n_runs
              << " operations averaged across " << avg_over << " executions." 
              << std::endl;
      for (int i = 0; i < n_runs; i++) {
        results << "cycle " << i+1 << " ran in " << sums[i] / avg_over 
                << " ms." << std::endl;
      };
      return results.str();
    };

    private:
    std::array<std::array<double, n_runs>, avg_over> times;
    double time() {
      return std::chrono::time_point_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    };
  };

};

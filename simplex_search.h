#ifndef SIMPLEX_SEARCH_H
#define SIMPLEX_SEARCH_H

#include <vector>
#include <functional>

class SimplexSearch {
    public:
        std::vector<double> Search(std::function<double(std::vector<double>)> error_function, std::vector<double> initial_guess);
        int MaxNumberIterations = 200;
        double MinFunctionConvergence = 1e-10;
};

#endif // SIMPLEX_SEARCH_H

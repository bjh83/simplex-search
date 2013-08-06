#include "simplex_search.h"

#include <vector>
#include <functional>
#include <cmath>
#include <gtest/gtest.h>

namespace {
    using namespace std;

    double Function(vector<double> coordinates) {
        const double& x = coordinates.front();
        const double& y = coordinates.back();
        return cos(x) * cos(y);
    }

    TEST(SimplexSearchTest, IdentifiesMinimum) {
        const double pi = 3.141592654;
        SimplexSearch seeker;
        vector<double> initial = {1, 0};
        vector<double> minimum = seeker.Search(Function, initial);
        EXPECT_EQ(1, minimum.front());
        EXPECT_EQ(pi, minimum.back());
    }

} // namespace

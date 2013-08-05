#include "simplex_search.h"

using namespace std;

namespace {

    class Vertex {
        public:
            Vertex(function<double(vector<double>)> error_function, vector<double> input)
                : error_function_(error_function), function_input_(input) {}

            Vertex(int size) : function_input_(size, 0);

            double error_value() const {
                if(is_error_value_set_) {
                    return error_value_;
                } else {
                    error_value_ = error_function_(function_input_);
                    is_error_value_set_ = true;
                    return error_value_;
                }
            }

            const vector<double>& function_input() const { return function_input_; }

            int size() const { return function_input_.size(); }

            bool operator<(const Vertex& vertex) const {
                return error_value() < vertex.error_value();
            }

            Vertex& operator+=(const Vertex& vertex) {
                return apply_operator<+=>(vertex);
            }

            Vertex& operator+=(double value) {
                return apply_operator<+=>(value);
            }

            Vertex operator-=(const Vertex& vertex) const {
                return apply_operator<-=>(vertex);
            }

            Vertex operator-=(double value) const {
                return apply_operator<-=>(vertex);
            }

            Vertex operator*=(const Vertex& vertex) const {
                return apply_operator<*=>(vertex);
            }

            Vertex operator*=(double value) const {
                return apply_operator<*=>(value);
            }

            Vertex operator/=(const Vertex& vertex) {
                return apply_operator</=>(vertex);
            }

            Vertex operator/=(double value) const {
                return apply_operator</=>(vertex);
            }

        private:
            mutable bool is_error_value_set_ = false;
            mutable double error_value_;
            function<double(vector<double>)> error_function_;
            vector<double> function_input_;

            template<OP>
            Vertex& apply_operator(const Vertex& vertex) {
                for(int i = 0; i < size(); i++) {
                    function_input_[i] OP vertex.function_input()[i];
                }

                is_error_value_set_ = false;
                return *this;
            }

            template<OP>
            Vertex& apply_operator(double value) {
                for(int i = 0; i < size(); i++) {
                    function_input_[i] OP vertex.function_input()[i];
                }

                is_error_value_set_ = false;
                return *this;
            }
    };

    Vertex operator+(Vertex left, const Vertex& right) {
        left += right;
        return left;
    }

    Vertex operator+(Vertex left, double right) {
        left += right;
        return left;
    }

    Vertex operator+(double left, Vertex right) {
        right += left;
        return right;
    }

    Vertex operator-(Vertex left, const Vertex& right) {
        left -= right;
        return left;
    }

    Vertex operator-(Vertex left, double right) {
        left -= right;
        return left;
    }

    Vertex operator*(Vertex left, const Vertex& right) {
        left *= right;
        return left;
    }

    Vertex operator*(Vertex left, double right) {
        left *= right;
        return left;
    }

    Vertex operator*(double left, Vertex right) {
        right *= left;
        return right;
    }

    Vertex operator/(Vertex left, const Vertex& right) {
        left /= right;
        return left;
    }

    Vertex operator/(Vertex left, double right) {
        left /= right;
        return left;
    }

    vector<Vertex> InitializeSimplex(function<double(vector<double>)> error_function, vector<double> input) {
        const int dimensions = input.size();
        vector<Vertex> simplex(dimensions + 1);
        simplex[0] = Vertex(error_function, input);

        for(int i = 0; i < dimensions; i++) {
            vector<double> next_value = input;
            next_value[i] *= i;
            simplex[i + 1] = Vertex(error_function, next_value);
        }

        sort(simplex.begin(), simplex.end());
        return simplex;
    }

    Vertex ComputeCentroid(const vector<Vertex>& best_values) {
        Vertex centroid(best_values[0].size());
        for(const Vertex& vertex: best_values) {
            centroid += vertex;
        }
        centroid /= centroid.size();
        return centroid;
    }

    Vertex Reflect(const Vertex& centroid, const Vertex& worst_point) {
        return 2 * centroid - worst_point;
    }

    Vertex Expand(const Vertex& centroid, const Vertex& reflected_point) {
        return 2 * reflected_point - centroid;
    }

    Vertex Contract(const Vertex& centroid, const Vertex& worst_point) {
        Vertex distance = (centroid - worst_point) / 2;
        Vertex point1 = centroid + distance;
        Vertex point2 = centroid - distance;
        return point1 < point2 ? point1 : point2;
    }

    void Shrink(vector<Vertex>& simplex) {
        const Vertex& best_point = simplex[0];
        for(int i = 1; i < simplex.size(); i++) {
            simplex[i] -= (simplex[i] - best_point) / 2;
        }
    }

} // namespace

vector<double> SimplexSearch::Search(function<double(vector<double>)> error_function, vector<double> input) {
    vector<Vertex> simplex = InitializeSimplex(error_function, input);
    double error = simplex.front().error_value();
    const double kErrorThreshold = 10;
    while(error > kErrorThreshold) {
        const Vertex& best_point = simplex.front();
        Vertex& worst_point = simplex.back();
        vector<Vertex> best_points(simplex.begin(), simplex.end() - 1);
        Vertex centroid = ComputeCentroid(best_points);
        Vertex reflected_point = Reflect(centroid, worst_point);

        if(new_point < worst_point) {
            if(best_point < reflected_point) {
                worst_point = reflected_point;
            } else {
                Vertex expanded_point = Expand(centroid, reflected_point);
                if(expanded_point < best_point) {
                    worst_point = expanded_point;
                } else {
                    worst_point = reflected_point;
                }
            }
        } else {
            Vertex contracted_point = Contract(centroid, worst_point);
            if(contracted_point < worst_point) {
                worst_point = contracted_point;
            } else {
                Shrink(simplex);
            }
        }

        sort(simplex.begin(), simplex.end());
    }

    return simplex.front().function_input();
}

#include "simplex_search.h"

#include <algorithm>
#include <iostream>

using namespace std;

namespace {

    class Vertex {
        public:
            Vertex(function<double(vector<double>)> error_function, vector<double> input)
                : error_function_(error_function), function_input_(input) {}

            Vertex(function<double(vector<double>)> error_function, int size) : error_function_(error_function), function_input_(size, 0) {}

            Vertex(const Vertex& vertex) : is_error_value_set_(vertex.is_error_value_set_), 
                error_value_(vertex.error_value_),
                error_function_(vertex.error_function_),
                function_input_(vertex.function_input_) {}

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

            const function<double(vector<double>)>& error_function() const { return error_function_; }

            int size() const { return function_input_.size(); }

            Vertex& operator=(Vertex vertex) {
                swap(vertex, *this);
                return *this;
            }

            bool operator<(const Vertex& vertex) const {
                return error_value() < vertex.error_value();
            }

            Vertex& operator+=(const Vertex& vertex) {
                return apply_operator(plus<double>(), vertex);
            }

            Vertex& operator+=(double value) {
                return apply_operator(plus<double>(), value);
            }

            Vertex& operator-=(const Vertex& vertex) {
                return apply_operator(minus<double>(), vertex);
            }

            Vertex& operator-=(double value) {
                return apply_operator(minus<double>(), value);
            }

            Vertex& operator*=(const Vertex& vertex) {
                return apply_operator(multiplies<double>(), vertex);
            }

            Vertex& operator*=(double value) {
                return apply_operator(multiplies<double>(), value);
            }

            Vertex& operator/=(const Vertex& vertex) {
                return apply_operator(divides<double>(), vertex);
            }

            Vertex& operator/=(double value) {
                return apply_operator(divides<double>(), value);
            }

            friend ostream& operator<<(ostream& stream, const Vertex& vertex) {
                return stream << "error value: " << vertex.error_value() << " x: " << vertex.function_input().front() << " y: " << vertex.function_input().back() << "\n";
            }

            friend void swap(Vertex& left, Vertex& right) {
                swap(left.is_error_value_set_, right.is_error_value_set_);
                swap(left.error_value_, right.error_value_);
                swap(left.error_function_, right.error_function_);
                swap(left.function_input_, right.function_input_);
            }

        private:
            mutable bool is_error_value_set_ = false;
            mutable double error_value_;
            function<double(vector<double>)> error_function_;
            vector<double> function_input_;

            inline Vertex& apply_operator(function<double(double, double)> op, const Vertex& vertex) {
                for(int i = 0; i < size(); i++) {
                    function_input_[i] = op(function_input_[i], vertex.function_input()[i]);
                }

                is_error_value_set_ = false;
                return *this;
            }

            inline Vertex& apply_operator(function<double(double, double)> op, double value) {
                for(int i = 0; i < size(); i++) {
                    function_input_[i] = op(function_input_[i], value);
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
        vector<Vertex> simplex(dimensions + 1, Vertex(error_function, input));

        for(int i = 0; i < dimensions; i++) {
            vector<double> next_value = input;
            next_value[i] += i;
            simplex[i + 1] = Vertex(error_function, next_value);
        }

        sort(simplex.begin(), simplex.end());
        return simplex;
    }

    Vertex ComputeCentroid(const vector<Vertex>& best_values) {
        Vertex centroid(best_values.front().error_function(), best_values.front().size());
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
        for(unsigned int i = 1; i < simplex.size(); i++) {
            simplex[i] -= (simplex[i] - best_point) / 2;
        }
    }

} // namespace

vector<double> SimplexSearch::Search(function<double(vector<double>)> error_function, vector<double> input) {
    vector<Vertex> simplex = InitializeSimplex(error_function, input);
    double new_error = simplex.front().error_value();
    double old_error;
    int i = 0;
    while(i < MaxNumberIterations) {
        i++;
        const Vertex& best_point = simplex.front();
        Vertex& worst_point = simplex.back();
        vector<Vertex> best_points(simplex.begin(), simplex.end() - 1);
        Vertex centroid = ComputeCentroid(best_points);
        Vertex reflected_point = Reflect(centroid, worst_point);

        if (worst_point.error_value() - best_point.error_value() < MinFunctionConvergence) {
            break;
        }

        if(reflected_point < worst_point) {
            if(best_point < reflected_point) {
                worst_point = reflected_point;
                cout << "Reflect\n";
            } else {
                Vertex expanded_point = Expand(centroid, reflected_point);
                if(expanded_point < best_point) {
                    worst_point = expanded_point;
                    cout << "Expand\n";
                } else {
                    worst_point = reflected_point;
                    cout << "Reflect\n";
                }
            }
        } else {
            Vertex contracted_point = Contract(centroid, worst_point);
            if(contracted_point < worst_point) {
                worst_point = contracted_point;
                cout << "Contract\n";
            } else {
                Shrink(simplex);
                cout << "Shrink\n";
                cout << "Best Point: " << best_point
                     << "Worst Point: " << worst_point
                     << "Reflected Point: " << reflected_point
                     << "Contracted Point: " << contracted_point
                     << "Centroid: " << centroid;
            }
        }

        sort(simplex.begin(), simplex.end());
        old_error = new_error;
        new_error = simplex.front().error_value();
    }

    return simplex.front().function_input();
}

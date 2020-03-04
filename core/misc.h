#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

int count_zeros(std::vector<bool> f);
int count_ones(std::vector<bool> f);
void print_formula(std::vector<bool> f);
void print_formula_int(std::vector<bool> f);
void print_err(std::string m);
bool ends_with(std::string const & value, std::string const & ending);
bool starts_with(std::string const & value, std::string const & start);
std::string exec(const char* cmd);
int random_number();
void shuffle_ints(std::vector<int> &v);

std::string convert_vec(std::vector<bool> node);

std::vector<int> subtract_int(std::vector<bool> &a, std::vector<bool> &b);
std::vector<bool> subtract(std::vector<bool> &a, std::vector<bool> &b);
std::vector<bool> intersection(std::vector<bool> &a, std::vector<bool> &b);
bool is_subset(std::vector<int> &a, std::vector<bool> &b);
bool is_subset(std::vector<bool> &a, std::vector<bool> &b);
bool is_disjoint(std::vector<bool> &a, std::vector<int> &b);
std::vector<bool> union_sets(std::vector<bool> &a, std::vector<bool> &b);
std::vector<bool> union_sets(std::vector<bool> &a, std::vector<int> &b);
std::vector<bool> complement(std::vector<bool> &a);
void trim(std::string &f);
bool is_hitting_pair(std::vector<int> cl1, std::vector<int> cl2);

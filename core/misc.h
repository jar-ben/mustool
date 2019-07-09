#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

int count_zeros(std::vector<bool> &f);
int count_ones(std::vector<bool> &f);
void print_formula(std::vector<bool> f);
void print_formula_int(std::vector<bool> f);
void print_err(std::string m);
bool ends_with(std::string const & value, std::string const & ending);
std::string exec(const char* cmd);
int random_number();
std::string convert_vec(std::vector<bool> node);

bool is_subset(std::vector<int> &a, std::vector<bool> &b);
bool is_subset(std::vector<bool> &a, std::vector<bool> &b);
bool is_disjoint(std::vector<bool> &a, std::vector<int> &b);
std::vector<bool> union_sets(std::vector<bool> &a, std::vector<bool> &b);
std::vector<bool> union_sets(std::vector<bool> &a, std::vector<int> &b);
void trim(std::string &f);
bool is_hitting_pair(std::vector<int> cl1, std::vector<int> cl2);



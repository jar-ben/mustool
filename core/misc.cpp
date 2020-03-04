#include "misc.h"
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <random>
#include <regex>

using namespace std;

void trim(string &f){
	f.erase( remove(f.begin(), f.end(), '\n'), f.end() );
	regex pattern("  ");
	f = regex_replace(f, pattern, "");	
}

int count_zeros(vector<bool> f){
	return std::count(f.begin(), f.end(), false);
}

int count_ones(vector<bool> f){
	return std::count(f.begin(), f.end(), true);
}

void print_formula(vector<bool> f){
	for(int i = 0; i < f.size(); i++)
		cout << (f[i])? "1" : "0";
	cout << endl;
}

void print_formula_int(vector<bool> f){
	for(int i = 0; i < f.size(); i++)
	       if(f[i]) cout << i << " ";
	cout << endl;	
}

void print_err(string m){
	cout << endl;
	cout << "########## ########## ########## ########## ########## " << endl;
	cout << m << endl;
	cout << "########## ########## ########## ########## ########## " << endl;
	cout << endl << endl << endl;
	exit(1);
}

bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

bool starts_with(std::string const & value, std::string const & start)
{
    if (start.size() > value.size()) return false;
    return value.rfind(start, 0) == 0;
}

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}

int random_number(){
	std::random_device rd; 
	std::mt19937 eng(rd());
	std::uniform_int_distribution<> distr(1, 1000000000); 

	return distr(eng);
}

void shuffle_ints(std::vector<int> &shuffled){
	auto rng = std::default_random_engine {};
	std::shuffle(std::begin(shuffled), std::end(shuffled), rng);
}

std::string convert_vec(std::vector<bool> node){
	string result = "";
	for(auto b: node)
		result += (b)? "1" : "0";
	return result;
}

bool is_hitting_pair(std::vector<int> cl1, std::vector<int> cl2){
	for(auto l1: cl1){
		for(auto l2: cl2){
			if(l1 == (-1 * l2)) return true;
		}
	}
	return false;
}

bool is_subset(std::vector<int> &a, std::vector<bool> &b){
	for(auto &constraint: a)
		if(!b[constraint])
			return false;
	return true;
}

bool is_subset(std::vector<bool> &a, std::vector<bool> &b){
	int dimension = a.size();
	for(int i = 0; i < dimension; i++)
		if(a[i] && !b[i])
			return false;
	return true;
}

bool is_disjoint(std::vector<bool> &a, std::vector<int> &b){
	for(auto c: b)
		if(a[c])
			return false;
	return true;
}

std::vector<bool> intersection(std::vector<bool> &a, std::vector<bool> &b){
	std::vector<bool> inter = a;
	for(int i = 0; i < a.size(); i++)
		if(!b[i]) inter[i] = false;
	return inter;
}

//subtracts a from b
std::vector<int> subtract_int(std::vector<bool> &a, std::vector<bool> &b){
	std::vector<int> inter;
	for(int i = 0; i < a.size(); i++)
		if(b[i] && !a[i]) inter.push_back(i);
	return inter;
}

//subtracts a from b
std::vector<bool> subtract(std::vector<bool> &a, std::vector<bool> &b){
	std::vector<bool> inter = b;
	for(int i = 0; i < a.size(); i++)
		if(a[i]) inter[i] = false;
	return inter;
}

std::vector<bool> complement(std::vector<bool> &a){
	std::vector<bool> b(a.size(), false);
	for(int i = 0; i < a.size(); i++)
		if(!a[i]) b[i] = true;
	return b;
}

std::vector<bool> union_sets(std::vector<bool> &a, std::vector<bool> &b){
	int dimension = a.size();
	std::vector<bool> res(dimension, false);
	for(int i = 0; i < dimension; i++)
		if(a[i] || b[i])	res[i] = true;
	return res;
}

std::vector<bool> union_sets(std::vector<bool> &a, std::vector<int> &b){
	int dimension = a.size();
	std::vector<bool> res = a;
	for(auto &l: b)
		res[l] = true;
	return res;
}



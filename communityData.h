//
//  communityData.h
//  transitionMatrix
//
//  Created by Will Pearse on 18/11/2011.
//  Copyright 2011 Imperial College London. All rights reserved.
//

#ifndef transitionMatrix_communityData_h
#define transitionMatrix_communityData_h


#include <functional>
#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <numeric>
#include <boost/bind.hpp>
#include "communitySimulator.h"
#include <math.h>
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>

class Community{
private:
    std::vector<int> years;
    std::vector<std::string> species_names, sp_names_repro_death;
    std::string community_name;
    int n_species,row,column;
    void set_e_m(void);
    void set_t_m(double stationary_prob=0.6);
    boost::numeric::ublas::matrix<double> t_m;
    std::vector<boost::numeric::ublas::matrix<int> > e_m;
    void initialise(void);
    std::vector<boost::numeric::ublas::matrix<int> > real_e_m;
    
public:
    std::vector< std::vector <std::string> > communities;
    Community(std::string species, std::string abundance, std::string year, std::string name);
    Community(boost::numeric::ublas::matrix<double> transition_matrix, std::vector<std::string> sp_names, int n_communities, int starting_length, std::string name, int seed=0);
    
    void add_species(std::string species, std::string abundance, std::string year);
    std::string name(void);
    std::vector<int>::size_type n_years(void);
    std::vector<std::string> extract_year(int year);
    std::vector<std::string> extract_index(int index);
    double calc_likelihoods(void);
    double optimise(int max_iter=100);
    boost::numeric::ublas::matrix<double> print_transition_matrix(void);
    boost::numeric::ublas::matrix<int> print_event_matrix(int index);
    boost::numeric::ublas::matrix<int> print_real_event_matrix(int index);
    std::vector<std::string> print_community(int index, int width=10);
    
    friend class DataSet;
    friend double likelihood_parameter_com(double prob, int row, int column, std::vector<boost::numeric::ublas::matrix<int> > e_m, boost::numeric::ublas::matrix<double> t_m, int minimising);
    friend double likelihood_parameter_data(double prob, int row, int column, std::vector<Community> communities, boost::numeric::ublas::matrix<double> t_m, int minimising);
};

class DataSet{
private:
    void set_t_m(double stationary_prob=0.6);
public:
    DataSet(const char *file, double prob=0.6);
    double internal_optimise(double prob);
    void optimise(int max_iter=100);
    void partition(int max_iter=100);
    boost::numeric::ublas::matrix<double> print_transition_matrix(int width=4);
    std::vector<Community> communities;
    
    int n_communities,row,column;
    std::vector<int> n_years;
    std::vector<std::string> community_names, species_names;
    boost::numeric::ublas::matrix<double> t_m;

};

#endif

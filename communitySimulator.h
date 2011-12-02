//
//  communitySimulator.h
//  transitionMatrix
//
//  Created by Will Pearse on 31/10/2011.
//  Copyright 2011 Imperial College London. All rights reserved.
//

#ifndef transitionMatrix_communitySimulator_h
#define transitionMatrix_communitySimulator_h

//HEADERS
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/tools/minima.hpp>
#include <vector>
#include <string>


//PROTOTYPES
std::vector<std::string> next_step(boost::numeric::ublas::matrix<double> transition_matrix, std::vector<std::string> current_community, std::vector<std::string> species_names);

std::vector<double> likelihood(boost::numeric::ublas::matrix<double> transition_matrix, std::vector<std::string> first_community, std::vector<std::string> second_community, std::vector<std::string> species_names);

std::vector<double> likelihood(boost::numeric::ublas::matrix<double> transition_matrix, boost::numeric::ublas::matrix<int> event_matrix);

int find_max(std::vector<double> vec);

int best_transition_to_sp(boost::numeric::ublas::matrix<double> t_m, int sp);

std::vector<std::vector<double> > matrix_to_vector(boost::numeric::ublas::matrix<double> transition_matrix);

std::vector<int> extract_transitions(std::string first_state, std::string second_state, std::vector<std::string> first_community, std::vector<std::string> second_community);

boost::numeric::ublas::matrix<int> likely_transitions(boost::numeric::ublas::matrix<double> transition_matrix, const std::vector<std::string>&first_community, const std::vector<std::string>& second_community, std::vector<std::string> species_names);

boost::numeric::ublas::matrix<double> optimise(boost::numeric::ublas::matrix<int> event_matrix, boost::numeric::ublas::matrix<double> guess_matrix);

std::vector<boost::numeric::ublas::matrix<int> > make_event_matrix(std::vector<std::vector<std::string> > communities, boost::numeric::ublas::matrix<double> transition_matrix, std::vector<std::string> species_names);

class BrentWrapper{
    int c_l,n_e;
public:
    BrentWrapper(unsigned community_length)
    {
        c_l = community_length;
    }
    void set_events(unsigned no_events)
    {
        n_e = no_events;
    }
    
    double operator()(double prob)
    {
        double bin_coeff;
        
        bin_coeff = (boost::math::factorial<double>(c_l) / (boost::math::factorial<double>(n_e) * boost::math::factorial<double>(c_l - n_e)));
        
        return (1 - (bin_coeff * pow(prob, n_e) * pow((1 - prob), (c_l - n_e))));
        
    }

};
#endif
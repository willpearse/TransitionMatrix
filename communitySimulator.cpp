//
//  communitySimulator.cpp
//  transitionMatrix
//
//  Created by Will Pearse on 04/11/2011.
//  Copyright 2011 Imperial College London. All rights reserved.
//

#include "communitySimulator.h"
#include "communityData.h"

//Namespaces
using namespace std;
using namespace boost::numeric;

vector<double> likelihood(ublas::matrix<double> transition_matrix, ublas::matrix<int> event_matrix)
{
    //Setup
    //int n_events=0,x=0;
    /*for(i=0; i<event_matrix.size1(); ++i)
        for(k=0; k<event_matrix.size2(); ++k)
            n_events += event_matrix(i,k);*/
    vector<double> likelihoods;//(n_events);
    
    //Loop through events and calculate likelihood
    // - check to see if double looping to pre-allocate is worth it
    for(int i=0; i<event_matrix.size1(); ++i)
        for(int k=0; k<event_matrix.size2(); ++k)
        {
            int to_do=event_matrix(i,k);
            while(to_do)
            {
                likelihoods.push_back(transition_matrix(i,k));
                --to_do;
            }
        }
    //Return and check
    //assert(x == n_events);
    return likelihoods;
}
vector<double> log_likelihood(ublas::matrix<double> transition_matrix, ublas::matrix<int> event_matrix)
{
    vector<double> likelihoods;
    for(int i=0; i<event_matrix.size1(); ++i)
        for(int k=0; k<event_matrix.size2(); ++k)
            while(event_matrix(i,k)--)
                likelihoods.push_back(log(transition_matrix(i,k)));
    return likelihoods;
}
vector<double> likelihood_null(vector<double> transition_null, vector<string> community, vector<string> species_names)
{
    vector<double> likelihoods(community.size());
    
    //Add each species' log likelihood to the likelihoods list
    for(int i=0; i<community.size(); ++i)
        for(int j=0; j<species_names.size(); ++j)
            if(species_names[j]==community[i])
            {
                likelihoods[i] = transition_null[j];
                break;
            }
    return likelihoods;
}
//Needed because you can't pass iterators around like indices
// - very dirty hack
int find_max(vector<double> vec)
{
    double curr_max = -1.0;
    int max_pos = -1;
    for(int i=0; i<vec.size(); i++)
        if(vec[i] > curr_max)
        {
            max_pos = i;
            curr_max = vec[i];
        }
    assert(max_pos != -1 && curr_max != -1); //We have found a new maximum
    return max_pos;
}

//What's the best way of getting to a certain species in a given transition matrix?
int best_transition_to_sp(ublas::matrix<double> t_m, int sp)
{
    double curr_max = -1.0;
    int max_pos = -1,i=0;
    for(; i<t_m.size1(); ++i)
    {
        if(t_m(i, sp) > curr_max)
        {
            max_pos = i;
            curr_max = t_m(i, sp);
        }
    }
    //We only want to consider reproduction as a last resort
    // - the value of 0.0 comes from other code that sets the value at that if we've run out of those species to use in the transition
    if(curr_max == 0.0)
    {
        //If we're not going to reproduction (a nonsense) or dying, check repro.
        if(sp < t_m.size1())
            if(t_m(sp, ++i) > curr_max)
            {
                max_pos = i;
                curr_max = t_m(sp, i);
            }
    }
    assert(max_pos != -1 && curr_max != -1); //We have found a new maximum
    return max_pos;
}

vector<vector<double> > matrix_to_vector(ublas::matrix<double> transition_matrix)
{
    //Make vector of vectors of transition_matrix, and store the best 
    // - this will be unnecessary when you remove boost
    int NSPECIES = transition_matrix.size1();
    vector<vector<double> > species_probs(NSPECIES);
    vector<int> best_prob(NSPECIES);
    for(int i=0; i<NSPECIES; i++)
    {
        species_probs[i].resize(transition_matrix.size2());
        for(int j=0; j<transition_matrix.size2(); j++)
            species_probs[i][j] = transition_matrix(i,j);
    }
    return species_probs;
}

ublas::matrix<int> likely_transitions(ublas::matrix<double> transition_matrix, const vector<string>& first_community, const vector<string>& second_community, vector<string> species_names)
{
    //Number of species in the list
    assert(species_names.size() == transition_matrix.size2());
    assert(species_names.size() == (transition_matrix.size1() +2));
    const unsigned REPRODUCTION = (species_names.size()-1);
    const int NSPECIES = (species_names.size() - 2);
    
    //Make an events matrix (to be returned)
    ublas::matrix<int> events_matrix(NSPECIES, NSPECIES+2);
    events_matrix = ublas::zero_matrix<int>(NSPECIES, NSPECIES+2);
    
    //Make some looping indices and an assertion checker
    int i,j;
    
    //Make a new, padded, second community
    vector<string> pad_second_community;
    if (second_community.size() > first_community.size())
        pad_second_community.resize(second_community.size());
    else
        pad_second_community.resize(first_community.size());
    
    //Asign second community
    if (second_community.size() > first_community.size())
        pad_second_community = second_community;
    else
    {
        for(i=0; i<first_community.size(); ++i)
            if(i < second_community.size())
                pad_second_community[i] = second_community[i];
            else
                pad_second_community[i] = "DEATH";
    }
    
    //Count of remaining species and no. species just staying the same
    vector<int> remaining_sp(NSPECIES), stationary_sp(NSPECIES);
    for(i=0; i<first_community.size(); i++)
        for(j=0; j<=NSPECIES; j++)
            if(species_names[j] == first_community[i])
                ++remaining_sp[j];
    //To handle death and reproduction (hack?)
    remaining_sp.push_back(pad_second_community.size());
    remaining_sp.push_back(pad_second_community.size());
    fill(stationary_sp.begin(), stationary_sp.end(), 0);
    
    //What're the most likely ways of getting to a species?
    vector<int> most_likely(NSPECIES+2);
    for(i=0; i<(NSPECIES+2); ++i)
        most_likely[i] = best_transition_to_sp(transition_matrix, i);
    
    //Loop through second padded community
    for(i=0; i<pad_second_community.size(); ++i)
    {
        //What species are we dealing with?
        for(j=0; j<NSPECIES; ++j)
            if(pad_second_community[i] == species_names[j])
                break;
        
        //Do something, but not infinitely
        // - should re-write the order of the nesting...
        int max_iter = 0,finished=0;
        while(max_iter < NSPECIES+2 && finished == 0)
        {
            //Are we reproducing?
            if(most_likely[j] == REPRODUCTION)
            {
                //Do we have enough free slots?
                if(remaining_sp[j] > 0)
                {
                    ++events_matrix(j,REPRODUCTION);
                    --remaining_sp[most_likely[j]];
                    finished = 1;
                }
                else
                {
                    //No, so can we use a previous stable transition?
                    if(stationary_sp[j] > 0)
                    {
                        ++events_matrix(j,REPRODUCTION);
                        //--events_matrix(j,j); Not under the new method - reproduction and stationarity are separate
                        --stationary_sp[j]; 
                        finished = 1;
                    }
                }
            }
            else
            {
                //Not reproducing; can we do a transition?
                if(remaining_sp[most_likely[j]] > 0)
                {
                    ++events_matrix(most_likely[j],j);
                    --remaining_sp[most_likely[j]];
                    //Stable transition?
                    if(most_likely[j] == j)
                        ++stationary_sp[j];
                    finished = 1;
                }
            }
            //Find a new minimum if we've had no luck this time
            if(!finished)
            {
                if(most_likely[j] == REPRODUCTION)
                    transition_matrix(j,REPRODUCTION) = 0.0;
                else
                    transition_matrix(most_likely[j],j) = 0.0;
                most_likely[j] = best_transition_to_sp(transition_matrix, j);
            }
            ++max_iter;
        }
        
        //Has there been a massive influx of species?
        // - these should be dealt with better than this!
        if(!finished)
        {
            //Do a reproduction
            ++events_matrix(j,REPRODUCTION);
            --remaining_sp[most_likely[j]];
            finished = 1;
        }
        
        //We have done something
        //assert(finished);
    }
        
    return events_matrix;
}

vector<int> extract_transitions(string first_state, string second_state, vector<string> first_community, vector<string> second_community)
{
    //Make output and assert
    vector<int> target(first_community.size(), 0);
    assert(first_community.size() == second_community.size());//We've already reconstructed the transition history
    
    //Initialise the loop
    vector<string>::size_type i = 0;
    vector<int>::size_type j = 0;
    
    for(; i < first_community.size(); ++i, ++j)
        if(first_community[i] == first_state && second_community[i] == second_state)
            target[j] = 1;
    
    return target;
}
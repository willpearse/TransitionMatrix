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

vector<string> next_step(ublas::matrix<double> transition_matrix, vector<string> current_community, vector<string> species_names)
{
    vector<string> next_community = current_community; //Make a copy of the current_community
    boost::mt19937 rnd_generator(time(NULL)); //Make a random number generator
    int current_random; //Index we'll be using to figure out what species to add
    vector<string> new_slots;
    for(int i=0; i < current_community.size(); ++i) //Loop through each member of the current commnuity
    {
        for(int k=0; k < species_names.size(); ++k) //Loop through the possible species
        {
            if(current_community[i] == species_names[k])
            {
                boost::numeric::ublas::matrix_row<ublas::matrix<double> > matrix_row (transition_matrix, k);
                boost::random::discrete_distribution<> discrete_dist(matrix_row); //Make a discrete distribution
                current_random = discrete_dist(rnd_generator); //Bind the generator distribution to a random number generator
                if(species_names[current_random] != "REPRODUCE") //Assign the next slot, checking to handle reproduced taxa
                    next_community[i] = species_names[current_random];
                else
                    new_slots.push_back(current_community[i]); //Add the extra species we should be adding in
                    
                break; //Break out of the for loop
            }
        }
    }
    next_community.insert(next_community.end(), new_slots.begin(), new_slots.end()); //Add the newly created species into the community
    return next_community;
}

vector<double> likelihood(ublas::matrix<double> transition_matrix, vector<string> first_community, vector<string> second_community, vector<string> species_names)
{
    assert(species_names.size() == transition_matrix.size2()); //As many columns as species_names (that includes DEAD and REPRODUCE)
    assert(species_names.size() == (transition_matrix.size1() +2)); //Two fewer rows, because dead and reproduced species don't count
    
    int first_index, second_index, sp_index; //Initialise the indices we'll be using
    const unsigned REPRODUCTION = (species_names.size() - 2); //Reproduction is the penultiamte element in the matrix
    const unsigned DEATH = (species_names.size() - 1); //Death is the last element in the matrix
    vector<double> slot_likelihood(second_community.size(), -1), species_count(species_names.size()); //Make holders for each slot, and number of additions to be 'corrected' later when dealing with reproduction
    for(int i=(second_community.size()-1); i >= 0; --i) //Go along the second community - HACK!!!
    {
        if(i > (first_community.size()-1)) //We're dealing with reproduction - HACK!!!
        {
            for(sp_index = 0; sp_index < species_names.size(); ++sp_index) //Go along the species_names
                if(second_community[i] == species_names[sp_index]) //...if they match what we've got
                    break; //...stop looking
            assert(sp_index != -1); //We have found a match (nothing should be DEAD here)
            slot_likelihood[i] = transition_matrix(sp_index, REPRODUCTION);
            ++species_count[sp_index];
        }
        else //Now handle everything else
        {
            if(second_community[i] == "DEATH")
                slot_likelihood[i] = transition_matrix(sp_index, DEATH); //Death is easy to handle immediately
            else //Anything else requires a search
            {
                first_index=-1;
                second_index=-1;
                for(int k=0; k < species_names.size(); k++) //Go along the species' names; don't bother stopping the search early
                {
                    if(first_community[i] == species_names[k]) //Found the first species
                        first_index = k;
                    if(second_community[i] == species_names[k]) //Found the second species
                        second_index = k;
                }
                assert(first_index != -1 && second_index != -1); //We've found both
                if(species_names[first_index]==species_names[second_index] && species_count[first_index] > 0) //Do we have some reproduction to handle?
                {
                    slot_likelihood[i] = 1; //Set the likelihood to something safe
                    --species_count[first_index]; //Change the handle count
                } else //Otherwise store the right likelihood
                    slot_likelihood[i] = transition_matrix(first_index, second_index);
            }
        }
    }
    return slot_likelihood;
}

vector<double> likelihood(ublas::matrix<double> transition_matrix, ublas::matrix<int> event_matrix)
{
    //Setup
    int n_events=0,i,k,x=0;
    for(i=0; i<event_matrix.size1(); ++i)
        for(k=0; k<event_matrix.size2(); ++k)
            n_events += event_matrix(i,k);
    vector<double> likelihoods(n_events);
    
    //Loop through events and calculate likelihood
    // - check to see if double looping to pre-allocate is worth it
    for(i=0; i<event_matrix.size1(); ++i)
        for(k=0; k<event_matrix.size2(); ++k)
        {
            while(event_matrix(i,k)>0)
            {
                likelihoods[x++] = transition_matrix(i,k);
                --event_matrix(i,k);
            }
        }
    //Return and check
    assert(x == n_events);
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
    //If we're not going to reproduction (a nonsense) or dying, check repro.
    if(sp < t_m.size1())
        if(t_m(sp, ++i) > curr_max)
        {
            max_pos = i;
            curr_max = t_m(sp, i);
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
    const unsigned REPRODUCTION = (species_names.size() - 1);
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
                        --events_matrix(j,j);
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
        assert(finished);
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
//
//  communityData.cpp
//  transitionMatrix
//
//  Created by Will Pearse on 18/11/2011.
//  Copyright 2011 Imperial College London. All rights reserved.
//
//  Contains classes for community data loading and intialisation
//  Does *not* contain any likelihood calculation functions
//  These are in "communitySimulator" files
//  
//  Class 'Community' represents a single community (with multiple entries for multiple years)
//  Class 'Data' represents multiple communities
//  The file is sorted according to the type of operation, not by class
//  TO-DO: use polymorphism, etc., to neaten things a little
//

#include "communityData.h"

using namespace std;
using namespace boost::numeric;
typedef vector<string> str_vec;

//HOUSEKEEPING
static int string_to_int(string var)
{
    int var_i;
    stringstream convert;
    convert << var;
    convert >> var_i;
    return var_i;
}
static vector<string> fill_community(string species, int abundance)
{
    vector<string> output;
    while(abundance > 0)
    {
        output.push_back(species);
        --abundance;
    }
    return output;
}
static double log_add(double x, double y)
{
    return(x + log(y));
}
void Community::set_e_m(void)
{   
    //Clean this bullshit up and learn to preallocate with pointers...
    ublas::matrix<int> temp(species_names.size(), species_names.size()+2);
    for(int i = 0; i<(communities.size()-1); ++i)
        e_m.push_back(likely_transitions(t_m, communities[i], communities[i+1], sp_names_repro_death));
}
void Community::set_t_m(double stationary_prob)
{
    t_m.resize(species_names.size(), species_names.size()+2);
    double leftover_prob = (1 - stationary_prob) / (species_names.size() + 1);
    for(int i=0; i<t_m.size1(); ++i)
        for(int j=0; j<t_m.size2(); ++j)
            if(i==j)
                t_m(i,j) = stationary_prob;
            else
                t_m(i,j) = leftover_prob;
}
void Community::set_t_null(void)
{
    t_null.resize(species_names.size());
    vector<double> species_count(species_names.size(), 0);
    double total_count;
    
    for(int i=0; i<communities.size(); ++i)
        for(int j=0; j<communities[i].size(); ++j)
            for(int k=0; k<species_names.size(); ++k)
                if(species_names[k] == communities[i][j])
                {
                    ++species_count[k];
                    break;
                }
    
    total_count = accumulate(species_count.begin(), species_count.end(), 0);
    for(int i=0; i<species_count.size(); ++i)
        t_null[i] = species_count[i] / total_count;
}
void DataSet::set_t_m(double stationary_prob)
{
    t_m.resize(species_names.size(), species_names.size()+2);
    double leftover_prob = (1 - stationary_prob) / (species_names.size() + 1);
    for(int i=0; i<t_m.size1(); ++i)
        for(int j=0; j<t_m.size2(); ++j)
            if(i==j)
                t_m(i,j) = stationary_prob;
            else
                t_m(i,j) = leftover_prob;
}
void DataSet::set_t_null(void)
{
    t_null.resize(species_names.size());
    vector<double> species_count(species_names.size(), 0);
    double total_count;
    for(int x=0; x<communities.size(); ++x)
        for(int i=0; i<communities[x].communities.size(); ++i)
            for(int j=0; j<communities[x].communities[i].size(); ++j)
                for(int k=0; k<species_names.size(); ++k)
                    if(species_names[k] == communities[x].communities[i][j])
                    {
                        ++species_count[k];
                        break;
                    }
    
    total_count = accumulate(species_count.begin(), species_count.end(), 0);
    for(int i=0; i<species_count.size(); ++i)
        t_null[i] = species_count[i] / total_count;
}

////////////////
//CONSTRUCTORS//
////////////////

Community::Community(string species, string abundance, string year, string name)
{
    //Setup
    int abundance_int=string_to_int(abundance), year_int=string_to_int(year);
    vector<string> community=fill_community(species, abundance_int);
    
    //Make the starting values
    communities.push_back(community);
    n_species = 1;
    years.push_back(year_int);
    species_names.push_back(species);
    community_name = name;
}
Community::Community(ublas::matrix<double> transition_matrix, vector<string> sp_names, int n_communities, int starting_length, string name, int seed)
{
    //Setup
    assert(starting_length >= 1);
    vector<string> community_vec;
    int current_vec = 0;
    
    //Make first community from uniform distribution
    boost::mt19937 rnd_generator(seed);
    vector<double> uniform_weights(sp_names.size());
    for(int i=0; i<sp_names.size(); ++i)
        uniform_weights[i] = 1.0 / sp_names.size();
    boost::random::discrete_distribution<> uniform(uniform_weights.begin(), uniform_weights.end());
    for(int i=0; i<starting_length; ++i){
        int x = uniform(rnd_generator);
        community_vec.push_back(sp_names[x]);
    }
    
    sort(community_vec.begin(), community_vec.end());
    communities.push_back(community_vec);
    
    //Make generators for all the species
    vector<boost::random::discrete_distribution<> > species_distributions;
    vector<double> current_rates(sp_names.size()+2);
    for(int i=0; i<sp_names.size(); ++i)
    {
        for(int j=0; j<current_rates.size(); ++j)
            current_rates[j] = transition_matrix(i,j);
        boost::random::discrete_distribution<> curr_dist(current_rates.begin(), current_rates.end());
        species_distributions.push_back(curr_dist);
    }
    
    //Add in death and reproduction for the species
    species_names = sp_names;
    n_species = sp_names.size();
    sp_names.push_back("REPRODUCE");
    sp_names.push_back("DEATH");
    assert(sp_names.size() == transition_matrix.size2());
    sp_names_repro_death = sp_names;
    set_t_m();
    
    //Make real transition matrix
    real_t_m = transition_matrix;
    
    //Make next load of communities
    while(++current_vec < n_communities)
    {
        //Setup
        ublas::matrix<int> event_matrix = ublas::zero_matrix<int>(species_names.size(), species_names.size()+2);
        vector<string> current_com;
        //Go along all species
        for(int i=0; i<communities[current_vec-1].size(); ++i)
        {
            //Find current species
            int j=0;
            for(; j<sp_names.size(); ++j)
                if(sp_names[j]==communities[current_vec-1][i])
                    break;
            
            //We should have handled reproduction earlier
            assert(sp_names[j] != "REPRODUCE");
            
            //Choose next step and record it in the event matrix
            int next_step = species_distributions[j](rnd_generator);
            ++event_matrix(j,next_step);
            //Add the new species in the list, two if it's a reproduction, and nothing if it's death
            if(sp_names[next_step] == "REPRODUCE")
            {
                current_com.push_back(sp_names[j]);
                current_com.push_back(sp_names[j]);
            }
            else
                if(sp_names[next_step] != "DEATH")
                {
                    current_com.push_back(sp_names[next_step]);
                }
        }
        real_e_m.push_back(event_matrix);
        sort(current_com.begin(), current_com.end());
        communities.push_back(current_com);
        e_m.push_back(likely_transitions(t_m, communities[current_vec-1], current_com, sp_names_repro_death));
    }
    //Construct a community from all this data
    // - there's a lot of commonality between the constructors...
    for(int i=0; i<n_communities; ++i)
    {
        years.push_back(i);
        communities[i].erase(communities[i].begin());
    }
    community_name = name;
    set_t_null();
}
DataSet::DataSet(const char *file, double stationary_prob)
{
    //Setup
    ifstream str(file);
    string line,cell;
    vector<string> splits;
    string species,abundance,year,community_name;
    int i;
    
    //For each line
    while(getline(str, line))
    {   
        //Get abundance, species, community name and year
        stringstream stream(line);
        getline(stream, abundance, ',');
        getline(stream, species, ',');
        getline(stream, community_name, ',');
        getline(stream, year, ',');
        
        //Do we already have this community?
        for(i=0; i<community_names.size(); ++i)
            if(community_name == community_names[i])
                break;
        if(i == community_names.size())
        {
            //No, so make a new one and store its attributes
            Community temp(species, abundance, year, community_name);
            communities.push_back(temp);
            community_names.push_back(community_name);
        }
        else //Yes, so add to that particular community
            communities[i].add_species(species, abundance, year);
        
        //Add species to list if it's new
        for(i=0; i<species_names.size(); ++i)
            if(species == species_names[i])
                break;
        if(i == species_names.size())
            species_names.push_back(species);
    }
    
    //Initialise all communities
    for(i=0; i<communities.size(); ++i)
        communities[i].initialise();
    
    set_t_m(stationary_prob);
    set_t_null();
}

//INITIALISATION
void Community::initialise(void)
{
    //A dirty hack: load everything into a map, then pull it all out again and the map will have ordered it for us!
    //Setup
    map<int, vector<string> > temp_map;
    map<int, vector<string> >::iterator iter;
    int i;
    
    //Map assingment
    for(int i=0; i<communities.size(); ++i)
        temp_map[years[i]] = communities[i];
    
    //(Inefficiently) put everything in place
    i = 0;
    for(iter=temp_map.begin(); iter!=temp_map.end(); ++iter)
    {
        years[i] = iter->first;
        communities[i] = iter->second;
        ++i;
    }
    
    //Set species names for calculations
    sp_names_repro_death = species_names;
    sp_names_repro_death.push_back("REPRODUCE");
    sp_names_repro_death.push_back("DEATH");
    
    //Make a rough transition matrix
    set_t_m();
    
    //Make event matrices for everything
    set_e_m();
    
    //Make the null transition matrix
    set_t_null();
}
void Community::add_species(string species, string abundance, string year)
{
    //Setup
    int abundance_int=string_to_int(abundance), year_int=string_to_int(year),i=0;
    vector<string> community=fill_community(species, abundance_int);
    vector<string>::iterator iter;
    
    //Do we already have a community for this year?
    // - can't use find without an iterator, can't use interators as indices --> write own search function
    for(; i<years.size(); ++i)
    {
        if(years[i] == year_int)
        {
            communities[i].insert(communities[i].end(), community.begin(), community.end());
            break;
        }
    }
    
    //If we didn't find a match, we should make a new year
    if(i==years.size())
    {
        communities.push_back(community);
        years.push_back(year_int);
    }
    
    //Insert any new names
    iter = find(species_names.begin(), species_names.end(), species);
    if(iter == species_names.end())
        species_names.push_back(species);
}

//OPTIMISATION
double likelihood_parameter_com(double prob, int row, int column, vector<ublas::matrix<int> > e_m, ublas::matrix<double> t_m, int minimising)
{
    //Setup
    vector<double> likelihoods(e_m.size()), current_likelihood(e_m.size());
    int i;
    double total_likelihood,fudge_factor,leftover;
    
    //Set the new probability in the transition_matrix
    leftover = 1 - prob;
    fudge_factor = 1 - t_m(row,column);
    for(i=0; i<t_m.size2(); ++i)
        t_m(row,i) = (t_m(row,i) / fudge_factor) * leftover;
    t_m(row,column) = prob;
    
    //Calculate likelihood
    for(i=0; i < e_m.size(); ++i)
    {
        current_likelihood = likelihood(t_m, e_m[i]);
        
        for(int j=0; j<current_likelihood.size(); ++j)
            current_likelihood[j] = log(current_likelihood[j]);
        
        likelihoods[i] = accumulate(current_likelihood.begin(), current_likelihood.end(), 0.0);
        current_likelihood.clear();
    }
    
    //Return inverse of log likelihood if passing to a minimising function
    total_likelihood = accumulate(likelihoods.begin(), likelihoods.end(), 0.0);
    if(minimising)
        return (0-total_likelihood);
    else
        return total_likelihood;
}
double likelihood_parameter_data(double prob, int row, int column, vector<Community> communities, ublas::matrix<double> t_m, int minimising)
{
    //Setup
    vector<double> likelihoods(communities.size());
    int i;
    double total_likelihood,fudge_factor,leftover;
    
    //Set the new probability in the transition_matrix
    leftover = 1 - prob;
    fudge_factor = 1 - t_m(row,column);
    for(i=0; i<t_m.size2(); ++i)
        t_m(row,i) = (t_m(row,i) / fudge_factor) * leftover;
    t_m(row,column) = prob;
    
    //Calculate likelihood for each community
    for(int com=0; com<communities.size(); ++com)
    {
        vector<double> com_likelihoods(communities[com].e_m.size());
        for(i=0; i < communities[com].e_m.size(); ++i)
        {
            vector<double> current_likelihood;
            current_likelihood = likelihood(t_m, communities[com].e_m[i]);
            
            for(int j=0; j<current_likelihood.size(); ++j)
                current_likelihood[j] = log(current_likelihood[j]);
            
            com_likelihoods[i] = accumulate(current_likelihood.begin(), current_likelihood.end(), 0.0);
        }
        likelihoods[com] = accumulate(com_likelihoods.begin(), com_likelihoods.end(), 0.0);
    }
    
    //Return inverse of log likelihood if passing to a minimising function
    total_likelihood = accumulate(likelihoods.begin(), likelihoods.end(), 0.0);
    if(minimising)
        return (0-total_likelihood);
    else
        return total_likelihood;
}
double likelihood_null_parameter_com(double prob, int species, vector<vector<string> > communities, vector<double> t_null, vector<string> species_names, int minimising)
{
    //Setup
    vector<double> likelihoods(communities.size()), current_likelihood(communities.size());
    double total_likelihood,fudge_factor,leftover;
    
    //Set the new probability in the transition_matrix
    leftover = 1 - prob;
    fudge_factor = 1 - t_null[species];
    for(int i=0; i<t_null.size(); ++i)
        t_null[i] = (t_null[i] / fudge_factor) * leftover;
    t_null[species] = prob;
    
    //Calculate likelihood
    for(int i=0; i < communities.size(); ++i)
    {
        current_likelihood = likelihood_null(t_null, communities[i], species_names);
        
        for(int j=0; j<current_likelihood.size(); ++j)
            current_likelihood[j] = log(current_likelihood[j]);
        
        likelihoods[i] = accumulate(current_likelihood.begin(), current_likelihood.end(), 0.0);
        current_likelihood.clear();
    }
    
    //Return inverse of log likelihood if passing to a minimising function
    total_likelihood = accumulate(likelihoods.begin(), likelihoods.end(), 0.0);
    if(minimising)
        return (0-total_likelihood);
    else
        return total_likelihood;
}
double likelihood_null_parameter_data(double prob, int species, vector<Community> communities, vector<double> t_null, vector<string> species_names, int minimising)
{
    //Setup
    vector<double> likelihoods(communities.size()), current_likelihood(communities.size());
    double total_likelihood,fudge_factor,leftover;
    
    //Set the new probability in the transition_matrix
    leftover = 1 - prob;
    fudge_factor = 1 - t_null[species];
    for(int i=0; i<t_null.size(); ++i)
        t_null[i] = (t_null[i] / fudge_factor) * leftover;
    t_null[species] = prob;
    
    //Calculate likelihood for each community
    for(int com=0; com<communities.size(); ++com)
    {
        vector<double> com_likelihoods(communities[com].communities.size());
        for(int i=0; i < communities[com].communities.size(); ++i)
        {
            vector<double> current_likelihood;
            current_likelihood = likelihood_null(t_null, communities[com].communities[i], species_names);
            
            for(int j=0; j<current_likelihood.size(); ++j)
                current_likelihood[j] = log(current_likelihood[j]);
            
            com_likelihoods[i] = accumulate(current_likelihood.begin(), current_likelihood.end(), 0.0);
        }
        likelihoods[com] = accumulate(com_likelihoods.begin(), com_likelihoods.end(), 0.0);
    }
    
    //Return inverse of log likelihood if passing to a minimising function
    total_likelihood = accumulate(likelihoods.begin(), likelihoods.end(), 0.0);
    if(minimising)
        return (0-total_likelihood);
    else
        return total_likelihood;
}
double Community::optimise(int max_iter)
{
    //Setup
    double sum_of_parameters=0.0;
    //Loop through each parameter
    for(row=0; row<t_m.size1(); ++row)
    {
        for(int column=0; column<t_m.size2(); ++column)
        {
            //Are we dealing with the last parameter?
            if(column == (t_m.size2()-1))
            {
                // - a dodgy way to set the last parameter!
                t_m(row,column) = 1.0 - sum_of_parameters;
                sum_of_parameters = 0.0;
            }
            else
            {
                t_m(row,column) = boost::math::tools::brent_find_minima(boost::bind(likelihood_parameter_com, _1, row,column,e_m,t_m,1), 0.0, 1.0, max_iter).first;
                sum_of_parameters += t_m(row,column);
            }
        }
    }
    return calc_likelihoods();
}
double Community::optimise_null(int max_iter)
{
    //Loop through each parameter
    for(species=0; species<t_null.size(); ++species)
        t_null[species] = boost::math::tools::brent_find_minima(boost::bind(likelihood_null_parameter_com, _1, species,communities,t_null,species_names, 1), 0.0, 1.0, max_iter).first;
    return calc_null_likelihoods();
}
void DataSet::optimise(int max_iter)
{
    //Setup
    double sum_of_parameters=0.0;
    vector<double> likelihoods(communities.size());
    //Loop through each parameter
    for(row=0; row<t_m.size1(); ++row)
    {
        for(column=0; column<t_m.size2(); ++column)
        {
            //Are we dealing with the last parameter?
            if(column == (t_m.size2()-1))
            {
                // - a dodgy way to set the last parameter!
                t_m(row,column) = 1.0 - sum_of_parameters;
                sum_of_parameters = 0.0;
            }
            else
            {
                t_m(row,column) = boost::math::tools::brent_find_minima(boost::bind(likelihood_parameter_data, _1, row,column,communities,t_m,1), 0.0, 1.0, max_iter).first;
                cout << t_m(row,column) << ",";
                sum_of_parameters += t_m(row,column);
            }
        }
    }
}
void DataSet::optimise_null(int max_iter)
{
    vector<double> likelihoods(communities.size());
    for(species=0; species<t_null.size(); ++species)
        t_null[species] = boost::math::tools::brent_find_minima(boost::bind(likelihood_null_parameter_data, _1, species,communities,t_null,species_names,1), 0.0, 1.0, max_iter).first;
}

/////////////////
//LIKELIHOOD/////
/////////////////
double Community::calc_likelihoods(void)
{
    //Setup
    double total_likelihood;
    vector<double> likelihoods(e_m.size()), current_likelihood;
    
    //Loop through and calculate likelihoods
    for(int i =0; i < e_m.size(); ++i)
    {
        current_likelihood = likelihood(t_m, e_m[i]);
        likelihoods[i] = accumulate(current_likelihood.begin(), current_likelihood.end(), 0.0, log_add);
        current_likelihood.clear();
    }
    
    //Accumulate and return total likelihood
    total_likelihood = accumulate(likelihoods.begin(), likelihoods.end(), 0.0);
    return total_likelihood;
}
double DataSet::calc_likelihoods(void)
{
    //Setup
    double total_likelihood;
    vector<double> likelihoods(communities.size());
    
    //Loop through and calculate likelihoods
    for(int i =0; i < communities.size(); ++i)
        likelihoods[i] = communities[i].calc_likelihoods();
    
    //Accumulate and return total likelihood
    total_likelihood = accumulate(likelihoods.begin(), likelihoods.end(), 0.0);
    return total_likelihood;
}
double Community::calc_null_likelihoods(void)
{
    //Setup
    double total_likelihood;
    vector<double> likelihoods(communities.size()), current_likelihood;
    
    //Loop through and calculate likelihoods
    for(int i =0; i < communities.size(); ++i)
    {
        current_likelihood = likelihood_null(t_null, communities[i], species_names);
        likelihoods[i] = accumulate(current_likelihood.begin(), current_likelihood.end(), 0.0, log_add);
        current_likelihood.clear();
    }
    
    //Accumulate and return total likelihood
    total_likelihood = accumulate(likelihoods.begin(), likelihoods.end(), 0.0);
    return total_likelihood;
}
double DataSet::calc_null_likelihoods(void)
{
    //Setup
    double total_likelihood;
    vector<double> likelihoods(communities.size());
    
    //Loop through and calculate likelihoods
    for(int i =0; i < communities.size(); ++i)
        likelihoods[i] = communities[i].calc_null_likelihoods();
    
    //Accumulate and return total likelihood
    total_likelihood = accumulate(likelihoods.begin(), likelihoods.end(), 0.0);
    return total_likelihood;
}
//ACCESS
vector<string> Community::extract_year(int year)
{
    int i=0;
    for(; i<years.size(); ++i)
        if(years[i] == year)
            break;
    assert(i < years.size());
    return communities[i];
}
vector<string> Community::extract_index(int index)
{
    assert(index < communities.size());
    return communities[index];
}
string Community::name(void){
    return community_name;}
vector<int>::size_type Community::n_years(void){
    return years.size();}

//DISPLAY
boost::numeric::ublas::matrix<double> Community::print_transition_matrix(int width)
{
    //Setup
    int i,j;
    
    //Header
    cout << endl << setw(width) << "" ;
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        cout << setw(width) << *iter;
    cout << setw(width) << "Repro." << setw(width) << "Death" << endl;
    
    //Looping through
    for(i = 0; i<t_m.size1(); ++i)
    {
        cout << setw(width) << species_names[i];
        for(j=0; j<t_m.size2(); ++j)
            if(t_m(i,j)>0.0001)
                cout << setw(width) << setprecision(4) << t_m(i,j);
            else
                cout << setw(width) << setprecision(4) << 0;
        cout << endl;
    }
    
    return t_m;
}
boost::numeric::ublas::matrix<int> Community::print_real_transition_matrix(int width)
{
    //Header
    cout << endl << setw(width) << "" ;
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        cout << setw(width) << *iter;
    cout << setw(width) << "Repro." << setw(width) << "Death" << endl;
    
    //Looping through
    for(int i = 0; i<real_t_m.size1(); ++i)
    {
        cout << setw(width) << species_names[i];
        for(int j=0; j<real_t_m.size2(); ++j)
            cout << setw(width) << real_t_m(i,j);
        cout << endl;
    }
    
    return real_t_m;
}
boost::numeric::ublas::matrix<double> DataSet::print_transition_matrix(int width)
{
    //Header
    cout << endl << setw(width) << "" ;
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        cout << setw(width) << *iter;
    cout << setw(width) << "Repro." << setw(width) << "Death" << endl;
    
    //Looping through
    for(int i = 0; i<t_m.size1(); ++i)
    {
        cout << setw(width) << species_names[i];
        for(int j=0; j<t_m.size2(); ++j)
            if(t_m(i,j)>0.0001)
                cout << setw(width) << setprecision(4) << t_m(i,j);
            else
                cout << setw(width) << setprecision(4) << 0;
        cout << endl;
    }
    
    return t_m;
}
vector<double> Community::print_null_transition_vector(int width)
{
    //Header
    cout << endl;
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        cout << setw(width) << *iter;
    cout << endl;
    
    //Looping through
    for(int i = 0; i<t_null.size(); ++i)
    {
        if(t_null[i]>0.0001)
            cout << setw(width) << setprecision(4) << t_null[i];
        else
            cout << setw(width) << setprecision(4) << 0;
    }
    cout << endl;   
    return t_null;
}
vector<double> DataSet::print_null_transition_vector(int width)
{
    //Header
    cout << endl;
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        cout << setw(width) << *iter;
    cout << endl;
    
    //Looping through
    for(int i = 0; i<t_null.size(); ++i)
    {
        if(t_null[i]>0.0001)
            cout << setw(width) << setprecision(4) << t_null[i];
        else
            cout << setw(width) << setprecision(4) << 0;
    }
    cout << endl;   
    return t_null;
}
vector<string> Community::print_community(int index, int width)
{
    assert(index < communities.size());
    cout << endl;
    for(vector<string>::const_iterator iter = communities[index].begin(); iter != communities[index].end(); ++iter)
        cout << setw(width) << *iter;
    cout <<endl;
    return communities[index];
}
boost::numeric::ublas::matrix<int> Community::print_event_matrix(int index, int width)
{
    //Setup
    assert(index < e_m.size());
    boost::numeric::ublas::matrix<int> curr_e_m = e_m[index];
    
    //Header
    cout << endl << setw(width) << "" ;
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        cout << setw(width) << *iter;
    cout << setw(width) << "Repro." << setw(width) << "Death" << endl;
    
    //Looping through
    for(int i = 0; i<curr_e_m.size1(); ++i)
    {
        cout << setw(width) << species_names[i];
        for(int j=0; j<curr_e_m.size2(); ++j)
            cout << setw(width) << curr_e_m(i,j);
        cout << endl;
    }
    
    return curr_e_m;
}
void Community::write_transition_matrix(const char* file_name)
{
    //Setup
    ofstream file(file_name);
    
    //Header
    file << calc_likelihoods();
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        file << "," << *iter;
    file << ",Reproduction,Death" << endl;
    
    //Looping through
    for(int i = 0; i<t_m.size1(); ++i)
    {
        file << species_names[i];
        for(int j=0; j<t_m.size2(); ++j)
            if(t_m(i,j)>0.0001)
                file << "," << setprecision(4) << t_m(i,j);
            else
                file << "," << 0;
        file << endl;
    }
    file.close();
}
void DataSet::write_transition_matrix(const char* file_name)
{
    //Setup
    ofstream file(file_name);
    
    //Header
    file << calc_likelihoods();
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        file << "," << *iter;
    file << ",Reproduction,Death" << endl;
    
    //Looping through
    for(int i = 0; i<t_m.size1(); ++i)
    {
        file << species_names[i];
        for(int j=0; j<t_m.size2(); ++j)
            if(t_m(i,j)>0.0001)
                file << "," << setprecision(4) << t_m(i,j);
            else
                file << "," << 0;
        file << endl;
    }
    file.close();
}
void DataSet::write_null_transition_vector(const char* file_name)
{
    //Setup
    ofstream file(file_name);
    
    //Header
    file << calc_null_likelihoods();
    //Header
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        file << "," << *iter;
    file << endl;
    
    //Looping through
    file << setprecision(4) << t_null[0];
    for(int i = 1; i<t_null.size(); ++i)
    {
        if(t_null[i]>0.0001)
            file << "," << setprecision(4) << t_null[i];
        else
            file << "," << setprecision(4) << 0;
    }
}
boost::numeric::ublas::matrix<int> Community::print_real_event_matrix(int index, int width)
{
    //Setup
    assert(index < real_e_m.size());
    boost::numeric::ublas::matrix<int> curr_e_m = real_e_m[index];
    
    //Header
    cout << endl << setw(width) << "" ;
    for(vector<string>::const_iterator iter = species_names.begin(); iter != species_names.end(); ++iter)
        cout << setw(width) << *iter;
    cout << setw(width) << "Repro." << setw(width) << "Death" << endl;
    
    //Looping through
    for(int i = 0; i<curr_e_m.size1(); ++i)
    {
        cout << setw(width) << species_names[i];
        for(int j=0; j<curr_e_m.size2(); ++j)
            cout << setw(width) << curr_e_m(i,j);
        cout << endl;
    }
    
    return curr_e_m;
}

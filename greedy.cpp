/***************************************************************************
                          greedy.cpp  -  description
                             -------------------
    begin                : March 27 2024
    copyright            : (C) 2024 by Christian Blum
    email                : christian.blum@iiia.csic.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

using namespace std;

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <limits>
#include <random>
#include <chrono>

struct Option {
    int vertex;
    double value;
};

struct Solution {
    set<int> vertices;
    int score;
};

// instance data
int n_of_vertices;
vector< set<int> > neigh;                                   //contains neighbors of each vertex
vector< set<int> > no_neigh;

string input_file;



inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

int produce_random_integer(int max, double& rnum) {

    int num = int(double(max)*rnum);
    if (num == max) num = num - 1;
    return num;
}

int get_random_element(const set<int>& s, double& rnum) {                 //for tie-breaking max_vertices

    double r = produce_random_integer(int(s.size()), rnum);
    set<int>::iterator it = s.begin();
    advance(it, r);
    return *it;
}

void read_parameters(int argc, char **argv) {

    int iarg=1;
    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) input_file = argv[++iarg];
        iarg++;
    }
}

void generate_greedy_solution(Solution& greedy_sol, default_random_engine& generator, uniform_real_distribution<double>& standard_distribution) {

    greedy_sol.score = 0;
    vector<bool> covered(n_of_vertices, false);
    vector<bool> complement_covered(n_of_vertices, false);
    int num_nodes_uncovered = n_of_vertices;
    int num_nodes_complement_uncovered = n_of_vertices;
    vector< set<int> > uncovered_neighbors = neigh;
    vector< set<int> > complement_uncovered_neighbors = no_neigh;
    
    set<int> candidates;
    for (int i = 0; i < n_of_vertices; ++i) candidates.insert(i);

    while (num_nodes_uncovered > 0 or num_nodes_complement_uncovered > 0) {
        double max_val = -1.0;
        set<int> max_vertices;
        
        //generate all options
        for (set<int>::iterator cit = candidates.begin(); cit != candidates.end(); ++cit) {
            Option opt;
            opt.vertex = *cit;
            double min = double(uncovered_neighbors[*cit].size());
            if (double(complement_uncovered_neighbors[*cit].size()) < min) min = double(complement_uncovered_neighbors[*cit].size());
            opt.value = min;
            if (opt.value >= max_val) {
                if (opt.value > max_val) {
                    max_val = opt.value;
                    max_vertices.clear();
                }
                max_vertices.insert(*cit);
            }
        }
        
        double rnum = standard_distribution(generator);
        int chosen_vertex = get_random_element(max_vertices, rnum);
        (greedy_sol.vertices).insert(chosen_vertex);
        greedy_sol.score += 1;

        if (not covered[chosen_vertex]) {
            covered[chosen_vertex] = true;
            --num_nodes_uncovered;
        }
        if (not complement_covered[chosen_vertex]) {
            complement_covered[chosen_vertex] = true;
            --num_nodes_complement_uncovered;
        }

        num_nodes_uncovered -= int(uncovered_neighbors[chosen_vertex].size());
        set<int> to_cover = uncovered_neighbors[chosen_vertex];
        for (set<int>::iterator sit = to_cover.begin(); sit != to_cover.end(); ++sit) {
            covered[*sit] = true;
            for (set<int>::iterator ssit = neigh[*sit].begin(); ssit != neigh[*sit].end(); ssit++) uncovered_neighbors[*ssit].erase(*sit);
        }
        uncovered_neighbors[chosen_vertex].clear();

        num_nodes_complement_uncovered -= int(complement_uncovered_neighbors[chosen_vertex].size());
        to_cover = complement_uncovered_neighbors[chosen_vertex];
        for (set<int>::iterator sit = to_cover.begin(); sit != to_cover.end(); ++sit) {
            complement_covered[*sit] = true;
            for (set<int>::iterator ssit = no_neigh[*sit].begin(); ssit != no_neigh[*sit].end(); ssit++) complement_uncovered_neighbors[*ssit].erase(*sit);
        }
        complement_uncovered_neighbors[chosen_vertex].clear();

        for (set<int>::iterator sit = neigh[chosen_vertex].begin(); sit != neigh[chosen_vertex].end(); sit++) uncovered_neighbors[*sit].erase(chosen_vertex);
        for (set<int>::iterator sit = no_neigh[chosen_vertex].begin(); sit != no_neigh[chosen_vertex].end(); sit++) complement_uncovered_neighbors[*sit].erase(chosen_vertex);
        candidates.erase(chosen_vertex);

        set<int> to_delete;
        for (set<int>::iterator cit = candidates.begin(); cit != candidates.end(); ++cit) {
            if (int(uncovered_neighbors[*cit].size()) == 0 and covered[*cit] and int(complement_uncovered_neighbors[*cit].size()) == 0 and complement_covered[*cit]) to_delete.insert(*cit);
        }
        for (set<int>::iterator sit = to_delete.begin(); sit != to_delete.end(); ++sit) candidates.erase(*sit);
    }
}


/**********
Main function
**********/

int main( int argc, char **argv ) {

    if ( argc < 3 ) {
        cout << "Use: cplex -i <input_file> ..." << endl;
        exit(1);
    }
    else read_parameters(argc,argv);

    std::cout << std::setprecision(2) << std::fixed;

    // initializes the random number generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> standard_distribution(0.0,1.0);

    ifstream indata;
    indata.open(input_file.c_str());
    if(!indata) { // file couldn't be opened
            cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_vertices;
    neigh = vector< set<int> >(n_of_vertices);          //neighbors data container

    set<int> complement;
    for (int i = 0; i < n_of_vertices; ++i) complement.insert(i);
    no_neigh = vector< set<int> >(n_of_vertices, complement);
    for (int i = 0; i < n_of_vertices; ++i) no_neigh[i].erase(i);

    int v1, v2;
    while (indata >> v1 >> v2) {
        neigh[v1].insert(v2);
        no_neigh[v1].erase(v2);
        neigh[v2].insert(v1);
        no_neigh[v2].erase(v1);
    }
    indata.close();

    // the computation time starts now
    clock_t start = clock();

    Solution greedy_sol;
    generate_greedy_solution(greedy_sol, generator, standard_distribution);
        
    // take the current time that has passed
    clock_t current = clock();
    double ctime = double(current - start) / CLOCKS_PER_SEC;

    cout << greedy_sol.score << "\t" << ctime << endl;

}


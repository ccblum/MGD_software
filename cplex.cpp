/***************************************************************************
                          cplex.cpp  -  description
                             -------------------
    begin                : Mar 19 2024
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <list>
#include <set>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <limits>
#include <random>
#include <chrono>
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN


vector< set<int> > neigh;
vector< set<int> > no_neigh;

double t_limit = 7200.0;

// instance data
int n_of_vertices;
vector<string> input_files;


inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg=1;
    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) input_files.push_back(argv[++iarg]);
        else if (strcmp(argv[iarg],"-t")==0) t_limit = atof(argv[++iarg]);
        iarg++;
    }
}

ILOSOLVECALLBACK6(loggingCallback,
                    clock_t&, start,
                    vector<double>&, results,
                    vector<double>&, times,
                    vector<double>&, gaps,
                    IloNum,         lastIncumbent,
                    int, iter)
{
    if (hasIncumbent()) {
        IloNum nv = getIncumbentObjValue();

        // take the current time that has passed
        clock_t current = clock();
        double newTime = double(current - start) / CLOCKS_PER_SEC;
        double newGap = 100.0*getMIPRelativeGap();
        if (lastIncumbent > nv) {
            cout << "value " << nv << "\ttime " << newTime << "\tgap " << newGap << endl;
            results[iter] = double(nv);
            times[iter] = newTime;
            gaps[iter] = newGap;
        }
        lastIncumbent = nv;
    }
}

void compute_optimal(set<int>& cpl_sol, clock_t& start, vector<double>& results, vector<double>& gaps, vector<double>& times, int iter) {

    IloEnv env;
    env.setOut(env.getNullStream());

    IloModel model(env);

    IloNumVarArray x(env, n_of_vertices, 0, 1, ILOINT);

    IloExpr obj(env);
    for (int i = 0; i < n_of_vertices; ++i) obj += x[i];
    model.add(IloMinimize(env, obj));
    obj.end();

    for (int i = 0; i < n_of_vertices; ++i) {
        IloExpr expr(env);
        expr += x[i];
        for (set<int>::iterator sit = neigh[i].begin(); sit != neigh[i].end(); ++sit) expr += x[*sit];
        model.add(expr >= 1);
        expr.end();

        IloExpr expr2(env);
        expr2 += x[i];
        for (set<int>::iterator sit = no_neigh[i].begin(); sit != no_neigh[i].end(); ++sit) expr2 += x[*sit];

        model.add(expr2 >= 1);
        expr2.end();
    }

    IloCplex cpl(model);

    cpl.setParam(IloCplex::TiLim, t_limit);
    cpl.setParam(IloCplex::EpGap, 0.0);
    cpl.setParam(IloCplex::EpAGap, 0.0);
    cpl.setParam(IloCplex::Threads, 1);
    cpl.setWarning(env.getNullStream());
    IloNum lastObjVal = std::numeric_limits<double>::max();
    cpl.use(loggingCallback(env, start, results, times, gaps, lastObjVal, iter));

    cpl.solve();

    if (cpl.getStatus() == IloAlgorithm::Optimal or cpl.getStatus() == IloAlgorithm::Feasible) {
        double lastVal = double(cpl.getObjValue());

        clock_t current = clock();
        double newTime = double(current - start) / CLOCKS_PER_SEC;
        double lastGap = 100.0*cpl.getMIPRelativeGap();

        if (lastVal < results[iter] or ((lastVal == results[iter]) and (lastGap < gaps[iter]))) {
            cout << "value " << lastVal << "\ttime " << newTime << "\tgap " << lastGap << endl;
            results[iter] = lastVal;
            times[iter] = newTime;
            gaps[iter] = lastGap;
        }
        if (cpl.getStatus() == IloAlgorithm::Optimal) {
            cout << "optimality proven" << endl;
        }

        IloNumArray x_val(env, n_of_vertices);
        cpl.getValues(x_val, x);

        int cplexsumnodes = 0;
        for (int i = 0; i < n_of_vertices; ++i) {
            if (x_val[i] > 0.9) {
                cpl_sol.insert(i);
                //cout << i << " ";
                ++cplexsumnodes;
            }
        }
    }

    env.end();

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

    std::cout << std::setprecision(3) << std::fixed;

    int n_fil = int(input_files.size());

    vector<double> results(n_fil, std::numeric_limits<double>::max());
    vector<double> times(n_fil, 0.0);
    vector<double> gaps(n_fil, 0.0);

    for (int iter = 0; iter < n_fil; ++iter) {

        ifstream indata;
        indata.open(input_files[iter].c_str());
        if(!indata) { // file couldn't be opened
            cout << "Error: file could not be opened" << endl;
        }

        indata >> n_of_vertices;

        neigh = vector< set<int> >(n_of_vertices);
        
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
        
        cout << "start file " << input_files[iter] << endl;

        set<int> cpl_sol;

        compute_optimal(cpl_sol, start, results, gaps, times, iter);

        cout << "end file " << input_files[iter] << endl;
    }

    double r_mean = 0.0;
    double g_mean = 0.0;
    double t_mean = 0.0;
    for (int i = 0; i < results.size(); i++) {
        r_mean = r_mean + results[i];
        g_mean = g_mean + gaps[i];
        t_mean = t_mean + times[i];
    }
    r_mean = r_mean / ((double)results.size());
    g_mean = g_mean / ((double)gaps.size());
    t_mean = t_mean / ((double)times.size());
    double rsd = 0.0;
    double gsd = 0.0;
    double tsd = 0.0;
    for (int i = 0; i < results.size(); i++) {
        rsd = rsd + pow(results[i]-r_mean,2.0);
        gsd = gsd + pow(gaps[i]-g_mean,2.0);
        tsd = tsd + pow(times[i]-t_mean,2.0);
    }
    rsd = rsd / ((double)(results.size()-1.0));
    if (rsd > 0.0) {
        rsd = sqrt(rsd);
    }
    gsd = gsd / ((double)(gaps.size()-1.0));
    if (gsd > 0.0) {
        gsd = sqrt(gsd);
    }
    tsd = tsd / ((double)(times.size()-1.0));
    if (tsd > 0.0) {
        tsd = sqrt(tsd);
    }

    cout << r_mean << "\t" << rsd << "\t" << t_mean << "\t" << tsd << "\t" << g_mean << "\t" << gsd<< endl;
}


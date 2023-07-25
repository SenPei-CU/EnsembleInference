#include <iostream>
#include <deque>
#include <time.h>
#include<random>
#include <algorithm>
#include <math.h>
#include<fstream>
#include<stdio.h>
#include<cassert>
#include <stdlib.h>
#pragma warning(disable : 4996)
using namespace std;
// Generate a fixed range of random numbers
double MyRand(double minValue, double maxValue)
{
    double Value;
    Value = (maxValue - minValue) * ((double)rand() / RAND_MAX) + minValue;
    return Value;
}
//Get edge set
void Get_Real_Network(int* NodeNum, int* EdgeNum, deque< deque<int> >* Link) {
    int i;
    int max = 0;
    int edge[200000];
    ifstream fin("hamster_1.txt");
    i = 0;
    while (fin >> edge[i]) {
        i++;
    }
    if (i % 2 == 0) {
        (*EdgeNum) = i / 2;;
    }
    else {
        cout << "Real Network has error\n" << endl;
        return;
    }
    for (int j = 0; j < (*EdgeNum); j++) {
        deque<int> unit;
        if (edge[2 * j] != edge[2 * j + 1]) {
            unit.push_back(edge[2 * j]);
            unit.push_back(edge[2 * j + 1]);
            if (unit[0] > max) {
                max = unit[0];
            }
            if (unit[1] > max) {
                max = unit[1];
            }
            (*Link).push_back(unit);
        }
    }
    (*NodeNum) = max + 1;
    (*EdgeNum) = (int)(*Link).size();
    return;
}
//Get network structure
void Real_Network_Transition(int NodeNum, int* EdgeNum, deque< deque<int> >* Neighbor, deque< deque<int> >* Link) {
    int inter;
    for (int i = 0; i < (*EdgeNum); i++) {
        if ((*Link)[i][0] > (*Link)[i][1]) {
            inter = (*Link)[i][0];
            (*Link)[i][0] = (*Link)[i][1];
            (*Link)[i][1] = inter;
        }
    }
    for (int i = 0; i < (*EdgeNum); i++) {
        for (int j = i + 1; j < (*EdgeNum); j++) {
            if ((*Link)[j][0] < (*Link)[i][0]) {
                (*Link)[i].swap((*Link)[j]);
            }
        }
    }
    int iter1 = 0;
    int iter2 = 0;
    int rep_Edge = 0;
    while (iter1 != (*EdgeNum) && iter2 != (*EdgeNum)) {//
        do {
            iter2++;
            if (iter2 == (*EdgeNum)) {
                break;
            }
        } while ((*Link)[iter2][0] == (*Link)[iter1][0]);

        for (int i = iter1; i < iter2; ++i) {
            for (int j = i + 1; j < iter2; ++j) {
                if ((*Link)[j][1] < (*Link)[i][1]) {
                    inter = (*Link)[i][1];
                    (*Link)[i][1] = (*Link)[j][1];
                    (*Link)[j][1] = inter;
                }
                else if ((*Link)[j][1] == (*Link)[i][1]) {
                    (*Link)[j][1] = -1;
                    rep_Edge++;
                }
            }
        }
        iter1 = iter2;
    }

    for (int i = 0; i < NodeNum; i++) {
        deque<int> iter;
        (*Neighbor).push_back(iter);
    }
    int first, second;
    for (int i = 0; i < (*EdgeNum); i++) {
        if ((*Link)[i][1] != -1) {
            first = (*Link)[i][0];
            second = (*Link)[i][1];
            (*Neighbor)[first].push_back(second);
            (*Neighbor)[second].push_back(first);
        }
    }
    (*EdgeNum) = (*EdgeNum) - rep_Edge;
    return;
}
//Get propagation probability for agents
void Get_beta_set(deque<double>* beta_set, int N, double* beta) {
    int i = 0;
    ifstream fin("beta_hum874_b.txt");//1*N
    double message[20000];
    while (fin >> message[i]) {
        i++;
    }
    if (i != N) {
        cout << "beta set.txt has error\n" << endl;
        cout << "i=" << i << endl;
        exit(100);
    }
    double sum = 0;
    for (int j = 0; j < i; j++) {
        (*beta_set).push_back(message[j]);
        sum += message[j];
    }
    (*beta) = sum / N;
    return;
}
// Randomly extract a fixed number of seeds
deque<int>Get_seed(int seed_num, deque< deque<int> > Neighbor, int N) {
    deque<int> seed_set;
    deque<int>check(N, 0);
    for (int j = 0; j < seed_num; j++) {
    Lab: int v = rand() % N;
        if (Neighbor[v].size() <= 1 || check[v] == 1)
            goto Lab;
        seed_set.push_back(v);
        check[v] = 1;
    }
    return seed_set;
}
// Simulate propagation based on seed and network structure
deque<deque<int>>Simulation(deque< deque<int> > Neighbor, int N, deque<double>beta_set, double gama, int tstop, deque<int>seed_set) {

    deque<int> infect_prior;
    deque<int> infect_current;
    deque <int> State_t;//1*N
    deque<deque<int>> State;//tstop*N

    //Initialization
    infect_prior = seed_set;
    for (int i = 0; i < N; i++) {
        State_t.push_back(0);
    }
    for (int i = 0; i < seed_set.size(); i++) {
        State_t[seed_set[i]] = 1;
    }
    State.push_back(State_t);

    //Simulation
    for (int i = 0; i < tstop; i++) {
        for (int j = 0; j < infect_prior.size(); j++) {
            for (int k = 0; k < Neighbor[infect_prior[j]].size(); k++) {
                if (State_t[Neighbor[infect_prior[j]][k]] == 0) {
                    double r1 = (double)rand() / RAND_MAX;
                    if (r1 < beta_set[Neighbor[infect_prior[j]][k]]) {
                        State_t[Neighbor[infect_prior[j]][k]] = 1;//Infection
                        infect_current.push_back(Neighbor[infect_prior[j]][k]);
                    }
                }
            }
            double r2 = (double)rand() / RAND_MAX;
            if (r2 < gama)
                State_t[infect_prior[j]] = 2;//Recover
            else
                infect_current.push_back(infect_prior[j]);
        }
        infect_prior = infect_current;
        infect_current.clear();
        State.push_back(State_t);
    }
    return State;
}
//Count the number of infected and recovered individuals on the last day
int Get_Infect_num(deque<deque<int>> State, int N) {
    int infect_num = 0;
    for (int j = 0; j < N; j++) {
        if (State[State.size() - 1][j] > 0)
            infect_num++;
    }
    return infect_num;
}
// Randomly extract diagnostic records based on three fixed probabilities
deque<deque<int>>Get_Observ(deque< deque<int> > Neighbor, deque<deque<int>>State, int N, int T_observ, double Test_S, double Test_I, double Test_R) {
    deque<deque<int>>State_in_observ;
    for (int i = 0; i <= T_observ; i++) {
        deque<int>iter;
        for (int j = 0; j < State[i].size(); j++) {
            iter.push_back(State[i][j]);
        }
        State_in_observ.push_back(iter);
        iter.clear();
    }
    deque<deque<int>>diagno;
    int S_num = 0, I_num = 0, R_num = 0;
    for (int day = 0; day < State_in_observ.size(); day++) {
        for (int node = 0; node < N; node++) {
            double v = (double)rand() / RAND_MAX;
            if (State[day][node] == 0 && v < Test_S) {
                deque<int>diag_detail(3, 0);
                diag_detail[0] = day;
                diag_detail[1] = node;
                diag_detail[2] = 0;
                diagno.push_back(diag_detail);
                S_num++;
                diag_detail.clear();
            }
            if (State[day][node] == 1 && v < Test_I) {
                deque<int>diag_detail(3, 0);
                diag_detail[0] = day;
                diag_detail[1] = node;
                diag_detail[2] = 1;
                diagno.push_back(diag_detail);
                I_num++;
                diag_detail.clear();
            }
            if (State[day][node] == 2 && v < Test_R) {
                deque<int>diag_detail(3, 0);
                diag_detail[0] = day;
                diag_detail[1] = node;
                diag_detail[2] = 2;
                diagno.push_back(diag_detail);
                R_num++;
                diag_detail.clear();
            }
        }
    }
    return diagno;
}

int main() {
    double S_u = 1;//up bound of susceptible
    double S_l = 0.95;//lower bound of susceptible
    int D = 7;// Average infectious period
    double gama = 1.0 / D;//Recovers probability
    srand((unsigned)time(NULL));
    int N = 0;//Nodes number
    int M = 0;//Edges number
    FILE* F;
    deque< deque<int> > Neighbor;//Neighbor list
    deque< deque<int> > Link;
    //Get real net
    Get_Real_Network(&N, &M, &Link);
    Real_Network_Transition(N, &M, &Neighbor, &Link);
    printf("N=%d  M=%d\n", N, M);

    ///////////initialization
    int T_end = 7;
    deque<double>Beta_simulation;
    double sim_aver_beta = 0;
    Get_beta_set(&Beta_simulation, N,&sim_aver_beta);
    
    deque<int>seed_set;
    int n = 45;//Seed number
    seed_set = Get_seed(n, Neighbor, N);
    cout << "seed_set.size=" << seed_set.size() << endl;

    ///////////Simulation
    deque< deque<int> > State = Simulation(Neighbor, N, Beta_simulation, gama, T_end, seed_set);
    cout << "State.size=" << State.size() << endl;

    for (int i = 0; i < State.size(); i++) {
        F = fopen("State.txt", "a");
        for (int j = 0; j < N; j++) {
            fprintf(F, "%d ", State[i][j]);
        }
        fprintf(F, "\n");
        fclose(F);
    }
    int Infect_num = Get_Infect_num(State, N);
    cout << "Infect_num=" << Infect_num << endl;
    double Test_S = 0.08 / D;//Daily testing probability of Susceptible
    double Test_I = 0.8 / D;//Daily testing probability of Infected
    double Test_R = 0.1 / D;//Daily testing probability of  Recovered
    deque<deque<int>>diagno = Get_Observ(Neighbor, State, N, T_end - 1, Test_S, Test_I, Test_R);
    ofstream diagno_out("diagno.txt");
    for (int patient_num = 0; patient_num < diagno.size(); patient_num++) {
        diagno_out << diagno[patient_num][0] << ' ' << diagno[patient_num][1] << ' ' << diagno[patient_num][2] << '\t' << endl;
    }
    diagno_out.close();
    return 0;
}




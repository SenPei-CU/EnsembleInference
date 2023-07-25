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
double myRand(double minValue, double maxValue)
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
    double message[100000];
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
//Get diagno info (date & node & real state)
void Get_diagno(deque<deque<int>>* diagno_list) {
    int i = 0;
    int diagno_num;
    ifstream fin("diagno.txt");
    int message[200000];
    while (fin >> message[i]) {
        i++;
    }
    if (i % 3 == 0) {
        diagno_num = i / 3;
        cout << "diagno_num= " << diagno_num << endl;
    }
    else {
        cout << "Diagno.txt has error\n" << endl;
        exit(100);
    }
    for (int j = 0; j < diagno_num; j++) {
        deque<int> iter;
        (*diagno_list).push_back(iter);
        (*diagno_list)[j].push_back(message[3 * j]);
        (*diagno_list)[j].push_back(message[3 * j + 1]);
        (*diagno_list)[j].push_back(message[3 * j + 2]);
    }
    return;
}
// Dynamic Message passing functions
deque<deque<double>> Get_PI(int N, deque<deque<int>> Neighbor, deque<deque<int>> diagno_list, deque<double> beta_set, double gamma, int time, int seed_candi) {

    deque<deque<deque<double>>> theta;//N*deg*(T_end+1)
    deque<deque<deque<double>>> phi;//N*deg*(T_end+1)
    deque<deque<deque<double>>> P_S_ij;//N*deg*(T_end+1)
    deque<double> PS0(N, 1 - 1.0 / seed_candi);
    deque<double> PI0(N, 1.0 / seed_candi);
    deque<double> PR0(N, 0);
    deque<deque<int>> check_S;//N*t 1
    //Initialization

    for (int i = 0; i < N; i++) {
        deque<deque<double>> iter1;
        theta.push_back(iter1);
        deque<deque<double>> iter2;
        phi.push_back(iter2);
        deque<deque<double>> iter3;
        P_S_ij.push_back(iter3);
        for (int j = 0; j < Neighbor[i].size(); j++) {
            deque<double> iter4;
            theta[i].push_back(iter4);
            deque<double> iter5;
            phi[i].push_back(iter5);
            deque<double> iter6;
            P_S_ij[i].push_back(iter6);
            //t=0
            theta[i][j].push_back(1.0);
            // phi[i][j].push_back(1.0 / seed_candi - myRand(0.0005, 0.001));
           //phi[i][j].push_back(rand()/RAND_MAX);
            double v = 1.0 / seed_candi - myRand(0.1 * 1.0 / seed_candi, 0.3 * 1.0 / seed_candi);
            phi[i][j].push_back(v);// if j is S,phi=0
            P_S_ij[i][j].push_back(1.0);
            //t>0
            for (int t = 0; t < time; t++) {
                theta[i][j].push_back(1.0);
                phi[i][j].push_back(1.0);
                P_S_ij[i][j].push_back(1.0);
            }
        }
        deque<int> iter7;
        check_S.push_back(iter7);
        for (int j = 0; j < time + 1; j++) {
            check_S[i].push_back(0);
        }
    }

    int node = 0;
    for (int i = 0; i < diagno_list.size(); i++) {
        if (diagno_list[i][2] == 0) {
            for (int j = 0; j < diagno_list[i][0] + 1; j++) {
                check_S[diagno_list[i][1]][j] = 1;

            }
            node = diagno_list[i][1];
            for (int k = 0; k < Neighbor[node].size(); k++) {
                for (int l = 0; l < Neighbor[Neighbor[node][k]].size(); l++) {
                    if (Neighbor[Neighbor[node][k]][l] == node) {
                        phi[Neighbor[node][k]][l][0] = 0.0;
                    }
                }
            }
            PS0[diagno_list[i][1]] = 1;
            PI0[diagno_list[i][1]] = 0;
        }
    }

    int t = 0;
    double multiply;
    do {
        //equation S7  ----theta(t+1)
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < Neighbor[j].size(); i++) {
                if (check_S[j][t + 1] == 1 || check_S[Neighbor[j][i]][t + 1] == 1) {
                    theta[j][i][t + 1] = 1;
                }
                else {
                    theta[j][i][t + 1] = theta[j][i][t] - beta_set[j] * phi[j][i][t];
                }
            }
        }
        //equation S6  ----P_S_ij(t+1)
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < Neighbor[j].size(); i++) {
                if (check_S[Neighbor[j][i]][t + 1] == 1) {
                    P_S_ij[j][i][t + 1] = 1;
                }
                else {
                    multiply = 1.0;
                    for (int k = 0; k < Neighbor[Neighbor[j][i]].size(); k++) {
                        if (Neighbor[Neighbor[j][i]][k] != j) {
                            multiply *= theta[Neighbor[j][i]][k][t + 1];
                        }
                    }
                    P_S_ij[j][i][t + 1] = PS0[Neighbor[j][i]] * multiply;
                }
            }
        }
        //equation S8  ----phi(t+1)
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < Neighbor[j].size(); i++) {
                if (check_S[Neighbor[j][i]][t + 1] == 1) {
                    phi[j][i][t + 1] = 0;
                }
                else {
                    phi[j][i][t + 1] = (1 - beta_set[j]) * (1 - gamma) * phi[j][i][t] - P_S_ij[j][i][t + 1] + P_S_ij[j][i][t];
                }
            }
        }
        t++;
    } while (t < time);

    //Marginal Probability
    deque<deque<double>>P_S;
    deque<deque<double>>P_I;
    deque<deque<double>>P_R;
    for (int t = 0; t < time + 1; t++) {
        deque<double>iter1(N, 1 - 1.0 / seed_candi);
        P_S.push_back(iter1);
        deque<double>iter2(N, 1.0 / seed_candi);
        P_I.push_back(iter2);
        deque<double>iter3(N, 0);
        P_R.push_back(iter3);
    }
    for (int i = 0; i < diagno_list.size(); i++) {
        if (diagno_list[i][2] == 0) {
            P_S[0][diagno_list[i][1]] = 1;
            P_I[0][diagno_list[i][1]] = 0;
        }
    }

    for (int t = 0; t < time; t++) {
        for (int i = 0; i < N; i++) {
            if (check_S[i][t + 1] == 1) {
                P_S[t + 1][i] = 1;
            }
            else {
                multiply = 1.0;
                for (int k = 0; k < Neighbor[i].size(); k++) {
                    multiply *= theta[i][k][t + 1];
                }
                P_S[t + 1][i] = P_S[0][i] * multiply;//equation S11
            }
            P_R[t + 1][i] = P_R[t][i] + gamma * P_I[t][i];//equation S12
            P_I[t + 1][i] = 1 - P_S[t + 1][i] - P_R[t + 1][i];//equation S13
        }
    }
    theta.clear();
    phi.clear();
    P_S_ij.clear();
    PS0.clear(); PI0.clear(); PR0.clear();
    check_S.clear(); P_R.clear(); P_S.clear();
    return P_I;
}

int main() {
    int D = 7;// Average infectious period
    double gama = 1.0 / D;//Recovers probability
    double beta_u = 0.075;//up bound of beta in samples
    double beta_l = 0.04;//lower bound of beta in samples
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
    int T_start = 0;
    int T_end = 7;
    deque<double>Beta_simulation;
    double sim_aver_beta = 0;
    Get_beta_set(&Beta_simulation, N, &sim_aver_beta);
    
    deque<deque<int>>diagno;
    Get_diagno(&diagno);
    cout << "diagno.size=" << diagno.size() << '\n' << endl;
    cout << "Start at day " << T_start << " and extrapolate the probability of disease on day " << T_end << endl;
    
    ///////////DMP uniform beta
    int seed_candi = 0;
    for (int i = 0; i < diagno.size(); i++) {
        if (diagno[i][2] > 0) {
            seed_candi++;
        }
    }
    deque<double>beta_set_DMP_u;
    for (int i = 0; i < N; i++)
        beta_set_DMP_u.push_back(myRand(beta_l, beta_u));
    deque<deque<double>> P_I_uni = Get_PI(N, Neighbor, diagno, beta_set_DMP_u, gama, T_end, seed_candi);
    F = fopen("DMP_uniform.txt", "a");
    for (int i = 0; i < N; i++) {
        fprintf(F, "%f    ", P_I_uni[T_end][i]);
    }
    fclose(F);

    ////////DMP average beta
    deque<double>beta_set_DMP_a(N, sim_aver_beta);
    deque<deque<double>> P_I_ave = Get_PI(N, Neighbor, diagno, beta_set_DMP_a, gama,T_end, seed_candi);
    F = fopen("DMP_aver.txt", "a");
    for (int i = 0; i < N; i++) {
        fprintf(F, "%f    ", P_I_ave[T_end][i]);
    }
    fclose(F);

    return 0;
}



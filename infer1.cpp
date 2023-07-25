
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
void Get_beta_set(deque<double>* beta_set, int N) {
    int i = 0;
    ifstream fin("beta_hum874_b.txt");//vector 1*N
    double message[10000];
    while (fin >> message[i]) {
        i++;
    }
    if (i != N) {
        cout << "beta.txt has error\n" << endl;
        exit(100);
    }
    for (int j = 0; j < i; j++) {
        (*beta_set).push_back(message[j]);
    }
    return;
}
//Get diagno info (date & node & real state)
void Get_diagno(deque<deque<int>>* diagno_list) {
    int i = 0;
    int diagno_num;
    ifstream fin("diagno_1.txt");
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
//Update Model for one step (Master equation)
deque<double> Model_Integ(deque<double> state, deque<deque<int>> Neighbor, int N, double gama, deque<double>Beta_sam) {
    deque<double> state_new = state;
    for (int i = 0; i < N; i++) {
        double temp_S = 1.0;
        for (int j = 0; j < Neighbor[i].size(); j++) {
            temp_S *= 1 - Beta_sam[i] * state[Neighbor[i][j] + N];
        }
        state_new[i] += -state[i] * (1 - temp_S);
        state_new[i + N] += state[i] * (1 - temp_S) - gama * state[i + N];
        state_new[i + 2 * N] += gama * state[i + N];
        
    }

    for (int k = 0; k < N; k++) {
        double sum = state_new[k] + state_new[k + N] + state_new[k + 2 * N];
        state_new[k] = state_new[k] / sum;
        state_new[k + N] = state_new[k + N] / sum;
        state_new[k + 2 * N] = state_new[k + 2 * N] / sum;
    }
    return state_new;
}
// Generate samples
deque<deque<double>> Get_Sample(int sample_num, double S_u, double S_l, int N) {
    deque<deque<double>> Sample;
    for (int i = 0; i < sample_num; i++) {
        deque <double> Sample_i(3 * N);
        Sample.push_back(Sample_i);
    }

    for (int i = 0; i < sample_num; i++) {//The first N position in each sample is the susceptibility probability, the second N position is the infection probability, and the third N position is the recovery probability
        for (int j = 0; j < N; j++) {
            double v1 = MyRand(S_l, S_u);
            Sample[i][j] = v1;
            Sample[i][j + N] = 1 - v1;
            Sample[i][j + 2 * N] = 0;
        }
    }
    return Sample;
}

// To get unique diagnostic time and location index
deque<deque<int>>Get_repeat_time(deque<int>itime) {
    deque<deque<int>>repeat;
    deque<int>non_repeat_day;
    deque<int>repeat_position;
    for (int t = 0; t < itime.size() - 1; t++) {
        if (itime[t] != itime[t + 1]) {
            repeat_position.push_back(t);
            non_repeat_day.push_back(itime[t]);
        }
    }
    non_repeat_day.push_back(itime[itime.size() - 1]);
    repeat_position.push_back(itime.size()-1);
    repeat.push_back(non_repeat_day);//The first row is the time, the second row is the position
    repeat.push_back(repeat_position);
    non_repeat_day.clear();
    repeat_position.clear();
    return repeat;
}

//Capture the currently available diagnostic record
deque<deque<int>>Get_diagno_current(deque<deque<int>> diagno, int curr_time, int N) {
    deque<deque<int>>diagno_current;
    deque<int>flag(N, 0);
    int length = 0;
    for (int a = 0; a < diagno.size(); a++) {
        if (diagno[a][0] >= curr_time && flag[diagno[a][1]] == 0) {
            flag[diagno[a][1]] = 1;
            diagno_current.push_back(diagno[a]);
        }
    }
    flag.clear();
    return diagno_current;
}

//////////The prior probabilities are obtained by backward inference
deque<double>Get_Bayesisan(deque<deque<int>> diagno, int diag_num, deque<deque<int>> Neighbor, deque<double>Sample_t_sam, int infer_time, int N, double gama, deque<double>Beta_sam) {
    deque<double>P_SIR(3, 0);
    int nodenum = diagno[diag_num][1];
    int State = diagno[diag_num][2];
    deque<int>diagno_new_time;
    diagno_new_time.push_back(infer_time);
    for (int z = 0; z < diagno.size(); z++)
        diagno_new_time.push_back(diagno[z][0]);
    deque<deque<int>>repeat_current = Get_repeat_time(diagno_new_time);
    double lh_S = 1, lh_I = 1;//Likelihood
    deque<double>Sample_t_sam_S = Sample_t_sam;//1*3N
    deque<double>Sample_t_sam_I = Sample_t_sam;
    Sample_t_sam_S[nodenum] = 1;
    Sample_t_sam_S[nodenum + N] = 0;
    Sample_t_sam_S[nodenum + 2 * N] = 0;
    Sample_t_sam_I[nodenum] = 0;
    Sample_t_sam_I[nodenum + N] = 1;
    Sample_t_sam_I[nodenum + 2 * N] = 0;
    if (State == 1) {
        //The first row is the time, the second row is the position
        for (int k = 0; k < repeat_current[0].size() - 1; k++) {
            for (int t = repeat_current[0][k]; t < repeat_current[0][k + 1]; t++) {
                Sample_t_sam_S = Model_Integ(Sample_t_sam_S, Neighbor, N, gama, Beta_sam);
                Sample_t_sam_I = Model_Integ(Sample_t_sam_I, Neighbor, N, gama, Beta_sam);
            }
            for (int ii = repeat_current[1][k]; ii < repeat_current[1][k + 1]; ii++) {
                lh_S *= Sample_t_sam_S[diagno[ii][1] + diagno[ii][2] * N];
                lh_I *= Sample_t_sam_I[diagno[ii][1] + diagno[ii][2] * N];

                if (diagno[ii][2] == 0) {
                    Sample_t_sam_S[diagno[ii][1]] = 1;
                    Sample_t_sam_S[diagno[ii][1] + N] = 0;
                    Sample_t_sam_S[diagno[ii][1] + 2 * N] = 0;
                    Sample_t_sam_I[diagno[ii][1]] = 1;
                    Sample_t_sam_I[diagno[ii][1] + N] = 0;
                    Sample_t_sam_I[diagno[ii][1] + 2 * N] = 0;
                }
                if (diagno[ii][2] == 1) {
                    Sample_t_sam_S[diagno[ii][1]] = 0;
                    Sample_t_sam_S[diagno[ii][1] + N] = 1;
                    Sample_t_sam_S[diagno[ii][1] + 2 * N] = 0;
                    Sample_t_sam_I[diagno[ii][1]] = 0;
                    Sample_t_sam_I[diagno[ii][1] + N] = 1;
                    Sample_t_sam_I[diagno[ii][1] + 2 * N] = 0;
                }
                if (diagno[ii][2] == 2) {
                    Sample_t_sam_S[diagno[ii][1]] = 0;
                    Sample_t_sam_S[diagno[ii][1] + N] = 0;
                    Sample_t_sam_S[diagno[ii][1] + 2 * N] = 1;
                    Sample_t_sam_I[diagno[ii][1]] = 0;
                    Sample_t_sam_I[diagno[ii][1] + N] = 0;
                    Sample_t_sam_I[diagno[ii][1] + 2 * N] = 1;
                }
            }
        }
        double sum = Sample_t_sam[nodenum] * lh_S + Sample_t_sam[nodenum + N] * lh_I;
        P_SIR[0] = Sample_t_sam[nodenum] * lh_S / sum;
        P_SIR[1] = Sample_t_sam[nodenum + N] * lh_I / sum;
        P_SIR[2] = 0;
    }
    if (State == 2) {
        double lh_R = 1;
        deque<double>Sample_t_sam_R = Sample_t_sam;
        Sample_t_sam_R[nodenum] = 0;
        Sample_t_sam_R[nodenum + N] = 0;
        Sample_t_sam_R[nodenum + 2 * N] = 1;
        for (int k = 0; k < repeat_current[0].size() - 1; k++) {
            for (int t = repeat_current[0][k]; t < repeat_current[0][k + 1]; t++) {
                Sample_t_sam_S = Model_Integ(Sample_t_sam_S, Neighbor, N, gama, Beta_sam);
                Sample_t_sam_I = Model_Integ(Sample_t_sam_I, Neighbor, N, gama, Beta_sam);
                Sample_t_sam_R = Model_Integ(Sample_t_sam_R, Neighbor, N, gama, Beta_sam);
            }
            for (int ii = repeat_current[1][k]; ii < repeat_current[1][k + 1]; ii++) {
                lh_S *= Sample_t_sam_S[diagno[ii][1] + diagno[ii][2] * N];
                lh_I *= Sample_t_sam_I[diagno[ii][1] + diagno[ii][2] * N];
                lh_R *= Sample_t_sam_R[diagno[ii][1] + diagno[ii][2] * N];

                if (diagno[ii][2] == 0) {
                    Sample_t_sam_S[diagno[ii][1]] = 1;
                    Sample_t_sam_S[diagno[ii][1] + N] = 0;
                    Sample_t_sam_S[diagno[ii][1] + 2 * N] = 0;
                    Sample_t_sam_I[diagno[ii][1]] = 1;
                    Sample_t_sam_I[diagno[ii][1] + N] = 0;
                    Sample_t_sam_I[diagno[ii][1] + 2 * N] = 0;
                    Sample_t_sam_R[diagno[ii][1]] = 1;
                    Sample_t_sam_R[diagno[ii][1] + N] = 0;
                    Sample_t_sam_R[diagno[ii][1] + 2 * N] = 0;
                }
                if (diagno[ii][2] == 1) {
                    Sample_t_sam_S[diagno[ii][1]] = 0;
                    Sample_t_sam_S[diagno[ii][1] + N] = 1;
                    Sample_t_sam_S[diagno[ii][1] + 2 * N] = 0;
                    Sample_t_sam_I[diagno[ii][1]] = 0;
                    Sample_t_sam_I[diagno[ii][1] + N] = 1;
                    Sample_t_sam_I[diagno[ii][1] + 2 * N] = 0;
                    Sample_t_sam_R[diagno[ii][1]] = 0;
                    Sample_t_sam_R[diagno[ii][1] + N] = 1;
                    Sample_t_sam_R[diagno[ii][1] + 2 * N] = 0;
                }
                if (diagno[ii][2] == 2) {
                    Sample_t_sam_S[diagno[ii][1]] = 0;
                    Sample_t_sam_S[diagno[ii][1] + N] = 0;
                    Sample_t_sam_S[diagno[ii][1] + 2 * N] = 1;
                    Sample_t_sam_I[diagno[ii][1]] = 0;
                    Sample_t_sam_I[diagno[ii][1] + N] = 0;
                    Sample_t_sam_I[diagno[ii][1] + 2 * N] = 1;
                    Sample_t_sam_R[diagno[ii][1]] = 0;
                    Sample_t_sam_R[diagno[ii][1] + N] = 0;
                    Sample_t_sam_R[diagno[ii][1] + 2 * N] = 1;
                }
            }
        }
        double sum = Sample_t_sam[nodenum] * lh_S + Sample_t_sam[nodenum + N] * lh_I + Sample_t_sam[nodenum + 2 * N] * lh_R;
        P_SIR[0] = Sample_t_sam[nodenum] * lh_S / sum;
        P_SIR[1] = Sample_t_sam[nodenum + N] * lh_I / sum;
        P_SIR[2] = Sample_t_sam[nodenum + 2 * N] * lh_R / sum;
    }
    return P_SIR;
}
//// Backward temporal propagation
deque<deque<deque<double>>>Backward(deque<deque<int>> Neighbor, deque<deque<int>> diagno, deque<deque<double>> Sample_t, int sample_num, int N, int curr_time, double gama, deque<deque<double>>Sample_beta) {
    //Intialize
    deque<deque<deque<double>>>P_post;//100*dia_num_after_t*3
    for (int i = 0; i < sample_num; i++) {
        deque<deque<double>>iter_diagno_list;
        for (int j = 0; j < diagno.size(); j++) {
            deque<double>iter_patient(3, 0);
            iter_diagno_list.push_back(iter_patient);
        }
        P_post.push_back(iter_diagno_list);
    }

    for (int i = 0; i < sample_num; i++) {
        for (int diag_num = 0; diag_num < diagno.size(); diag_num++) {
            if (diagno[diag_num][2] != 0)
                P_post[i][diag_num] = Get_Bayesisan(diagno, diag_num, Neighbor, Sample_t[i], curr_time, N, gama, Sample_beta[i]);
            else {
                P_post[i][diag_num][0] = 1;
                P_post[i][diag_num][1] = 0;
                P_post[i][diag_num][2] = 0;
            }
        }
    }
    return P_post;
}
//Update SIR probability of diagnostic individuals
deque<deque<double>> Update_Sample(deque<deque<deque<double>>>post_Inference, deque<deque<double>>Sample, int sample_num, deque<deque<int>>diagno, int N) {
    deque<deque<double>>Sample_new = Sample;
    for (int i = 0; i < sample_num; i++) {
        for (int j = 0; j < diagno.size(); j++) {
            Sample_new[i][diagno[j][1] + N] = post_Inference[i][j][1];
            Sample_new[i][diagno[j][1]] = post_Inference[i][j][0];
            Sample_new[i][diagno[j][1] + 2 * N] = post_Inference[i][j][2];
        }
    }
    return Sample_new;
}

//Cocorrelation adjustment
double average(deque<double>temp) {//get mean
    double sum = 0;
    for (int i = 0; i < temp.size(); i++)
        sum += temp[i];
    double average = sum / temp.size();
    return average;
}
//get the covariance, variance
double variance(deque<double>temp_1, deque<double>temp_2) {
    double average_1 = average(temp_1);
    double average_2 = average(temp_2);
    double sum = 0;
    for (int i = 0; i < temp_1.size(); i++) {
        sum += (temp_1[i] - average_1) * (temp_2[i] - average_2);
    }
    return sum / temp_1.size();
}
// Cross-ensemble covariability 
deque<deque<double>> Cross_ensemble(deque<deque<int>> diagno_current, deque<deque<int>> Neighbor, int sample_num, int N, deque<deque<double>>Sample_t, deque<deque<deque <double>>>post_Inference, int t) {
    deque<deque<double>>Sample_t_new = Sample_t;
    for (int j = 0; j < diagno_current.size(); j++) {
        deque<double> I_j_post; deque<double> S_j_post; deque<double> R_j_post;
        deque<double> I_j_prior; deque<double> S_j_prior; deque<double> R_j_prior;
        deque<double> I_k_prior; deque<double> S_k_prior; deque<double> R_k_prior;
        for (int k = 0; k < Neighbor[diagno_current[j][1]].size(); k++) {
            for (int i = 0; i < sample_num; i++) {
                I_j_prior.push_back(Sample_t[i][diagno_current[j][1] + N]);
                S_j_prior.push_back(Sample_t[i][diagno_current[j][1]]);
                R_j_prior.push_back(Sample_t[i][diagno_current[j][1] + 2 * N]);
                I_j_post.push_back(post_Inference[i][j][1]);
                S_j_post.push_back(post_Inference[i][j][0]);
                R_j_post.push_back(post_Inference[i][j][2]);
                S_k_prior.push_back(Sample_t[i][Neighbor[diagno_current[j][1]][k]]);
                I_k_prior.push_back(Sample_t[i][Neighbor[diagno_current[j][1]][k] + N]);
                R_k_prior.push_back(Sample_t[i][Neighbor[diagno_current[j][1]][k] + 2 * N]);

            }
            double weight_S, weight_I, weight_R;
            if (isnan(variance(I_j_prior, I_k_prior) / variance(I_j_prior, I_j_prior)) == 1)
                weight_I = MyRand(0, 0.1);
            else
                weight_I = variance(I_j_prior, I_k_prior) / variance(I_j_prior, I_j_prior);
            if (isnan(variance(S_j_prior, S_k_prior) / variance(S_j_prior, S_j_prior)) == 1)
                weight_S = MyRand(0, 0.1);
            else
                weight_S = variance(S_j_prior, S_k_prior) / variance(S_j_prior, S_j_prior);
            if (isnan(variance(R_j_prior, R_k_prior) / variance(R_j_prior, R_j_prior)) == 1)
                weight_R = MyRand(0, 0.1);
            else
                weight_R = variance(R_j_prior, R_k_prior) / variance(R_j_prior, R_j_prior);
            deque<int>flag_S(sample_num, 0); deque<int>flag_R(sample_num, 0); deque<int>flag_I(sample_num, 0);
            double I_k_post = 0, S_k_post = 0, R_k_post = 0;
            int I_effect_num = 0, S_effect_num = 0, R_effect_num = 0;
            for (int b = 0; b < sample_num; b++) {
                Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + N] += weight_I * (I_j_post[b] - I_j_prior[b]);
                Sample_t_new[b][Neighbor[diagno_current[j][1]][k]] += weight_S * (S_j_post[b] - S_j_prior[b]);
                Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + 2 * N] += weight_R * (R_j_post[b] - R_j_prior[b]);

                if (isnan(Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + N]) == 1)
                    flag_I[b] = 1;
                else {
                    if (Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + N] > 1)
                        flag_I[b] = 2;
                    else {
                        if (Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + N] < 0)
                            flag_I[b] = 3;
                        else {
                            I_k_post += Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + N];
                            I_effect_num++;
                        }


                    }

                }
                if (isnan(Sample_t_new[b][Neighbor[diagno_current[j][1]][k]]) == 1)
                    flag_S[b] = 1;
                else {
                    if (Sample_t_new[b][Neighbor[diagno_current[j][1]][k]] > 1)
                        flag_S[b] = 2;
                    else {
                        if (Sample_t_new[b][Neighbor[diagno_current[j][1]][k]] < 0)
                            flag_S[b] = 3;
                        else {
                            S_k_post += Sample_t_new[b][Neighbor[diagno_current[j][1]][k]];
                            S_effect_num++;
                        }
                    }

                }

                if (isnan(Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + 2 * N]) == 1)
                    flag_R[b] = 1;
                else {
                    if (Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + 2 * N] > 1)
                        flag_R[b] = 2;
                    else {
                        if (Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + 2 * N] < 0)
                            flag_R[b] = 3;
                        else {
                            R_k_post += Sample_t_new[b][Neighbor[diagno_current[j][1]][k] + 2 * N];
                            R_effect_num++;
                        }
                    }

                }
            }
            //normalize
            for (int g = 0; g < sample_num; g++) {
                if (flag_I[g] == 1) {
                    Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + N] = S_k_post / S_effect_num;
                }
                if (flag_S[g] == 1) {
                    Sample_t_new[g][Neighbor[diagno_current[j][1]][k]] = I_k_post / I_effect_num;
                }
                if (flag_R[g] == 1) {
                    Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + 2 * N] = R_k_post / R_effect_num;
                }
                if (flag_I[g] == 2) {
                    Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + N] = 1 - MyRand(0, 0.1);
                }
                if (flag_S[g] == 2) {
                    Sample_t_new[g][Neighbor[diagno_current[j][1]][k]] = 1 - MyRand(0, 0.1);
                }

                if (flag_R[g] == 2) {
                    Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + 2 * N] = 1 - MyRand(0, 0.1);
                }
                if (flag_I[g] == 3) {
                    Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + N] = MyRand(0, 0.1);
                }
                if (flag_S[g] == 3) {
                    Sample_t_new[g][Neighbor[diagno_current[j][1]][k]] = MyRand(0, 0.1);
                }
                if (flag_R[g] == 3) {
                    Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + 2 * N] = MyRand(0, 0.1);
                }
            }
            for (int g = 0; g < sample_num; g++) {
                double sum = Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + N] + Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + 2 * N] + Sample_t_new[g][Neighbor[diagno_current[j][1]][k]];
                Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + N] = Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + N] / sum;
                Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + 2 * N] = Sample_t_new[g][Neighbor[diagno_current[j][1]][k] + 2 * N] / sum;
                Sample_t_new[g][Neighbor[diagno_current[j][1]][k]] = Sample_t_new[g][Neighbor[diagno_current[j][1]][k]] / sum;
            }

            flag_S.clear();
            flag_I.clear();
            flag_R.clear();
            I_k_prior.clear();
            S_k_prior.clear();
            R_k_prior.clear();
            I_j_prior.clear(); S_j_prior.clear(); R_j_prior.clear();
            I_j_post.clear(); S_j_post.clear(); R_j_post.clear();

        }
    }
    //Update the probability of the patients on the diagnostic record
    Sample_t_new = Update_Sample(post_Inference, Sample_t_new, sample_num, diagno_current, N);
    return  Sample_t_new;
}

//Forward propagation, and the sample probability obtained by propagation is used as the prior probability of the next day
deque<deque<double>>Integration(deque<deque<double>>Sample_t, deque<deque<int>> Neighbor, int N, int sample_num, double gama, deque<deque<double>>Sample_beta) {
    deque<deque<double>>next_day_State;
    for (int i = 0; i < sample_num; i++) {
        deque<double>iter6 = Model_Integ(Sample_t[i], Neighbor, N, gama, Sample_beta[i]);
        next_day_State.push_back(iter6);
        iter6.clear();
    }
    return    next_day_State;
}
int main() {
    int sample_num = 100;//number of sample
    double S_u = 1;//up bound of susceptible
    double S_l = 0.95;//lower bound of susceptible
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
    Get_beta_set(&Beta_simulation, N);

    deque<deque<int>>diagno;
    Get_diagno(&diagno);
    cout << "diagno.size=" << diagno.size() << '\n' << endl;
    cout << "Start at day " << T_start << " and extrapolate the probability of disease on day " << T_end << endl;
 
 ///////////inference
    deque<deque<double>>Sample_beta;//Initialize transmission rate in sample
    for (int i = 0; i < sample_num; i++) {
        deque<double>iter_beta(N, 0);
        Sample_beta.push_back(iter_beta);
    }
    for (int i = 0; i < sample_num; i++) {
        for (int j = 0; j < N; j++)
            Sample_beta[i][j] = MyRand(beta_l, beta_u);
    }

    deque<deque<deque<double>>>Sample;//Store the sample changes corresponding to different time t
    deque<deque<double>>Sample_0 = Get_Sample(sample_num, S_u, S_l, N);//100*(3N)
    Sample.push_back(Sample_0);

    for (int i = T_start; i < T_end; i++) {
        deque<deque<double>>iter_time;
        for (int j = 0; j < sample_num; j++) {
            deque<double>iter_samp(3 * N);
            iter_time.push_back(iter_samp);
        }
        Sample.push_back(iter_time);
    }


    for (int t = T_start; t < T_end; t++) {
        deque<deque<int>>diagno_current = Get_diagno_current(diagno, t, N);
        deque<int>current_time_diagno_set;
        for (int patient = 0; patient < diagno_current.size(); patient++) {
            if (diagno_current[patient][0] == t) {
                current_time_diagno_set.push_back(patient);
            }
        }
        if (current_time_diagno_set.size() > 0) {
            for (int i = 0; i < sample_num; i++) {
                for (int ii = 0; ii < current_time_diagno_set.size(); ii++) {
                    if (diagno_current[current_time_diagno_set[ii]][2] == 1) {
                        Sample[t][i][diagno_current[current_time_diagno_set[ii]][1]] = 0;
                        Sample[t][i][diagno_current[current_time_diagno_set[ii]][1] + N] = 1;
                        Sample[t][i][diagno_current[current_time_diagno_set[ii]][1] + 2 * N] = 0;
                    }
                    if (diagno_current[current_time_diagno_set[ii]][2] == 0) {
                        Sample[t][i][diagno_current[current_time_diagno_set[ii]][1]] = 1;
                        Sample[t][i][diagno_current[current_time_diagno_set[ii]][1] + N] = 0;
                        Sample[t][i][diagno_current[current_time_diagno_set[ii]][1] + 2 * N] = 0;
                    }
                    if (diagno_current[current_time_diagno_set[ii]][2] == 2) {
                        Sample[t][i][diagno_current[current_time_diagno_set[ii]][1]] = 0;
                        Sample[t][i][diagno_current[current_time_diagno_set[ii]][1] + N] = 0;
                        Sample[t][i][diagno_current[current_time_diagno_set[ii]][1] + 2 * N] = 1;
                    }
                }
            }
            diagno_current.erase(diagno_current.begin(), diagno_current.begin() + (int)current_time_diagno_set.size());
        }


        deque<deque<deque<double>>>post_Inference = Backward(Neighbor, diagno_current, Sample[t], sample_num, N, t, gama, Sample_beta);//Inference

        Sample[t] = Cross_ensemble(diagno_current, Neighbor, sample_num, N, Sample[t], post_Inference, t);//Uptade_Neighbor

        Sample[t + 1] = Integration(Sample[t], Neighbor, N, sample_num, gama, Sample_beta);//Integration

        post_Inference.clear();
        diagno_current.clear();
        current_time_diagno_set.clear();
    }

    deque<double>inference_average(N);
    F = fopen("infer.txt", "a");
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < sample_num; j++) {
            sum += Sample[T_end][j][i + N];
        }
        inference_average[i] = sum / sample_num;
        fprintf(F, "%f ", inference_average[i]);
    }
    fclose(F);
    return 0;
}



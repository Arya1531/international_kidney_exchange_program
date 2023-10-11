/*
*    int_kidney_exchange
*    contribution+benefit.cpp
*    Purpose: computational study for Computing Balanced Solutions for Large International Kidney
*			  Exchange Schemes When Cycle Length Is Unbounded
*             using the benefit and contribution value as initial allocations with equal country sizes
*
*    @author Xin Ye
*    @version 1.0 09/10/2023
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version and the Gurobi License.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*    GNU General Public License for more details.
*
*/

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/matching.h>
#include <lemon/adaptors.h>
#include <lemon/core.h>
#include <lemon/base.cc>
#include <lemon/concepts/maps.h>
#include <time.h>
#include <iomanip>
#include <glpk.h>
#include <math.h>
#include <stdio.h>
//#include "windows.h"
//#include "psapi.h"
#include <iterator>
#include "gurobi_c++.h"


using namespace lemon;
using namespace std;

double cpuTime();
bool is_next_char_digit(string& line, unsigned int l);
unsigned int char2uint(char& p);
void undi_lemon(unsigned int& m, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, vector<unsigned short int>& label_positions, ListGraph& g, ListDigraph& g_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, ListGraph::EdgeMap<double>& edge_card_weight, ListDigraph::ArcMap<unsigned short int>& arc_card_weight, unsigned short int& no_of_nodes);
void coop_game(ListGraph& g, vector<double>& v, vector<double>& v_impu, vector<double>& v_S, unsigned int& S, vector<unsigned short int>& s, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, ListGraph::EdgeMap<double>& edge_card_weight, bool& dispy, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, map<int, int>& numofMaxSolution, unsigned short int& Q, bool& arbitray_maximum, vector<pair<int, int>>& cycle_distri, map<int, int>& cycle_dis, double& game_generation, std::map<int, std::map<int, int>>& cycle_dis_period);
void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, unsigned short int& k, ListGraph& g, ListDigraph& g_original, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, unsigned int& m, unsigned short int& no_of_nodes);

void min_d_1(vector<unsigned short int>& node_arrives, ListGraph& g, ListDigraph& g_original, vector<pair<int, int>>& arc_pair, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<double>& v, bool& target_omega, vector<double>& target, vector<double>& credit, ListGraph::EdgeMap<double>& edge_card_weight, double& t0, vector<vector<unsigned short int>>& actual_alloc, vector<double>& v_impu, vector<int>& nodeset, vector<pair<int, int>>& cycle_distri, vector<double>& d, double& M_total, double& d_total, bool& c_involved, map<int, int>& cycle_dis, map<int, int>& numofMaxSolution, bool& arbitray_maximum, unsigned short int& initialSize, vector<double>& v_S, unsigned int& S, double& core_100, long& negative_core, double& d_c_total, unsigned short int inst, bool lex_min, double& max_d, double& game_generation, double& solution_concept_time, double& scenario_time, std::map<int, std::map<int, int>>& cycle_dis_period, std::map<int, std::map<int, int>>& cycle_dis_arbitrary_period);
void arbitraryMaximum(vector<unsigned short int>& node_arrives, ListGraph& g, ListDigraph& g_original, vector<pair<int, int>>& arc_pair, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<double>& v, bool& target_omega, vector<double>& target, vector<double>& credit, ListGraph::EdgeMap<double>& edge_card_weight, double& t0, vector<vector<unsigned short int>>& actual_alloc, vector<double>& v_impu, vector<int>& nodeset, vector<pair<int, int>>& cycle_distri, vector<double>& d, double& M_total, double& d_total, bool& c_involved, map<int, int>& cycle_dis, map<int, int>& numofMaxSolution, bool& arbitray_maximum, unsigned short int& initialSize, vector<double>& v_S, unsigned int& S, double& core_100, long& negative_core, unsigned short int inst, double& max_d, double& game_generation, double& solution_concept_time, std::map<int, std::map<int, int>>& cycle_dis_arbitrary_period);
void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, unsigned short int& initialSize);
void period_0(unsigned short int& Q, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, vector<unsigned short int>& s, unsigned short int& Vp, vector<unsigned short int>& node_arrives, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, vector<double>& credit, unsigned short int& initialSize, vector<bool>& leaving);
void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods);
void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, vector<unsigned short int>& s, vector<double>& d, vector<double>& target);
void ILP_d1_gurobi(unsigned short int& Q, unsigned short int& N, ListDigraph& g_original, unsigned short int& Vp, vector<unsigned short int>& node_arrives, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<pair<int, int>>& arc_pair, vector<int>& nodeset, vector<vector<unsigned short int>>& actual_alloc, double& M, double& M_total, vector<unsigned short int>& s, vector<pair<int, int>>& cycle_distri, vector<double>& target, vector<bool>& leaving, vector<double>& d, double& d_total, bool& c_involved, vector<double>& credit, map<int, int>& cycle_dis, bool lex_min, unsigned short int inst, std::map<int, std::map<int, int>>& cycle_dis_period);
void pair_arcs(unsigned short int& Q, ListDigraph& g_original, vector<unsigned short int>& node_arrives, ListDigraph::NodeMap<bool>& active_nodes_original, vector<pair<int, int>>& arc_pair, vector<int>& nodeset);
void lex_min_d_star(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, vector<double>& credit, vector<GRBVar>& var_lexmin, vector<GRBModel>& vector_model, unsigned short int& track);
void lex_min_n_star(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, vector<GRBVar>& var_lexmin, vector<GRBModel>& vector_model, unsigned short int& track);
void sort_d_t(vector<double>& d_t, vector<GRBVar>& var_bi, long& col_num, unsigned short int& N, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, unsigned short int& t, vector<double>& credit, double& epsilon, vector<GRBVar>& var_lexmin, vector<unsigned short int>& N_star);
void lexmin_searching(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, vector<double>& credit, vector<GRBVar>& var_lexmin, unsigned short int inst, vector<GRBModel>& vector_model, unsigned short int& track);
void cycle_distribution(std::map<int, std::map<int, int>>& cycle_dis_period, map<int, int>& cycle_dis, vector<pair<int, int>>& cycle_distri, unsigned short int& N, unsigned short int& Q);
void epsilon_func(vector<double> target, vector<double> credit, double& epsilon, unsigned short int N);
double frac(double ori);
int main() {
	try {
		std::cout << "I solemnly swear that I am up to no good." << endl;
		bool target_omega = true; // true: benefit value, false: constribution value
		bool dispy = false; // true: information in terminal while running
		bool disp = false; // true: extremely detailed information while running, avoid with large graphs
		bool c_involved = false;// true: credits considred; false:without credits 
		bool arbitray_maximum = false; //true: arbitray maximum cycple packing
		bool lex_min = false;
		string solution_concept;
		string version;
		bool d1 = true;
		bool d_c = true;
		bool lexmin_call = true;
		bool lexmin_c_call = true;
		bool arbitrary = true;
		if (d1) {
			version = "d1";
		}
		else {
			if (lexmin_call) {
				version = "lexmin_call";
			}
			if (lexmin_c_call) {
				version = "lexmin_c_call";
			}
			if (d_c) {
				version = "d_c";
			}
			if (arbitrary) {
				version = "arbitrary";
			}
		}

		if (target_omega) {
			solution_concept = "benefit";
			cout << solution_concept << endl;
		}
		else {
			solution_concept = "contribution";
			cout << solution_concept << endl;
		}
		unsigned short int years = 6;
		unsigned short int periods_per_year = 4;
		// input parameters and data
		unsigned short int N;
		unsigned short int inst; // instance number, integer between 0 and 99
		map<int, int> cycle_dis;
		map<int, int> cycle_dis_d;
		map<int, int> cycle_dis_t_c;
		map<int, int> cycle_dis_arbitrary;
		map<int, int> cycle_dis_lexmin;
		map<int, int> cycle_dis_lexmin_c;
		map<int, map<int, int>> cycle_dis_d_period;
		map<int, map<int, int>> cycle_dis_t_c_period;
		map<int, map<int, int>> cycle_dis_arbitrary_period;
		map<int, map<int, int>> cycle_dis_lexmin_period;
		map<int, map<int, int>> cycle_dis_lexmin_c_period;
		map<int, int> numofMaxSolution;
		map<int, int> numofMaxSolution_t_d;
		map<int, int> numofMaxSolution_arbitrary;
		vector<double> relative_d1_N(12, 0);
		vector<double> relative_d1_N_c(12, 0);
		vector<double> relative_arbitrary_N(12, 0);
		vector<double> relative_lexmin(12, 0);
		vector<double> relative_lexmin_c(12, 0);
		vector<double> max_d1_N(12, 0);
		vector<double> max_d1_N_c(12, 0);
		vector<double> max_arbitrary_N(12, 0);
		vector<double> max_lexmin(12, 0);
		vector<double> max_lexmin_c(12, 0);
		vector<double> M_N(12, 0);
		vector<double> M_N_d_c(12, 0);
		vector<double> M_N_d_arbitrary(12, 0);
		vector<double> M_N_lex_min(12, 0);
		vector<double> M_N_lex_min_c(12, 0);
		vector<double> core_dis_N(12, 0);
		vector<double> core_dis_N_c(12, 0);
		vector<double> core_dis_N_arbitrary(12, 0);
		vector<long> out_of_core(12, 0);
		vector<long> out_of_core_c(12, 0);
		vector<long> out_of_core_arbitrary(12, 0);
		vector<double> data_preparation_N(12, 0);
		vector<double> graph_building_N(12, 0);
		vector<double> game_generation_arbitrary_N(12, 0);
		vector<double> game_generation_d1_N(12, 0);
		vector<double> game_generation_d1_c_N(12, 0);
		vector<double> game_generation_lexmin_N(12, 0);
		vector<double> game_generation_lexmin_c_N(12, 0);
		vector<double> time_arbitrary_N(12, 0);
		vector<double> time_d1_N(12, 0);
		vector<double> time_d1_c_N(12, 0);
		vector<double> time_lex_min_N(12, 0);
		vector<double> time_lex_min_c_N(12, 0);
		vector<double> total_time_arbitrary_N(12, 0);
		vector<double> total_time_d1_N(12, 0);
		vector<double> total_time_d1_c_N(12, 0);
		vector<double> total_time_lex_min_N(12, 0);
		vector<double> total_time_lex_min_c_N(12, 0);
		vector<double> solution_concept_time_arbitrary_N(12, 0);
		vector<double> solution_concept_time_d1_N(12, 0);
		vector<double> solution_concept_time_d1_c_N(12, 0);
		vector<double> solution_concept_time_lexmin_N(12, 0);
		vector<double> solution_concept_time_lexmin_c_N(12, 0);
		for (N = 4; N < 11; ++N) {
			cycle_dis.clear();
			cycle_dis_d.clear();
			cycle_dis_t_c.clear();
			cycle_dis_arbitrary.clear();
			cycle_dis_lexmin.clear();
			cycle_dis_lexmin_c.clear();
			double relative_d1 = 0;
			double relative_d1_c = 0;
			double relative_d1_arbitrary = 0;
			double relative_lexmin_0 = 0;
			double relative_lexmin_c_0 = 0;
			double data_preparation = 0;
			double graph_building = 0;
			double game_generation_arbitrary = 0;
			double game_generation_d1 = 0;
			double game_generation_d1_c = 0;
			double game_generation_lexmin = 0;
			double game_generation_lexmin_c = 0;
			double time_arbitrary = 0;
			double time_d1 = 0;
			double time_d1_c = 0;
			double time_lex_min = 0;
			double time_lex_min_c = 0;
			double total_time_arbitrary = 0;
			double total_time_d1 = 0;
			double total_time_d1_c = 0;
			double total_time_lex_min = 0;
			double total_time_lex_min_c = 0;
			double solution_concept_time_arbitrary = 0;
			double solution_concept_time_d1 = 0;
			double solution_concept_time_d1_c = 0;
			double solution_concept_time_lexmin = 0;
			double solution_concept_time_lexmin_c = 0;
			double max_d1 = 0;
			double max_d1_c = 0;
			double max_d1_arbitrary = 0;
			double max_lexmin_0 = 0;
			double max_lexmin_c_0 = 0;
			double M_100 = 0;
			double M_100_d_c = 0;
			double M_100_d_arbitrary = 0;
			double M_lex_min = 0;
			double M_lex_min_c = 0;
			double core_100 = 0;
			double core_d = 0;
			double core_d_c = 0;
			double core_d_arbitrary = 0;
			long negative_core = 0;
			long negative_core_d = 0;
			long negative_core_d_c = 0;
			long negative_core_d_arbitrary = 0;

			for (inst = 0; inst < 100; ++inst) {
				std::cout << N << "countries" << " " << "instance_" << inst << endl;
				string line;
				ifstream inp;
				unsigned short int graph_size = 2000;
				// read the data
				inp.open("genxml-" + to_string(inst) + ".xml"); // 1 out of the 100 instances generated by William Pettersson's web tool: https://wpettersson.github.io/kidney-webapp/#/
				getline(inp, line);
				inp.close();

				
				unsigned short int Vp = 4 * (unsigned short int)((graph_size / 4) / N);
				unsigned short int no_of_nodes = N * Vp;
				vector<unsigned int> arc_out(0, 0);
				vector<unsigned int> arc_in(0, 0);
				unsigned int m = 0;
				unsigned short int k = 0;
				unsigned short int M = 0;
				double M_total = 0;
				vector<unsigned short int> node_labels(no_of_nodes, 0);
				vector<unsigned short int> label_positions(graph_size, graph_size + 1);
				unsigned int S = pow(2, N) - 2;
				// biparite graph
				ListGraph g;
				// orginal compatibility graph
				ListDigraph g_original; 
				vector<ListGraph::Node> c(no_of_nodes);
				vector<ListGraph::Node> c_b(no_of_nodes); 
				vector<ListDigraph::Node> c_original(no_of_nodes);
				double t0 = cpuTime();
				// paste the data
				xml_parser(line, node_labels, label_positions, c, c_b, c_original, k, g, g_original, arc_in, arc_out, m, no_of_nodes);
				double t1 = cpuTime();
				data_preparation += t1 - t0;
				
				unsigned short int periods = years * periods_per_year;
				unsigned short int initialSize = Vp / 4;
				vector<unsigned short int> no_of_active_nodes(N, initialSize); //initial active nodes for the period 0
				ListGraph::NodeMap<bool> active_nodes(g);
				ListDigraph::NodeMap<bool> active_nodes_original(g_original);
				for (unsigned short int i = 0; i < N; i++) {
					for (unsigned short int j = 0; j < Vp; j++) {
						active_nodes[c[i * Vp + j]] = false;
						active_nodes[c_b[i * Vp + j]] = false;
						active_nodes_original[c_original[i * Vp + j]] = false;
					}
				}
				//read the seed
				string line_seed;
				ifstream seed_doc;
				std::cout << "start reading the seed doc" << endl;
				seed_doc.open("seeds/seeds/n" + to_string(N) + "inst" + to_string(inst) + ".txt");
				getline(seed_doc, line_seed);
				seed_doc.close();
				unsigned int seed = 0;
				seed = stoi(line_seed); 
				srand(seed);

				// determine starting pairs and arrival times of others
				initial_pairs(Vp, N, active_nodes, active_nodes_original, c, c_b, c_original, initialSize);
				vector<unsigned short int> node_arrives(no_of_nodes, 0);
				arrival_times(node_arrives, Vp, N, active_nodes, c, periods);
				t1 = cpuTime();
				double rand_time = t1 - t0;


				t0 = cpuTime();
				ListGraph::EdgeMap<double> edge_card_weight(g, 0);
				ListDigraph::ArcMap<unsigned short int> arc_card_weight(g_original, 0);
				// build the graph
				undi_lemon(m, arc_in, arc_out, label_positions, g, g_original, c, c_b, c_original, edge_card_weight, arc_card_weight, no_of_nodes);
				t1 = cpuTime();
				graph_building += t1 - t0;

				vector<double> v_impu(N, 0);
				vector<double> v(N + 1, 0);
				vector<double> v_S(S, 0);
				vector<unsigned short int> s(N, 0);
				vector<vector<unsigned short int>> actual_alloc_d1C(periods, vector<unsigned short int>(N, 0));
				vector<vector<unsigned short int>> actual_alloc_d1(periods, vector<unsigned short int>(N, 0));
				vector<vector<unsigned short int>> actual_alloc_rand(periods, vector<unsigned short int>(N, 0));
				double prec = pow(10, -7);
				vector<double> target(N, 0);
				vector<vector<double>> init_alloc_d1C(periods, vector<double>(N, 0));
				vector<vector<double>> init_alloc_d1(periods, vector<double>(N, 0));
				vector<vector<double>> init_alloc_rand(periods, vector<double>(N, 0));
				vector<double> credit(N, 0);
				vector<bool> leaving(no_of_nodes, false);
				unsigned short int Q = 0;


				//set the active nodes for the period 0
				period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit, initialSize, leaving);

				vector<pair<int, int>> arc_pair;
				vector<vector<unsigned short int>> actual_alloc;
				vector<int> nodeset(2000, 0);
				vector<pair<int, int>> cycle_distri;
				vector<double> d(N, 0);
				double d_total = 0;
				double d_c_total = 0;
				double max_d = 0;

				//------------arbitray maximum cycle packing----------------
				if (arbitrary) {
					std::cout << N << "countries" << " " << "instance_" << inst << "starts arbitrary" << endl;
					arbitray_maximum = true;
					period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit, initialSize, leaving);
					t0 = cpuTime();
					arbitraryMaximum(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, target_omega, target, credit, edge_card_weight, t0, actual_alloc, v_impu, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_arbitrary, numofMaxSolution, arbitray_maximum, initialSize, v_S, S, core_100, negative_core, inst, max_d, game_generation_arbitrary, solution_concept_time_arbitrary, cycle_dis_arbitrary_period);
					t1 = cpuTime();
					total_time_arbitrary += t1 - t0;
					arbitray_maximum = false;
					M_100_d_arbitrary += M_total;
					//core_d_arbitrary += core_100;
					//negative_core_d_arbitrary += negative_core;
					relative_d1_arbitrary += (d_total / M_total);
					max_d1_arbitrary += max_d;
					std::cout << N << "countries" << " " << "instance_" << inst << "arbitrary done...";
				}
				//---------d1----------
				if (d1) {
					std::cout << N << "countries" << " " << "instance_" << inst << "starts d1" << endl;
					period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit, initialSize, leaving);
					t0 = cpuTime();
					min_d_1(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, target_omega, target, credit, edge_card_weight, t0, actual_alloc, v_impu, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_d, numofMaxSolution, arbitray_maximum, initialSize, v_S, S, core_100, negative_core, d_c_total, inst, lex_min, max_d, game_generation_d1, solution_concept_time_d1, time_d1, cycle_dis_d_period, cycle_dis_arbitrary_period);
					t1 = cpuTime();
					total_time_d1 += t1 - t0;
					cout << "total_time_d1: " << total_time_d1 << endl;
					relative_d1 += (d_total / M_total);
					max_d1 += max_d;
					std::cout << "relative_d1: " << relative_d1 << endl;
					std::cout << "the number of countries: " << N << " " << "relative_d1" << " " << inst << " " << relative_d1 / (inst + 1) << endl;
					M_100 += M_total;
					std::cout << "the number of countries: " << N << " " << "relative_d1" << " " << inst << " " << M_100 / (inst + 1);
					std::cout << N << "countries" << " " << "instance_" << inst << "d1 done...";
				}
				
				// -----------d1+c----------------
				if (d_c) {
					std::cout << N << "countries" << " " << "instance_" << inst << "starts d1+c" << endl;
					c_involved = true;
					period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit, initialSize, leaving);
					t0 = cpuTime();
					min_d_1(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, target_omega, target, credit, edge_card_weight, t0, actual_alloc, v_impu, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_t_c, numofMaxSolution, arbitray_maximum, initialSize, v_S, S, core_100, negative_core, d_c_total, inst, lex_min, max_d, game_generation_d1_c, solution_concept_time_d1_c, time_d1_c, cycle_dis_t_c_period, cycle_dis_arbitrary_period);
					t1 = cpuTime();
					total_time_d1_c += t1 - t0;
					c_involved = false;
					M_100_d_c += M_total;
					relative_d1_c += (d_total / M_total);
					max_d1_c += max_d;
					std::cout << N << "countries" << " " << "instance_" << inst << "d1+c done...";
				}
				
				//-----------------lexmin-----------------
				if (lexmin_call) {
					std::cout << N << "countries" << " " << "instance_" << inst << "lexmin starts...";
					lex_min = true;
					period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit, initialSize, leaving);
					t0 = cpuTime();
					min_d_1(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, target_omega, target, credit, edge_card_weight, t0, actual_alloc, v_impu, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_lexmin, numofMaxSolution, arbitray_maximum, initialSize, v_S, S, core_100, negative_core, d_c_total, inst, lex_min, max_d, game_generation_lexmin, solution_concept_time_lexmin, time_lex_min, cycle_dis_lexmin_period, cycle_dis_arbitrary_period);
					t1 = cpuTime();
					total_time_lex_min += t1 - t0;
					relative_lexmin_0 += (d_total / M_total);
					max_lexmin_0 += max_d;
					lex_min = false;
					M_lex_min += M_total;
					std::cout << N << "countries" << " " << "instance_" << inst << "relative deviation" << d_total / M_total << "lexmin done...";
				}
				//-----------------lexmin+c----------------
				if (lexmin_c_call) {
					std::cout << N << "countries" << " " << "instance_" << inst << "lexmin+c starts...";
					lex_min = true;
					c_involved = true;
					period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit, initialSize, leaving);
					t0 = cpuTime();
					min_d_1(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, target_omega, target, credit, edge_card_weight, t0, actual_alloc, v_impu, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_lexmin_c, numofMaxSolution, arbitray_maximum, initialSize, v_S, S, core_100, negative_core, d_c_total, inst, lex_min, max_d, game_generation_lexmin_c, solution_concept_time_lexmin_c, time_lex_min_c, cycle_dis_lexmin_c_period, cycle_dis_arbitrary_period);
					t1 = cpuTime();
					total_time_lex_min_c += t1 - t0;
					relative_lexmin_c_0 += (d_total / M_total);
					max_lexmin_c_0 += max_d;
					lex_min = false;
					c_involved = false;
					M_lex_min_c += M_total;
					std::cout << N << "countries" << " " << "instance_" << inst << "relative deviation" << d_total / M_total << "lexmin+c done...";
				}
				
			}
			relative_d1_N[N - 4] = relative_d1 / 100;
			relative_d1_N_c[N - 4] = relative_d1_c / 100;
			relative_arbitrary_N[N - 4] = relative_d1_arbitrary / 100;
			relative_lexmin[N - 4] = relative_lexmin_0 / 100;
			relative_lexmin_c[N - 4] = relative_lexmin_c_0 / 100;
			max_d1_N[N - 4] = max_d1 / 100;
			max_d1_N_c[N - 4] = max_d1_c / 100;
			max_arbitrary_N[N - 4] = max_d1_arbitrary / 100;
			max_lexmin[N - 4] = max_lexmin_0 / 100;
			max_lexmin_c[N - 4] = max_lexmin_c_0 / 100;
			M_N[N - 4] = M_100 / 100;
			M_N_d_c[N - 4] = M_100_d_c / 100;
			M_N_d_arbitrary[N - 4] = M_100_d_arbitrary / 100;
			M_N_lex_min[N - 4] = M_lex_min / 100;
			M_N_lex_min_c[N - 4] = M_lex_min_c / 100;
			data_preparation_N[N - 4] = data_preparation / 100;
			graph_building_N[N - 4] = graph_building / 100;
			game_generation_arbitrary_N[N - 4] = game_generation_arbitrary / 100;
			game_generation_d1_N[N - 4] = game_generation_d1 / 100;
			game_generation_d1_c_N[N - 4] = game_generation_d1_c / 100;
			game_generation_lexmin_N[N - 4] = game_generation_lexmin / 100;
			game_generation_lexmin_c_N[N - 4] = game_generation_lexmin_c / 100;
			time_arbitrary_N[N - 4] = time_arbitrary / 100;
			time_d1_N[N - 4] = time_d1 / 100;
			time_d1_c_N[N - 4] = time_d1_c / 100;
			time_lex_min_N[N - 4] = time_lex_min / 100;
			time_lex_min_c_N[N - 4] = time_lex_min_c / 100;
			total_time_arbitrary_N[N - 4] = total_time_arbitrary / 100;
			total_time_d1_N[N - 4] = total_time_d1 / 100;
			total_time_d1_c_N[N - 4] = total_time_d1_c / 100;
			total_time_lex_min_N[N - 4] = total_time_lex_min / 100;
			total_time_lex_min_c_N[N - 4] = total_time_lex_min_c / 100;
			solution_concept_time_arbitrary_N[N - 4] = solution_concept_time_arbitrary / 100;
			solution_concept_time_d1_N[N - 4] = solution_concept_time_d1 / 100;
			solution_concept_time_d1_c_N[N - 4] = solution_concept_time_d1_c / 100;
			solution_concept_time_lexmin_N[N - 4] = solution_concept_time_lexmin / 100;
			solution_concept_time_lexmin_c_N[N - 4] = solution_concept_time_lexmin_c / 100;
			ofstream res;
			res.open(solution_concept+"/"+version+"/results.txt", ofstream::out | ofstream::trunc);
			for (unsigned short int i = 0; i < N - 3; i++) {
				res << i + 4 << "countries" << endl;
				res << "data preparation: " << data_preparation_N[i] << endl;
				res << "build graph: " << graph_building_N[i] << endl;
				res << "minimizing d_1: " << relative_d1_N[i] << endl;
				res << "minimizing max_d_1: " << max_d1_N[i] << endl;
				res << "average number of transplants: " << M_N[i] << endl;
				res << "total time: " << total_time_d1_N[i] << endl;
				res << "scenario time: " << time_d1_N[i] << endl;
				res << "game generation: " << game_generation_d1_N[i] << endl;
				res << "solution concept: " << solution_concept_time_d1_N[i] << endl;
				res << "minimizing d_1_c: " << relative_d1_N_c[i] << endl;
				res << "minimizing max_d_1_c: " << max_d1_N_c[i] << endl;
				res << "average number of transplants_c: " << M_N_d_c[i] << endl;
				res << "total time: " << total_time_d1_c_N[i] << endl;
				res << "scenario time: " << time_d1_c_N[i] << endl;
				res << "game generation: " << game_generation_d1_c_N[i] << endl;
				res << "solution concept: " << solution_concept_time_d1_c_N[i] << endl;
				res << "minimizing d_1_arbitrary: " << relative_arbitrary_N[i] << endl;
				res << "minimizing max_d_1_arbitrary: " << max_arbitrary_N[i] << endl;
				res << "average number of transplants_arbitrary: " << M_N_d_arbitrary[i] << endl;
				res << "total time: " << total_time_arbitrary_N[i] << endl;
				res << "scenario time: " << time_arbitrary_N[i] << endl;
				res << "game generation: " << game_generation_arbitrary_N[i] << endl;
				res << "solution concept: " << solution_concept_time_arbitrary_N[i] << endl;
				res << "lex min " << relative_lexmin[i] << endl;
				res << "lex min max d " << max_lexmin[i] << endl;
				res << "average number of transplants: " << M_N_lex_min[i] << endl;
				res << "total time: " << total_time_lex_min_N[i] << endl;
				res << "scenario time: " << time_lex_min_N[i] << endl;
				res << "game generation: " << game_generation_lexmin_N[i] << endl;
				res << "solution concept: " << solution_concept_time_lexmin_N[i] << endl;
				res << "lex min+c " << relative_lexmin_c[i] << endl;
				res << "lex min+c max " << max_lexmin_c[i] << endl;
				res << "average number of transplants: " << M_N_lex_min_c[i] << endl;
				res << "total time: " << total_time_lex_min_c_N[i] << endl;
				res << "scenario time: " << time_lex_min_c_N[i] << endl;
				res << "game generation: " << game_generation_lexmin_c_N[i] << endl;
				res << "solution concept: " << solution_concept_time_lexmin_c_N[i] << endl;
				res << endl;
			}

			res.close();

			// cycle distributions
			vector<long> value(5, 0);
			if (d1) {
				ofstream res_dis;
				res_dis.open(solution_concept+"/"+version+"/cycle_dis_d" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_d) {
					res_dis << elem.first << ": " << elem.second << endl;
					value[0] += elem.first * elem.second;
				}
				res_dis << endl;
				res_dis.close();
			}
			

			if (d_c) {
				ofstream res_dis_c;
				res_dis_c.open(solution_concept + "/" + version + "/cycle_dis_c" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_t_c) {
					res_dis_c << elem.first << ": " << elem.second << endl;
					value[1] += elem.first * elem.second;
				}
				res_dis_c << endl;
				res_dis_c.close();
			}
			
			if (arbitrary) {
				ofstream res_dis_arbitrary;
				res_dis_arbitrary.open(solution_concept + "/" + version + "/cycle_dis_arbitrary" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_arbitrary) {
					res_dis_arbitrary << elem.first << ": " << elem.second << endl;
					value[2]= elem.first * elem.second;
				}
				res_dis_arbitrary << endl;
				res_dis_arbitrary.close();
			}
			
			if (lexmin_call) {
				ofstream res_dis_lexmin;
				res_dis_lexmin.open(solution_concept + "/" + version + "/cycle_dis_lexmin" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_lexmin) {
					res_dis_lexmin << elem.first << ": " << elem.second << endl;
					value[3]+= elem.first * elem.second;
				}
				res_dis_lexmin << endl;
				res_dis_lexmin.close();
			}
			
			if (lexmin_c_call) {
				ofstream res_dis_lexmin_c;
				long value_lexmin_c = 0;
				res_dis_lexmin_c.open(solution_concept + "/" + version + "/cycle_dis_lexmin_c" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_lexmin_c) {
					res_dis_lexmin_c << elem.first << ": " << elem.second << endl;
					value[4]+= elem.first * elem.second;
				}
				res_dis_lexmin_c << endl;
				res_dis_lexmin_c.close();
			}
			

			// cycle distributions seperating by periods
			vector<long> check(5, 0);
			if (d1) {
				ofstream res_cycle_dis_d_period;
				for (unsigned short int i = 0; i < 24; ++i) {
					res_cycle_dis_d_period.open(solution_concept + "/" + version + "/cycle_dis_d_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
					for (const auto& elem : cycle_dis_d_period[(N) * (i + 1)]) {
						check[0] += elem.first * elem.second;
						res_cycle_dis_d_period << elem.first << ": " << elem.second << endl;
					}
					res_cycle_dis_d_period << endl;
					res_cycle_dis_d_period.close();
				}
				res_cycle_dis_d_period.close();
			}
			

			if (d_c) {
				ofstream res_cycle_dis_c_period;
				for (unsigned short int i = 0; i < 24; ++i) {
					res_cycle_dis_c_period.open(solution_concept + "/" + version + "/cycle_dis_c_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
					for (const auto& elem : cycle_dis_t_c_period[(N) * (i + 1)]) {
						check[1] += elem.first * elem.second;
						res_cycle_dis_c_period << elem.first << ": " << elem.second << endl;
					}
					res_cycle_dis_c_period << endl;
					res_cycle_dis_c_period.close();
				}
				res_cycle_dis_c_period.close();
			}
			
			if (arbitrary) {
				ofstream res_cycle_dis_arbitrary_period;
				for (unsigned short int i = 0; i < 24; ++i) {
					res_cycle_dis_arbitrary_period.open(solution_concept + "/" + version + "/cycle_dis_arbitrary_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
					for (const auto& elem : cycle_dis_arbitrary_period[(N) * (i + 1)]) {
						check[2] += elem.first * elem.second;
						res_cycle_dis_arbitrary_period << elem.first << ": " << elem.second << endl;
					}
					res_cycle_dis_arbitrary_period << endl;
					res_cycle_dis_arbitrary_period.close();
				}
				res_cycle_dis_arbitrary_period.close();
			}
			

			if (lexmin_call) {
				ofstream res_dis_lexmin_period;
				for (unsigned short int i = 0; i < 24; ++i) {
					res_dis_lexmin_period.open(solution_concept + "/" + version + "/lexmin/cycle_dis_lexmin_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
					for (const auto& elem : cycle_dis_lexmin_period[(N) * (i + 1)]) {
						check[3] += elem.first * elem.second;
						res_dis_lexmin_period << elem.first << ": " << elem.second << endl;
					}
					res_dis_lexmin_period << endl;
					res_dis_lexmin_period.close();
				}
				res_dis_lexmin_period.close();
			}
			
			if (lexmin_c_call) {
				ofstream res_cycle_dis_lexmin_c_period;
				for (unsigned short int i = 0; i < 24; ++i) {
					res_cycle_dis_lexmin_c_period.open(solution_concept + "/" + version + "/cycle_dis_lexmin_c_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
					for (const auto& elem : cycle_dis_lexmin_c_period[(N) * (i + 1)]) {
						check[4] += elem.first * elem.second;
						res_cycle_dis_lexmin_c_period << elem.first << ": " << elem.second << endl;
					}
					res_cycle_dis_lexmin_c_period << endl;
					res_cycle_dis_lexmin_c_period.close();
				}
				res_cycle_dis_lexmin_c_period.close();
			}
			
			for (unsigned short int i = 0; i < 5; ++i) {
				if (value[i] != check[i]) {
					cout << "Error in the number of transplants" << endl;
				}
			}
			


		}

	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << endl;
		std::cout << e.getMessage() << endl;
	}
	catch (...) {
		std::cout << "Exception during optimization" << endl;
	}
	return 0;
}




void coop_game(ListGraph& g, vector<double>& v, vector<double>& v_impu, vector<double>& v_S, unsigned int& S, vector<unsigned short int>& s, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, ListGraph::EdgeMap<double>& edge_card_weight, bool& dispy, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, map<int, int>& numofMaxSolution, unsigned short int& Q, bool& arbitray_maximum, vector<pair<int, int>>& cycle_distri, map<int, int>& cycle_dis, double& game_generation, std::map<int, std::map<int, int>>& cycle_dis_arbitrary_period) {
	vector<bool> a(N, false);
	double t0 = cpuTime();
	for (unsigned int i = 0; i < N; i++) {
		ListGraph::NodeMap<bool> coal1(g, false);
		//ListDigraph::NodeMap<bool> coal1_original(g_original, false);
		for (unsigned short int k = i * Vp; k < (i + 1) * Vp; k++) {
			if (active_nodes[c[k]]) {
				coal1[c[k]] = true;
				coal1[c_b[k]] = true;
			}
		}

		MaxWeightedPerfectMatching<FilterNodes<ListGraph>, ListGraph::EdgeMap<double>> coal_m1(FilterNodes<ListGraph>(g, coal1), edge_card_weight);
		coal_m1.run();
		v_impu[i] = coal_m1.matchingWeight();
		std::cout << "finish generating v_impu" << endl;
		ListGraph::NodeMap<bool> coal2(g, false);
		for (unsigned int j = 0; j < N; j++) {
			if (i != j) {
				for (unsigned short int k = j * Vp; k < (j + 1) * Vp; k++) {
					if (active_nodes[c[k]]) {
						coal2[c[k]] = true;
						coal2[c_b[k]] = true;
					}
				}
			}
		}


		//MaxWeightedPerfectMatching<FilterNodes<ListGraph>> coal_m2(g, coal2);
		MaxWeightedPerfectMatching<FilterNodes<ListGraph>, ListGraph::EdgeMap<double>> coal_m2(FilterNodes<ListGraph>(g, coal2), edge_card_weight);
		coal_m2.run();
		v[i] = coal_m2.matchingWeight();
		std::cout << "finish generating v_i" << endl;

	}


	std::cout << "finish generating the copy" << endl;
	FilterNodes<ListGraph> sg(g, active_nodes);
	MaxWeightedPerfectMatching<FilterNodes<ListGraph>, ListGraph::EdgeMap<double>> grand_coal(FilterNodes<ListGraph>(g, active_nodes), edge_card_weight);
	grand_coal.run();
	v[N] = grand_coal.matchingWeight();
	double t1 = cpuTime();
	std::cout << "Time in game generation " << t1 - t0 << endl;
	game_generation += t1 - t0;
	grand_coal.matchingMap();
	std::cout << "finish generating v[N]" << endl;
	if (arbitray_maximum) {
		unsigned short int a = 0;
		unsigned short int b = 0;
		for (FilterNodes<ListGraph>::NodeIt n(sg); n != INVALID; ++n) {
			cout << "n(sg): " << sg.id(n) << endl;
			if (!(grand_coal.matching(n) == INVALID) && edge_card_weight[grand_coal.matching(n)] > 0 && sg.id(n) % 2 == 0) {
				cycle_distri.push_back({ sg.id(n) / 2, (sg.id(grand_coal.mate(n)) - 1) / 2 });
				leaving[sg.id(n) / 2] = true;
				cout << "sg.id(n):" << sg.id(n) << endl;
				if (a == v[N]) {
					std::cout << sg.id(n) << endl;
				}
				cout << "sg.id(n): " << sg.id(n) << endl;
				for (unsigned short int i = 0; i < N; ++i) {
					if (i * Vp <= sg.id(n) / 2 && sg.id(n) / 2 < (i + 1) * Vp) {
						++s[i];
					}
				}

			}
			//cycle_distri.push_back(arc_pair[i - 1]);
		}
		cout << "finish cycple_distri" << endl;
		cycle_distribution(cycle_dis_arbitrary_period, cycle_dis, cycle_distri, N, Q);
	}



	if (dispy)
		std::cout << "grand coal: " << v[N] << endl;

	if (dispy) {
		std::cout << "s: ";
		for (unsigned short int i = 0; i < N; i++) {
			std::cout << s[i] << " " << endl;
		}
	}
	return;
}

void min_d_1(vector<unsigned short int>& node_arrives, ListGraph& g, ListDigraph& g_original, vector<pair<int, int>>& arc_pair, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<double>& v, bool& target_omega, vector<double>& target, vector<double>& credit, ListGraph::EdgeMap<double>& edge_card_weight, double& t0, vector<vector<unsigned short int>>& actual_alloc, vector<double>& v_impu, vector<int>& nodeset, vector<pair<int, int>>& cycle_distri, vector<double>& d, double& M_total, double& d_total, bool& c_involved, map<int, int>& cycle_dis, map<int, int>& numofMaxSolution, bool& arbitray_maximum, unsigned short int& initialSize, vector<double>& v_S, unsigned int& S, double& core_100, long& negative_core, double& d_c_total, unsigned short int inst, bool lex_min, double& max_d, double& game_generation, double& solution_concept_time, double& scenario_time, std::map<int, std::map<int, int>>& cycle_dis_period, std::map<int, std::map<int, int>>& cycle_dis_arbitrary_period) {
	Q = 0;
	d_total = 0;
	d_c_total = 0;
	M_total = 0;
	core_100 = 0;
	negative_core = 0;
	max_d = 0;
	if (dispy)
		std::cout << " --== Without lex min matching == -- " << endl;
	for (unsigned short int i = 0; i < N; i++) {
		d[i] = 0;
		s[i] = 0;
		credit[i] = 0;
		no_of_active_nodes[i] = initialSize;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (node_arrives[i * Vp + j] == 0) {
				active_nodes[c[i * Vp + j]] = true;
				active_nodes[c_b[i * Vp + j]] = true;
				active_nodes_original[c_original[i * Vp + j]] = true;
			}
			else {
				active_nodes[c[i * Vp + j]] = false;
				active_nodes[c_b[i * Vp + j]] = false;
				active_nodes_original[c_original[i * Vp + j]] = false;
			}
		}
	}
	while (Q < periods) {
		if (dispy) {
			std::cout << "--== PERIOD " << Q + 1 << " ==--" << endl;
		}
		if (dispy) {
			std::cout << "Number of active nodes: ";
			for (unsigned short int i = 0; i < N; i++)
				std::cout << no_of_active_nodes[i] << " ";
			std::cout << endl;
		}
		// cooperative game and target
		std::cout << "start generating values" << endl;
		coop_game(g, v, v_impu, v_S, S, s, c, c_b, edge_card_weight, dispy, Vp, N, active_nodes, leaving, numofMaxSolution, Q, arbitray_maximum, cycle_distri, cycle_dis, game_generation, cycle_dis_arbitrary_period);
		double t0 = cpuTime();
		double suma = 0;
		double sumimpu = 0;
		if (target_omega) {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i] - v_impu[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i] - v_impu[i]) / suma);
		}
		else {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i]) / suma);
		}
		double t1 = cpuTime();
		solution_concept_time += t1 - t0;

		if (suma>0) {
			cout << "denominator>0"<<endl;
		}
		else{
			ofstream res_denominator;
			res_denominator << N << " countries" << " " << "instance_" << inst << " round_" << Q << endl;
			res_denominator.close();
		}
		//init_alloc[Q] = target;
		if (dispy) {
			if (target_omega) {
				std::cout << "Benefit: ";
			}
			else {
				std::cout << "Contribution: ";
			}
			for (unsigned short int i = 0; i < N; i++) {
				std::cout << target[i] << " ";
			}
			std::cout << endl;
		}
		t0 = cpuTime();
		ILP_d1_gurobi(Q, N, g_original, Vp, node_arrives, active_nodes, active_nodes_original, arc_pair, nodeset, actual_alloc, v[N], M_total, s, cycle_distri, target, leaving, d, d_total, c_involved, credit, cycle_dis, lex_min, inst, cycle_dis_period);
		t1 = cpuTime();
		scenario_time += t1 - t0;
		std::string fileName = "output_credits_benefit_no_c.txt";
		std::ofstream outputFile(fileName);

		if (outputFile.is_open()) {
			// Write data to the file
			outputFile << "credits--Period: " << Q << endl;
			for (unsigned short int i = 0; i < N - 4; ++i) {
				outputFile << "credits" << credit[i] << endl;
			}

			// Close the file
			outputFile.close();

			std::cout << "Data has been written to the file." << std::endl;
		}
		else {
			std::cout << "Unable to open the file." << std::endl;
		}
		Q++;
		changing_nodes(active_nodes, active_nodes_original, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, c_b, c_original, s, d, target);
		if (dispy)
			cin.get();
	}
	for (unsigned short int i = 0; i < N; ++i) {
		d_total += abs(d[i]);
		if (c_involved) {
			d_c_total += abs(credit[i]);
		}
	}
	vector<double> max_deviation(N, 0);
	for (unsigned short int i = 0; i < N; ++i) {
		max_deviation[i] = abs(d[i]) / M_total;
	}
	std::sort(max_deviation.begin(), max_deviation.end());
	max_d = max_deviation[N - 1];
	return;
}

void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, vector<unsigned short int>& s, vector<double>& d, vector<double>& target) {
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		target[i] = 0;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (leaving[i * Vp + j]) {
				active_nodes[c[i * Vp + j]] = false;
				active_nodes[c_b[i * Vp + j]] = false;
				active_nodes_original[c_original[i * Vp + j]] = false;
				no_of_active_nodes[i]--;
				leaving[i * Vp + j] = false;
			}
			else {
				if (active_nodes[c[i * Vp + j]] && node_arrives[i * Vp + j] == Q - 4) {
					active_nodes[c[i * Vp + j]] = false;
					active_nodes[c_b[i * Vp + j]] = false;
					active_nodes_original[c_original[i * Vp + j]] = false;
					no_of_active_nodes[i]--;
				}
			}
			if (node_arrives[i * Vp + j] == Q) {
				active_nodes[c[i * Vp + j]] = true;
				active_nodes[c_b[i * Vp + j]] = true;
				active_nodes_original[c_original[i * Vp + j]] = true;
				no_of_active_nodes[i]++;
			}
		}
	}
	return;
}

void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, unsigned short int& initialSize) {
	unsigned short int coal = rand() % Vp;
	unsigned short int count = 0;
	for (unsigned short int i = 0; i < N; i++) {
		while (count < initialSize) {
			if (active_nodes[c[i * Vp + coal]]) {
				coal = rand() % Vp;
			}
			else {
				active_nodes[c[i * Vp + coal]] = true;
				active_nodes[c_b[i * Vp + coal]] = true;
				active_nodes_original[c_original[i * Vp + coal]] = true;
				count++;
				coal = rand() % Vp;
			}
		}
		count = 0;
	}
	return;
}

void period_0(unsigned short int& Q, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, vector<unsigned short int>& s, unsigned short int& Vp, vector<unsigned short int>& node_arrives, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, vector<double>& credit, unsigned short int& initialSize, vector<bool>& leaving) {
	Q = 0;
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		credit[i] = 0;
		no_of_active_nodes[i] = initialSize;
		for (unsigned short int j = 0; j < Vp; j++) {
			leaving[i * Vp + j] = false;
			if (node_arrives[i * Vp + j] == 0) {
				active_nodes[c[i * Vp + j]] = true;
				active_nodes[c_b[i * Vp + j]] = true;
				active_nodes_original[c_original[i * Vp + j]] = true;
			}
			else {
				active_nodes[c[i * Vp + j]] = false;
				active_nodes[c_b[i * Vp + j]] = false;
				active_nodes_original[c_original[i * Vp + j]] = false;
			}
		}
	}
	return;
}

void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods) {
	for (unsigned short int i = 0; i < N; i++) {
		for (unsigned short int j = 0; j < Vp; j++) {
			if (!(active_nodes[c[i * Vp + j]])) {
				node_arrives[i * Vp + j] = rand() % (periods - 1) + 1;
			}
		}
	}
	return;
}



void undi_lemon(unsigned int& m, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, vector<unsigned short int>& label_positions, ListGraph& g, ListDigraph& g_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, ListGraph::EdgeMap<double>& edge_card_weight, ListDigraph::ArcMap<unsigned short int>& arc_card_weight, unsigned short int& no_of_nodes) {
	bool halt = false;
	for (int i = 0; i < no_of_nodes; i++) {
		ListGraph::Edge e = g.addEdge(c[i], c_b[i]);
		edge_card_weight[e] = 0;
	}
	for (unsigned int i = 0; i < m; i++) {
		if (label_positions[arc_in[i]] < no_of_nodes) { //XY: filter 65535 positions
			ListGraph::Edge e = g.addEdge(c[label_positions[arc_out[i]]], c_b[label_positions[arc_in[i]]]);
			edge_card_weight[e] = 1;
			ListDigraph::Arc a_original = g_original.addArc(c_original[label_positions[arc_out[i]]], c_original[label_positions[arc_in[i]]]);
			arc_card_weight[a_original] = 1;
		}
	}
	return;
}

void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, unsigned short int& k, ListGraph& g, ListDigraph& g_original, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, unsigned int& m, unsigned short int& no_of_nodes) {
	unsigned int l = 6;
	unsigned short int n = 0; //XY: track the number of nodes
	while (l < line.size() - 7) {
		if (line[l] == '<' && line[l + 1] == 'e') {
			l = l + 17;
			n++;
			if (!is_next_char_digit(line, l)) {
				node_labels[n - 1] = char2uint(line[l]); //XY: donor id
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					node_labels[n - 1] = 10 * char2uint(line[l]) + char2uint(line[l + 1]);
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						node_labels[n - 1] = 100 * char2uint(line[l]) + 10 * char2uint(line[l + 1]) + char2uint(line[l + 2]);
						l = l + 2;
					}
					else {
						node_labels[n - 1] = 1000 * char2uint(line[l]) + 100 * char2uint(line[l + 1]) + 10 * char2uint(line[l + 2]) + char2uint(line[l + 3]);
						l = l + 3;
					}
				}
			}
			if (n + k - 1 == node_labels[n - 1]) {
				label_positions[n + k - 1] = n - 1;
			}
			else {
				while (n + k - 1 < node_labels[n - 1]) {
					label_positions[n + k - 1] = 65535;
					label_positions.push_back(0);
					k++;
				}
				label_positions[n + k - 1] = n - 1;
			}

			c[n - 1] = g.addNode();//XY: add donor ids
			c_b[n - 1] = g.addNode();//changed by XY, add patient ids
			c_original[n - 1] = g_original.addNode();//changed by XY, add patient-donor pairs to the original graph
			l = l + 9;
			if (!is_next_char_digit(line, l)) {
				////donor_ages.push_back(char2uint(line[l]));
				//donor_ages[n - 1] = char2uint(line[l]);
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					////donor_ages.push_back(10*char2uint(line[l])+char2uint(line[l+1]));
					//donor_ages[n - 1] = 10*char2uint(line[l])+char2uint(line[l+1]);
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						////donor_ages.push_back(100*char2uint(line[l])+10*char2uint(line[l+1])+char2uint(line[l+2]));
						//donor_ages[n - 1] = 100*char2uint(line[l])+10*char2uint(line[l+1])+char2uint(line[l+2]);
						l = l + 2;
					}
					else {
						////if (!is_next_char_digit(line, l + 3)){
						////donor_ages.push_back(1000*char2uint(line[l])+100*char2uint(line[l+1])+10*char2uint(line[l+2])+char2uint(line[l+3]));
						//donor_ages[n - 1] = 1000*char2uint(line[l])+100*char2uint(line[l+1])+10*char2uint(line[l+2])+char2uint(line[l+3]);
						l = l + 3;
						////}
					}
				}
			}
			l = l + 25;
			if (!is_next_char_digit(line, l)) {
				if (node_labels[n - 1] != char2uint(line[l]))
					std::cout << "ID ERROR!" << endl;
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					if (node_labels[n - 1] != 10 * char2uint(line[l]) + char2uint(line[l + 1]))
						std::cout << "ID ERROR!" << endl;
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						if (node_labels[n - 1] != 100 * char2uint(line[l]) + 10 * char2uint(line[l + 1]) + char2uint(line[l + 2]))
							std::cout << "ID ERROR!" << endl;
						l = l + 2;
					}
					else {
						//if (!is_next_char_digit(line, l + 3)){
						if (node_labels[n - 1] != 1000 * char2uint(line[l]) + 100 * char2uint(line[l + 1]) + 10 * char2uint(line[l + 2]) + char2uint(line[l + 3]))
							std::cout << "ID ERROR!" << endl;
						l = l + 3;
						//}
					}
				}
			}
			if (line[l + 21] == 'm')
				l = l + 29;
			else
				l = l + 28;
		}
		// XY: recipients
		while (line[l] == '<' && line[l + 1] == 'm' && line[l + 6] == '>') {
			m++;//number of compatibilities
			l = l + 18;
			arc_out.push_back(node_labels[n - 1]);
			if (!is_next_char_digit(line, l)) {
				arc_in.push_back(char2uint(line[l]));
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					arc_in.push_back(10 * char2uint(line[l]) + char2uint(line[l + 1]));
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						arc_in.push_back(100 * char2uint(line[l]) + 10 * char2uint(line[l + 1]) + char2uint(line[l + 2]));
						l = l + 2;
					}
					else {
						//if (!is_next_char_digit(line, l + 3)){
						arc_in.push_back(1000 * char2uint(line[l]) + 100 * char2uint(line[l + 1]) + 10 * char2uint(line[l + 2]) + char2uint(line[l + 3]));
						l = l + 3;
						//}
					}
				}
			}
			l = l + 20;
			if (!is_next_char_digit(line, l)) {
				//arc_weight.push_back(char2uint(line[l]));
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					//arc_weight.push_back(10*char2uint(line[l])+char2uint(line[l+1]));
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						//arc_weight.push_back(100*char2uint(line[l])+10*char2uint(line[l+1])+char2uint(line[l+2]));
						l = l + 2;
					}
					else {
						////if (!is_next_char_digit(line, l + 3)){
						//arc_weight.push_back(1000*char2uint(line[l])+100*char2uint(line[l+1])+10*char2uint(line[l+2])+char2uint(line[l+3]));
						l = l + 3;
						////}
					}
				}
			}
			l = l + 17;
		}
		if (!(line[l] == '<' && line[l + 1] == 'e')) {
			l = l + 18;
		}
		if (n == no_of_nodes)
			break;
	}
	std::cout << "the number of nodes" << n;
	std::cout << "m: " << m << "\n";
	std::cout << "arc_in.size(): " << arc_in.size() << "\n" << " arc_out.size(): " << arc_out.size() << endl;
	return;
}

bool is_next_char_digit(string& line, unsigned int l) {
	if (line[l + 1] == '0' || line[l + 1] == '1' || line[l + 1] == '2' || line[l + 1] == '3' || line[l + 1] == '4' || line[l + 1] == '5' || line[l + 1] == '6' || line[l + 1] == '7' || line[l + 1] == '8' || line[l + 1] == '9')
		return true;
	return false;
}

unsigned int char2uint(char& p) {
	if (p == '1')
		return 1;
	else
		if (p == '2')
			return 2;
		else
			if (p == '3')
				return 3;
			else
				if (p == '4')
					return 4;
				else
					if (p == '5')
						return 5;
					else
						if (p == '6')
							return 6;
						else
							if (p == '7')
								return 7;
							else
								if (p == '8')
									return 8;
								else
									if (p == '9')
										return 9;
									else
										return 0;
}

double cpuTime() {
	return (double)clock() / CLOCKS_PER_SEC;
}



void ILP_d1_gurobi(unsigned short int& Q, unsigned short int& N, ListDigraph& g_original, unsigned short int& Vp, vector<unsigned short int>& node_arrives, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<pair<int, int>>& arc_pair, vector<int>& nodeset, vector<vector<unsigned short int>>& actual_alloc, double& M, double& M_total, vector<unsigned short int>& s, vector<pair<int, int>>& cycle_distri, vector<double>& target, vector<bool>& leaving, vector<double>& d, double& d_total, bool& c_involved, vector<double>& credit, map<int, int>& cycle_dis, bool lex_min, unsigned short int inst, std::map<int, std::map<int, int>>& cycle_dis_period) {
	pair_arcs(Q, g_original, node_arrives, active_nodes_original, arc_pair, nodeset);
	FilterNodes<ListDigraph> sg(g_original, active_nodes_original);
	long col_num = countArcs(sg) + 1;
	std::cout << "column: " << col_num - 1 << endl;
	std::cout << "countArcs(g_original): " << countArcs(g_original) << endl;
	M_total += M;
	//int node_number = countNodes(g_original);
	const unsigned short int row_num = N + 2 * nodeset.size() + 1;
	// Create an environment
	GRBEnv env = GRBEnv(true);
	env.set("LogFile", "mip_lexmin.log");

	env.start();

	// Create an empty model
	GRBModel model = GRBModel(env);
	//model.set("Threads", "16");
	vector<GRBModel> vector_model(2 * N - 1, GRBEnv(env));

	// Create variables
	vector<GRBVar> var_bi(col_num + 2 * N);
	for (unsigned short int i = 0; i < col_num - 1; ++i) {
		var_bi[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" + to_string(i));
	}
	var_bi[col_num - 1] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x" + to_string(col_num - 1));
	// N difference variables
	for (unsigned short int i = col_num; i < col_num + N; ++i) {
		var_bi[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x" + to_string(i));
	}
	// N abs variables
	for (unsigned short int i = col_num + N; i < col_num + 2 * N; ++i) {
		var_bi[i] = model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x" + to_string(i));
	}
	model.setObjective(1 * var_bi[col_num - 1], GRB_MINIMIZE);
	vector<double> bound(row_num, 0);
	for (unsigned short int i = 1; i < N + 1; ++i) {
		bound[i - 1] = target[i - 1] + credit[i - 1];
	}
	
	for (unsigned short int i = N + 1; i < N + 1 + nodeset.size(); ++i) {
		bound[i - 1] = 0;
	}
	for (unsigned short int i = N + nodeset.size() + 1; i < N + 2 * (nodeset.size()) + 1; ++i) {
		bound[i - 1] = 1;
	}
	bound[row_num - 1] = M;
	std::cout << "size of max-weighted perfect matching: " << M << endl;
	long matrix_num = row_num * col_num;
	vector<int> ia(matrix_num + 1, 0);
	//cout << "ia.size(): " << ia.size();
	vector<int> ja(matrix_num + 1, 0);
	vector<double> ar(matrix_num + 1, 0);

	//cout << "Period" << to_string(Q) << endl;
	//cout << "N: " << N << "col_num: " << col_num << "N*col_num: " << N * col_num;
	long cnt_2 = 0;

	//cout << ":arc_pair.size(): " << arc_pair.size() << endl;
	/*for (unsigned short int k = 0; k < arc_pair.size(); ++k) {
		cout << "arc_pair[k].first: " << arc_pair[k].first << " " << "arc_pair[k].second: " << arc_pair[k].second << endl;
	}*/
	// sum_y[i]+d1>=target[i]
	for (int i = 1; i < N + 1; ++i) {
		for (int k = 1; k < arc_pair.size() + 1; ++k) {
			if ((i - 1) * Vp <= arc_pair[k - 1].second && arc_pair[k - 1].second < i * Vp) {
				++cnt_2;
				ia[cnt_2] = i;
				ja[cnt_2] = k;
				ar[cnt_2] = 1.0;
			}

		}
	}

	int cnt_2_row = N;
	//sum_yij-sum_yji=0;
	for (int i = 0; i < nodeset.size(); ++i) {
		for (int k = 1; k < arc_pair.size() + 2; ++k) {
			if (k < arc_pair.size() + 1) {
				if (arc_pair[k - 1].first == nodeset[i]) {
					++cnt_2;
					ia[cnt_2] = cnt_2_row + i + 1;
					ja[cnt_2] = k;
					ar[cnt_2] = 1.0;
				}
				if (arc_pair[k - 1].second == nodeset[i]) {
					++cnt_2;
					ia[cnt_2] = cnt_2_row + i + 1;
					ja[cnt_2] = k;
					ar[cnt_2] = -1.0;
				}
			}

		}
	}



	int cnt_3_row = N + nodeset.size();
	//sum_yij<=1;
	for (int i = 0; i < nodeset.size(); ++i) {
		for (int k = 1; k < arc_pair.size() + 2; ++k) {
			if (k < arc_pair.size() + 1) {
				if (arc_pair[k - 1].second == nodeset[i]) {
					++cnt_2;
					ia[cnt_2] = cnt_3_row + i + 1;
					ja[cnt_2] = k;
					ar[cnt_2] = 1.0;
				}
			}


		}
	}


	int cnt_4_row = N + 2 * nodeset.size();
	//sum_yij=M
	for (int k = 1; k < arc_pair.size() + 2; ++k) {
		if (k < arc_pair.size() + 1) {
			++cnt_2;
			ia[cnt_2] = row_num;
			ja[cnt_2] = k;
			ar[cnt_2] = 1.0;
		}

	}

	vector<GRBLinExpr> sum_row(row_num, 0);

	for (long j = 1; j < cnt_2 + 1; ++j) {
		sum_row[ia[j] - 1] += ar[j] * var_bi[ja[j] - 1];
	}
	std::cout << "finish loading efficiencies" << endl;
	for (unsigned short int i = 1; i < N + 1; ++i) {
		model.addConstr(var_bi[col_num + i - 1] == sum_row[i - 1] - bound[i - 1]);
		model.addGenConstrAbs(var_bi[col_num + N + i - 1], var_bi[col_num + i - 1]);
		model.addConstr(var_bi[col_num + N + i - 1] <= var_bi[col_num - 1]);
	}

	for (unsigned short int i = N + 1; i < N + 1 + nodeset.size(); ++i) {
		model.addConstr(sum_row[i - 1] == bound[i - 1]);
	}
	for (unsigned short int i = N + nodeset.size() + 1; i < N + 2 * (nodeset.size()) + 1; ++i) {
		model.addConstr(sum_row[i - 1] <= bound[i - 1]);
	}
	model.addConstr(sum_row[row_num - 1] == bound[row_num - 1]);


	model.optimize();

	//lexicographically minimal maximum cycle packing


	vector<double> d_t(N, 0);
	d_t[0] = var_bi[col_num - 1].get(GRB_DoubleAttr_X);
	std::cout << "d_t[0]" << d_t[0] << endl;
	vector<unsigned short int> N_star(N, 0);
	unsigned short int n_star = 0;
	unsigned short int t_star = 0;
	double epsilon = 0;
	unsigned short int track = 0;
	vector<GRBVar> var_lexmin(arc_pair.size());
	if (d_t[0] > 0.5 && lex_min) {
		epsilon_func(target, credit, epsilon, N);
		sort_d_t(d_t, var_bi, col_num, N, Vp, arc_pair, target, t_star, credit, epsilon, var_lexmin, N_star);
	}
	std::cout << "finish sorting" << "epsilon:" << epsilon << endl;
	std::cout << "start n_star_1" << endl;
	if (lex_min && d_t[0] > 0.5 && abs(epsilon) > pow(10, -4)) {
		lex_min_n_star(d_t, lex_min, t_star, N, col_num, epsilon, n_star, model, ia, ja, ar, row_num, cnt_2, bound, nodeset, N_star, var_lexmin, vector_model, track);
	}
	std::cout << "finish n_star_1" << endl;
	lexmin_searching(d_t, lex_min, t_star, N, col_num, epsilon, n_star, model, ia, ja, ar, row_num, cnt_2, bound, nodeset, N_star, Vp, arc_pair, target, credit, var_lexmin, inst, vector_model, track);

	cout << "finish lexmin_searching" << endl;

	unsigned short int t = 0;
	for (unsigned short int i = 1; i < arc_pair.size() + 1; ++i) {
		//cout << "{" << arc_pair[i].first << "," << arc_pair[i].second << "}" << endl;
		//cout << var_bi[i - 1].get(GRB_DoubleAttr_X);
		if (lex_min && d_t[0] > 0.5) {
			if (var_lexmin[i - 1].get(GRB_DoubleAttr_X) > pow(10, -4)) {
				++t;
				leaving[arc_pair[i - 1].first] = true;
				leaving[arc_pair[i - 1].second] = true;
				cycle_distri.push_back(arc_pair[i - 1]);
				for (unsigned short int j = 0; j < N; ++j) {
					if (j * Vp <= arc_pair[i - 1].second && arc_pair[i - 1].second < (j + 1) * Vp) {
						++s[j];
					}
				}
			}

		}
		else {
			if (var_bi[i - 1].get(GRB_DoubleAttr_X) > pow(10, -4)) {
				++t;
				leaving[arc_pair[i - 1].first] = true;
				leaving[arc_pair[i - 1].second] = true;
				cycle_distri.push_back(arc_pair[i - 1]);
				for (unsigned short int j = 0; j < N; ++j) {
					if (j * Vp <= arc_pair[i - 1].second && arc_pair[i - 1].second < (j + 1) * Vp) {
						++s[j];
					}
				}
			}
		}

	}
	model.reset();
	model.update();
	for (unsigned short int i = 0; i < track + 1; ++i) {
		vector_model[i].reset();
		vector_model[i].update();
	}

	/*for (unsigned short int i = 0; i < cycle_distri.size(); ++i) {
		cout << "{" << cycle_distri[i].first << "," << cycle_distri[i].second << "}" << endl;
	}*/

	cycle_distribution(cycle_dis_period, cycle_dis, cycle_distri, N, Q);
	std::cout << "maximum size: " << t;
	for (unsigned short int i = 0; i < N; ++i) {
		d[i] += target[i] - s[i];
		if (c_involved) {
			credit[i] += target[i] - s[i];
		}
		else {
			credit[i] = 0;
		}
		std::cout << "country " << to_string(i)<<" " << "initial allocation " <<i<<": " << target[i] << '/n' << "s[i] : " << s[i] << "d[i] : " << d[i] << "credit[i] : " << credit[i] << endl;
		//actual_alloc[Q].push_back(s[i]);
	}
	return;
}


void arbitraryMaximum(vector<unsigned short int>& node_arrives, ListGraph& g, ListDigraph& g_original, vector<pair<int, int>>& arc_pair, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<double>& v, bool& target_omega, vector<double>& target, vector<double>& credit, ListGraph::EdgeMap<double>& edge_card_weight, double& t0, vector<vector<unsigned short int>>& actual_alloc, vector<double>& v_impu, vector<int>& nodeset, vector<pair<int, int>>& cycle_distri, vector<double>& d, double& M_total, double& d_total, bool& c_involved, map<int, int>& cycle_dis, map<int, int>& numofMaxSolution, bool& arbitray_maximum, unsigned short int& initialSize, vector<double>& v_S, unsigned int& S, double& core_100, long& negative_core, unsigned short int inst, double& max_d, double& game_generation, double& solution_concept_time, std::map<int, std::map<int, int>>& cycle_dis_arbitrary_period) {
	Q = 0;
	d_total = 0;
	M_total = 0;
	core_100 = 0;
	negative_core = 0;
	max_d = 0;
	if (dispy)
		std::cout << " --== Without lex min matching == -- " << endl;
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		credit[i] = 0;
		no_of_active_nodes[i] = initialSize;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (node_arrives[i * Vp + j] == 0) {
				active_nodes[c[i * Vp + j]] = true;
				active_nodes[c_b[i * Vp + j]] = true;
				active_nodes_original[c_original[i * Vp + j]] = true;
			}
			else {
				active_nodes[c[i * Vp + j]] = false;
				active_nodes[c_b[i * Vp + j]] = false;
				active_nodes_original[c_original[i * Vp + j]] = false;
			}
		}
	}
	while (Q < periods) {
		if (dispy) {
			std::cout << "--== PERIOD " << Q + 1 << " ==--" << endl;
		}
		if (dispy) {
			std::cout << "Number of active nodes: ";
			for (unsigned short int i = 0; i < N; i++)
				std::cout << no_of_active_nodes[i] << " ";
			std::cout << endl;
		}
		// cooperative game and target
		std::cout << "start generating values" << endl;
		coop_game(g, v, v_impu, v_S, S, s, c, c_b, edge_card_weight, dispy, Vp, N, active_nodes, leaving, numofMaxSolution, Q, arbitray_maximum, cycle_distri, cycle_dis, game_generation, cycle_dis_arbitrary_period);
		double t0 = cpuTime();
		double suma = 0;
		double sumimpu = 0;
		if (target_omega) {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i] - v_impu[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i] - v_impu[i]) / suma);
		}
		else {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i]) / suma);
		}
		double t1 = cpuTime();
		solution_concept_time += t1 - t0;
		M_total += v[N];

		//compute deviations
		if (dispy) {
			if (target_omega) {
				std::cout << "Benefit: ";
			}
			else {
				std::cout << "Contribution: ";
			}
			for (unsigned short int i = 0; i < N; i++) {
				std::cout << target[i] << " ";
			}
			std::cout << endl;
		}


		for (unsigned short int i = 0; i < N; ++i) {
			d[i] += target[i] - s[i];
			if (c_involved) {
				credit[i] += target[i] - s[i];
			}
			else {
				credit[i] = 0;
			}
			std::cout << "country" << to_string(i) << "target[i]: " << target[i] << '/n' << "s[i]: " << s[i] << "d[i]: " << d[i] << "credit[i]: " << credit[i] << endl;
			//actual_alloc[Q].push_back(s[i]);
		}


		Q++;
		changing_nodes(active_nodes, active_nodes_original, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, c_b, c_original, s, d, target);
		if (dispy)
			cin.get();
	}
	for (unsigned short int i = 0; i < N; ++i) {
		d_total += abs(d[i]);
	}
	vector<double> max_deviation(N, 0);
	for (unsigned short int i = 0; i < N; ++i) {
		max_deviation[i] = abs(d[i]) / M_total;
	}
	std::sort(max_deviation.begin(), max_deviation.end());
	max_d = max_deviation[N - 1];
	return;
}


void pair_arcs(unsigned short int& Q, ListDigraph& g_original, vector<unsigned short int>& node_arrives, ListDigraph::NodeMap<bool>& active_nodes_original, vector<pair<int, int>>& arc_pair, vector<int>& nodeset) {
	FilterNodes<ListDigraph> sg(g_original, active_nodes_original);
	arc_pair.clear();
	nodeset.clear();
	for (FilterNodes<ListDigraph>::NodeIt n(sg); n != INVALID; ++n) {
		nodeset.push_back(sg.id(n));
	}
	sort(nodeset.begin(), nodeset.end());
	for (FilterNodes<ListDigraph>::ArcIt a(sg); a != INVALID; ++a) {
		arc_pair.push_back({ sg.id(sg.source(a)),sg.id(sg.target(a)) });
	}
	sort(arc_pair.begin(), arc_pair.end());
	return;

}

void sort_d_t(vector<double>& d_t, vector<GRBVar>& var_bi, long& col_num, unsigned short int& N, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, unsigned short int& t, vector<double>& credit, double& epsilon, vector<GRBVar>& var_lexmin, vector<unsigned short int>& N_star) {
	unsigned short int n_star = 0;
	vector<unsigned short int> s_copy(N, 0);
	vector<double> d_copy(N, 0);
	for (unsigned short int i = 1; i < arc_pair.size() + 1; ++i) {
		var_lexmin[i - 1] = var_bi[i - 1];
		//cout << "{" << arc_pair[i].first << "," << arc_pair[i].second << "}" << endl;
		//cout << var_bi[i - 1].get(GRB_DoubleAttr_X);
		if (var_bi[i - 1].get(GRB_DoubleAttr_X) > pow(10, -4)) {
			for (unsigned short int j = 0; j < N; ++j) {
				if (j * Vp <= arc_pair[i - 1].second && arc_pair[i - 1].second < (j + 1) * Vp) {
					++s_copy[j];
				}
			}
		}
	}
	for (unsigned short int i = 0; i < N; ++i) {
		d_copy[i] = abs(target[i] + credit[i] - s_copy[i]);
	}
	//sort
	std::sort(d_copy.begin(), d_copy.end());
	/*for (unsigned short int i = 0; i < N; ++i) {
		std::cout << "d_copy[i]" << d_copy[i] << endl;
		std::cout << "s_copy" << s_copy[i] << endl;
	}*/
	d_t[t] = d_copy[N - 1 - t];

	if (t < N - 1) {
		d_t[t + 1] = d_copy[N - 2 - t];
	}
	std::cout << "d" << to_string(t - 1) << d_t[t - 1] << endl;
	std::cout << "epsilon: " << epsilon << endl;
}

void lex_min_n_star(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, vector<GRBVar>& var_lexmin, vector<GRBModel>& vector_model, unsigned short int& track) {
	++track;
	std::cout << "t_star: " << t_star << endl;
	vector<GRBVar> var_bi(col_num + N * (t_star + 3));
	std::cout << "var_bi.size(): " << var_bi.size() << endl;
	//arc variables
	for (unsigned short int i = 0; i < col_num - 1; ++i) {
		var_bi[i] = vector_model[track].addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" + to_string(i));
	}
	// N difference variables
	for (unsigned short int i = col_num; i < col_num + N; ++i) {
		var_bi[i] = vector_model[track].addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x" + to_string(i));
	}
	// N abs variables
	for (unsigned short int i = col_num + N; i < col_num + N * 2; ++i) {
		var_bi[i] = vector_model[track].addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x" + to_string(i));
	}
	// z^t_p variables
	for (unsigned short i = col_num + N * 2; i < col_num + N * (t_star + 3); ++i) {
		var_bi[i] = vector_model[track].addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" + to_string(i));
	}

	//sum zp
	GRBLinExpr sum_zp = 0;
	for (unsigned short int i = col_num + N * (t_star + 2); i < col_num + N * (t_star + 3); ++i) {
		sum_zp += var_bi[i];
	}
	vector_model[track].setObjective(sum_zp, GRB_MINIMIZE);
	vector<GRBLinExpr> sum_row(row_num, 0);
	vector<GRBLinExpr> sum_zp_dp(N, 0);
	vector<GRBLinExpr> sum_zp_0(N, 0);
	vector<GRBLinExpr> sum_t_star(t_star + 1, 0);
	vector<GRBLinExpr> sum_p(N, 0);
	for (unsigned short int i = 0; i < t_star; ++i) {
		for (unsigned short int j = 0; j < N; ++j) {
			sum_zp_dp[j] += var_bi[i * N + col_num + N * 2 + j] * d_t[i];
			sum_zp_0[j] += var_bi[i * N + col_num + N * 2 + j];
			sum_t_star[i] += var_bi[i * N + col_num + N * 2 + j];
			sum_p[j] += var_bi[i * N + col_num + N * 2 + j];
		}
	}
	for (unsigned short int j = 0; j < N; ++j) {
		sum_p[j] += var_bi[t_star * N + col_num + N * 2 + j];
	}
	std::cout << "finish loading z^t_p" << endl;
	for (long j = 1; j < cnt_2 + 1; ++j) {
		sum_row[ia[j] - 1] += ar[j] * var_bi[ja[j] - 1];
	}
	std::cout << "finish loading coefficiencies" << endl;
	for (unsigned short int i = 1; i < N + 1; ++i) {
		vector_model[track].addConstr(var_bi[col_num + i - 1] == sum_row[i - 1] - bound[i - 1]);
		vector_model[track].addGenConstrAbs(var_bi[col_num + N + i - 1], var_bi[col_num + i - 1]);
		if (t_star == 0) {
			vector_model[track].addConstr(var_bi[col_num + N + i - 1] <= d_t[t_star] - epsilon * (1 - var_bi[col_num + 2 * N + t_star * N + i - 1]));
		}
		if (t_star > 0) {
			vector_model[track].addQConstr(var_bi[col_num + N + i - 1] <= (1 - sum_zp_0[i - 1]) * (d_t[t_star] - epsilon * (1 - var_bi[col_num + 2 * N + t_star * N + i - 1])) + sum_zp_dp[i - 1]);
		}
	}
	//std::cout << "finish loading constraint 1" << endl;
	for (unsigned short int i = N + 1; i < N + 1 + nodeset.size(); ++i) {
		vector_model[track].addConstr(sum_row[i - 1] == bound[i - 1]);
	}
	//std::cout << "finish loading constraint 2" << endl;
	for (unsigned short int i = N + nodeset.size() + 1; i < N + 2 * (nodeset.size()) + 1; ++i) {
		vector_model[track].addConstr(sum_row[i - 1] <= bound[i - 1]);
	}
	//std::cout << "finish loading constraint 3" << endl;
	vector_model[track].addConstr(sum_row[row_num - 1] == bound[row_num - 1]);
	//std::cout << "finish loading constraint 4" << endl;
	for (unsigned short int i = 0; i < N; ++i) {
		vector_model[track].addConstr(sum_p[i] <= 1);
	}
	//std::cout << "finish loading constraint 5" << endl;
	for (unsigned short int i = 0; i < t_star; ++i) {
		vector_model[track].addConstr(sum_t_star[i] == N_star[i]);
	}
	//std::cout << "finish loading constraint 6" << endl;
	std::cout << "finish loading constraints" << endl;
	vector_model[track].optimize();
	for (unsigned short int i = 0; i < t_star + 1; ++i) {
		for (unsigned short int j = 0; j < N; ++j) {
			std::cout << "t: " << i << " " << "country:" << j << ": " << var_bi[col_num + N * (2 + i) + j].get(GRB_DoubleAttr_X) << endl;
		}
	}
	for (unsigned short int i = 0; i < N; ++i) {
		if (var_bi[col_num + N * (2 + t_star) + i].get(GRB_DoubleAttr_X) > pow(10, -4)) {
			N_star[t_star] += 1;
		}
		//std::cout << "var_bi[col_num + N * (2 + t_star) + i].get(GRB_DoubleAttr_X): " << var_bi[col_num + N * (2 + t_star) + i].get(GRB_DoubleAttr_X) << endl;

	}
	std::cout << "epsilon: " << epsilon << "d_t[t_star: " << d_t[t_star] << endl;
	std::cout << "N_star[t_star]: " << N_star[t_star] << endl;
	std::cout << "finish n_" << to_string(t_star) << endl;
	n_star = 0;
	for (unsigned short int i = 0; i < N; ++i) {
		n_star += N_star[i];
	}
	t_star += N_star[t_star];
	for (unsigned short int i = 0; i < col_num - 1; ++i) {
		var_lexmin[i] = var_bi[i];
	}

}

void lex_min_d_star(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, vector<double>& credit, vector<GRBVar>& var_lexmin, vector<GRBModel>& vector_model, unsigned short int& track) {
	++track;
	vector<GRBVar> var_bi(col_num + N * (t_star + 3));
	//arc variables
	for (unsigned short int i = 0; i < col_num - 1; ++i) {
		var_bi[i] = vector_model[track].addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" + to_string(i));
	}
	// d_t variable
	var_bi[col_num - 1] = vector_model[track].addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "d" + to_string(t_star));
	// N difference variables
	for (unsigned short int i = col_num; i < col_num + N; ++i) {
		var_bi[i] = vector_model[track].addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x" + to_string(i));
	}
	// N abs variables
	for (unsigned short int i = col_num + N; i < col_num + N * 2; ++i) {
		var_bi[i] = vector_model[track].addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x" + to_string(i));
	}
	//zp variables
	for (unsigned short i = col_num + N * 2; i < col_num + N * (t_star + 2); ++i) {
		var_bi[i] = vector_model[track].addVar(0.0, 1.0, 0.0, GRB_BINARY, "x" + to_string(i));
	}
	vector_model[track].setObjective(1 * var_bi[col_num - 1], GRB_MINIMIZE);
	vector<GRBLinExpr> sum_row(row_num, 0);
	vector<GRBLinExpr> sum_zp_dp(N, 0);
	vector<GRBLinExpr> sum_zp_0(N, 0);
	vector<GRBLinExpr> sum_t_star(t_star + 1, 0);
	vector<GRBLinExpr> sum_p(N, 0);
	for (unsigned short int i = 0; i < t_star; ++i) {
		for (unsigned short int j = 0; j < N; ++j) {
			sum_zp_dp[j] += var_bi[i * N + col_num + N * 2 + j] * d_t[i];
			sum_zp_0[j] += var_bi[i * N + col_num + N * 2 + j];
			sum_t_star[i] += var_bi[i * N + col_num + N * 2 + j];
			sum_p[j] += var_bi[i * N + col_num + N * 2 + j];
		}
	}


	for (long j = 1; j < cnt_2 + 1; ++j) {
		sum_row[ia[j] - 1] += ar[j] * var_bi[ja[j] - 1];
	}

	for (unsigned short int i = 1; i < N + 1; ++i) {
		vector_model[track].addConstr(var_bi[col_num + i - 1] == sum_row[i - 1] - bound[i - 1]);
		vector_model[track].addGenConstrAbs(var_bi[col_num + N + i - 1], var_bi[col_num + i - 1]);
		if (t_star == 0) {
			vector_model[track].addConstr(var_bi[col_num + N + i - 1] <= var_bi[col_num - 1] + sum_zp_dp[i - 1]);
		}
		if (t_star > 0) {
			vector_model[track].addQConstr(var_bi[col_num + N + i - 1] <= (1 - sum_zp_0[i - 1]) * var_bi[col_num - 1] + sum_zp_dp[i - 1]);
		}
	}

	for (unsigned short int i = N + 1; i < N + 1 + nodeset.size(); ++i) {
		vector_model[track].addConstr(sum_row[i - 1] == bound[i - 1]);
	}
	for (unsigned short int i = N + nodeset.size() + 1; i < N + 2 * (nodeset.size()) + 1; ++i) {
		vector_model[track].addConstr(sum_row[i - 1] <= bound[i - 1]);
	}
	vector_model[track].addConstr(sum_row[row_num - 1] == bound[row_num - 1]);

	for (unsigned short int i = 0; i < N; ++i) {
		vector_model[track].addConstr(sum_p[i] <= 1);
	}
	for (unsigned short int i = 0; i < t_star; ++i) {
		vector_model[track].addConstr(sum_t_star[i] == N_star[i]);
	}

	vector_model[track].optimize();
	sort_d_t(d_t, var_bi, col_num, N, Vp, arc_pair, target, t_star, credit, epsilon, var_lexmin, N_star);
}

void lexmin_searching(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, vector<double>& credit, vector<GRBVar>& var_lexmin, unsigned short int inst, vector<GRBModel>& vector_model, unsigned short int& track) {
	while (lex_min && abs(d_t[t_star - 1]) > 0.5 && n_star < N) {
		std::cout << "begin search d_t" << endl;
		lex_min_d_star(d_t, lex_min, t_star, N, col_num, epsilon, n_star, model, ia, ja, ar, row_num, cnt_2, bound, nodeset, N_star, Vp, arc_pair, target, credit, var_lexmin, vector_model, track);
		std::cout << "inst: " << inst << endl;
		std::cout << "abs(d_t[t_star])" << abs(d_t[t_star]) << endl;
		std::cout << "t_star: " << t_star << endl;
		std::cout << "N-1: " << N - 1 << endl;
		if (t_star == N - 1) {
			std::cout << "N-1==t_star" << endl;
		}
		else {
			std::cout << "N-1!=t_star" << endl;
		}
		std::cout << "epsilon in the loop: " << epsilon << endl;
		std::cout << "absolute epsilon in the loop: " << abs(epsilon) << endl;
		if (abs(d_t[t_star]) > 0.5) {
			if (abs(epsilon) > pow(10, -4)) {
				if (t_star == N - 1) {
					std::cout << "congratulations t_star == N - 1" << endl;
					//break
					n_star = N;
					std::cout << "congratulations after break" << endl;
				}
				else
				{
					if (d_t[t_star + 1] < 0.5) {
						std::cout << "congratulations d_t[t_star + 1] < 0.5" << endl;
						//break;
						n_star = N;
					}
					else {
						std::cout << "congratulations d_t[t_star + 1] > 0.5" << endl;
						lex_min_n_star(d_t, lex_min, t_star, N, col_num, epsilon, n_star, model, ia, ja, ar, row_num, cnt_2, bound, nodeset, N_star, var_lexmin, vector_model, track);
					}
				}
			}
			if (abs(epsilon) < pow(10, -4)) {
				std::cout << "congratulations abs(epsilon) < pow(10, -7))" << endl;
				lexmin_searching(d_t, lex_min, t_star, N, col_num, epsilon, n_star, model, ia, ja, ar, row_num, cnt_2, bound, nodeset, N_star, Vp, arc_pair, target, credit, var_lexmin, inst, vector_model, track);
			}
		}
		else {
			std::cout << "congratulations break" << endl;
			//break
			n_star = N;
		}
	}
}


void cycle_distribution(std::map<int, std::map<int, int>>& cycle_dis_period, map<int, int>& cycle_dis, vector<pair<int, int>>& cycle_distri, unsigned short int& N, unsigned short int& Q) {
	unsigned short int i = 0;
	while (i < cycle_distri.size()) {
		int first, last;
		unsigned short int cycle_length = 1;
		first = cycle_distri[i].first;
		last = cycle_distri[i].second;
		cycle_distri.erase(cycle_distri.begin());
		//cout << cycle_distri.size() << '\n' << i << endl;
		for (unsigned short int j = i; j < cycle_distri.size(); ++j) {
			if (last == cycle_distri[j].first) {
				++cycle_length;
				if (first == cycle_distri[j].second) {
					//cout << "period " << j << " " << cycle_length << endl;
					++cycle_dis_period[N * (Q + 1)][cycle_length];
					//cout << "cycle_dis_period[N*(Q+1)][cycle_length]: " << cycle_dis_period[N * (Q + 1)][cycle_length] << endl;
					++cycle_dis[cycle_length];
					cycle_length = 1;
					cycle_distri.erase(cycle_distri.begin() + j);
					i = 0;
					//cout << "i:" << i << '\n' << "j: " << j << endl;
					break;
				}
				else {
					last = cycle_distri[j].second;
					//cout << last << endl;
					cycle_distri.erase(cycle_distri.begin() + j);
					j = -1;
				}

			}
		}
	}
}

void epsilon_func(vector<double> target, vector<double> credit, double& epsilon, unsigned short int N) {
	vector<double> target_credit(N, 0);
	vector<double> epsilon_sort(N * (N - 1), 0);
	for (unsigned short int i = 0; i < N; ++i) {
		target_credit[i] = target[i] + credit[i];
	}

	unsigned short int t = -1;
	for (unsigned short int i = 0; i < N - 1; ++i) {
		for (unsigned short j = i + 1; j < N; ++j) {
			++t;
			//cout << "target_credit[i]: " << target_credit[i] << "target_credit[j]: " << target_credit[j] << endl;
			epsilon_sort[t] = abs(frac(target_credit[i]) - frac(target_credit[j]));
			//cout << "a-b: " << epsilon_sort[t] << endl;
			++t;
			epsilon_sort[t] = abs(frac(target_credit[i]) - (1 - frac(target_credit[j])));
			//cout << "a-(1-b): " << epsilon_sort[t] << endl;
		}
	}
	cout << "t" << t << endl;


	auto newEnd = std::remove_if(epsilon_sort.begin(), epsilon_sort.end(), [](double num) {
		return num < 2 * pow(10, -4);
		});
	epsilon_sort.erase(newEnd, epsilon_sort.end());

	std::sort(epsilon_sort.begin(), epsilon_sort.end());
	epsilon = epsilon_sort[0];
	cout << "epsilon_sort[0]" << epsilon_sort[0] << endl;
}

double frac(double ori) {
	double abs_frac;
	abs_frac = abs(ori) - abs(int(ori));
	return abs_frac;
}

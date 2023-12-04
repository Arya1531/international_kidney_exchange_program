/*
*    int_kidney_exchange
*    banzhaf.cpp
*    Purpose: computational study for Computing Balanced Solutions for Large International Kidney
*			  Exchange Schemes When Cycle Length Is Unbounded
*             using the banzhaf value as initial allocations with equal country sizes
*
*
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
#include <lemon/list_graph.h> // needed for ListDigraph
#include <lemon/matching.h>
#include <lemon/adaptors.h>
#include <lemon/core.h>
#include <lemon/base.cc>
#include <time.h>
#include <glpk.h>
#include <iomanip>

#include<cfloat>
#include <math.h>
#include <stdio.h>
#include <iterator>
#include "gurobi_c++.h"

using namespace lemon;
using namespace std;

void coop_game(ListGraph& g, vector<double>& v, unsigned int& S, vector<unsigned short int>& s, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, ListGraph::EdgeMap<double>& edge_card_weight, bool& dispy, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, map<int, int>& numofMaxSolution, unsigned short int& Q, bool& arbitray_maximum, vector<pair<int, int>>& cycle_distri, map<int, int>& cycle_dis, bool& unique_impu, double& prec, bool& unique_imputation, bool& credit_adjusted, vector<double>& credit, double& M, double& game_generation, std::map<int, std::map<int, int >>& cycle_dis_arbitrary_period);
void norm_banzhaf(vector<double>& banz, vector<double>& v, unsigned short int& n, unsigned int& s);
void de2bi_card(unsigned int& k, vector<bool>& a, unsigned short int& n, unsigned short int& card);
void de2bi(unsigned int& k, vector<bool>& a, unsigned short int& n);
double core_distance(vector<double>& x, vector<double>& v, unsigned short int& n, unsigned int& S);
void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, vector<unsigned short int>& s, vector<double>& d);
void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, unsigned short int& initialSize);
void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods, bool& disp);
void undi_lemon(unsigned int& m, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, vector<unsigned short int>& label_positions, ListGraph& g, ListDigraph& g_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, ListGraph::EdgeMap<double>& edge_card_weight, ListDigraph::ArcMap<unsigned short int>& arc_card_weight, unsigned short int& no_of_nodes);
void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, unsigned short int& k, ListGraph& g, ListDigraph& g_original, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, unsigned int& m, unsigned short int& no_of_nodes);
bool is_next_char_digit(string& line, unsigned int l);
unsigned int char2uint(char& p);
double cpuTime();

void arbitraryMaximum(vector<unsigned short int>& node_arrives, ListGraph& g, ListDigraph& g_original, vector<pair<int, int>>& arc_pair, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<double>& v, vector<double>& credit, ListGraph::EdgeMap<double>& edge_card_weight, double& t0, vector<int>& nodeset, vector<pair<int, int>>& cycle_distri, vector<double>& d, double& M_total, double& d_total, bool& c_involved, map<int, int>& cycle_dis, map<int, int>& numofMaxSolution, bool& arbitray_maximum, unsigned short int& initialSize, unsigned int& S, double& core_100, long& negative_core, double& prec, vector<double>& init_alloc, unsigned short int inst, bool& unique_impu, bool& unique_imputation, bool& core_dist, vector<double>& init_alloc_accum, vector<double>& s_accum, vector<double>& v_accum, double& y_core_dist, double& s_core_dist, bool& credit_adjusted, double& M, double& max_d, double& game_generation, double& solution_concept_time, std::map<int, std::map<int, int>>& cycle_dis_arbitrary_period);
void period_0(unsigned short int& Q, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, vector<unsigned short int>& s, unsigned short int& Vp, vector<unsigned short int>& node_arrives, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, vector<double>& credit);
void min_d_1(vector<unsigned short int>& node_arrives, ListGraph& g, ListDigraph& g_original, vector<pair<int, int>>& arc_pair, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<double>& v, vector<double>& credit, ListGraph::EdgeMap<double>& edge_card_weight, double& t0, vector<int>& nodeset, vector<pair<int, int>>& cycle_distri, vector<double>& d, double& M_total, double& d_total, bool& c_involved, map<int, int>& cycle_dis, map<int, int>& numofMaxSolution, bool& arbitray_maximum, unsigned short int& initialSize, unsigned int& S, double& core_100, long& negative_core, double& prec, vector<double>& init_alloc, double& d_c_total, unsigned short int inst, bool& unique_impu, bool& unique_imputation, bool& core_dist, vector<double>& init_alloc_accum, vector<double>& s_accum, vector<double>& v_accum, double& y_core_dist, double& s_core_dist, bool& credit_adjusted, double& M, bool lex_min, double& max_d, double& game_generation, double& solution_concept_time, double& scenario_time, std::map<int, std::map<int, int>>& cycle_dis_period, std::map<int, std::map<int, int>>& cycle_dis_arbitrary_period);
void ILP_d1_gurobi(unsigned short int& Q, unsigned short int& N, ListDigraph& g_original, unsigned short int& Vp, vector<unsigned short int>& node_arrives, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<pair<int, int>>& arc_pair, vector<int>& nodeset, double& M, double& M_total, vector<unsigned short int>& s, vector<pair<int, int>>& cycle_distri, vector<bool>& leaving, vector<double>& d, double& d_total, bool& c_involved, vector<double>& credit, map<int, int>& cycle_dis, vector<double>& init_alloc, bool& credit_adjusted, bool lex_min, unsigned short int inst, std::map<int, std::map<int, int>>& cycle_dis_period);
void pair_arcs(unsigned short int& Q, ListDigraph& g_original, vector<unsigned short int>& node_arrives, ListDigraph::NodeMap<bool>& active_nodes_original, vector<pair<int, int>>& arc_pair, vector<int>& nodeset);
void cycle_distribution(std::map<int, std::map<int, int>>& cycle_dis_period, map<int, int>& cycle_dis, vector<pair<int, int>>& cycle_distri, unsigned short int& N, unsigned short int& Q);
void lex_min_d_star(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, vector<double>& credit, vector<GRBVar>& var_lexmin, vector<GRBModel>& vector_model, unsigned short int& track, bool& credit_adjusted);
void lex_min_n_star(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, vector<GRBVar>& var_lexmin, vector<GRBModel>& vector_model, unsigned short int& track);
void sort_d_t(vector<double>& d_t, vector<GRBVar>& var_bi, long& col_num, unsigned short int& N, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, unsigned short int& t, vector<double>& credit, double& epsilon, vector<GRBVar>& var_lexmin, vector<unsigned short int>& N_star, bool& credit_adjusted);
void lexmin_searching(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, vector<double>& credit, vector<GRBVar>& var_lexmin, unsigned short int inst, vector<GRBModel>& vector_model, unsigned short int& track, bool& credit_adjusted);
void epsilon_func(vector<double>& target, vector<double>& credit, double& epsilon, unsigned short int N, bool& credit_adjusted);
double frac(double ori);

int main() {
	cout << "I solemnly swear that I am up to no good." << endl;
	bool c_involved = false;// true: credits considred; false:without credits 
	bool arbitray_maximum = false; //true: arbitray maximum cycple packing
	bool lex_min = false;
	string solution_concept = "banzhaf";
	string version;
	bool d1 = false;
	bool d_c = false;
	bool arbitrary = true;
	bool lexmin_call = false;
	bool lexmin_c_call = false;
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
			version = "d1_c";
		}
		if (arbitrary) {
			version = "arbitrary";
		}
	}
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
	unsigned short int N; // number of countries/players
	unsigned short int inst; // instance number, integer between 0 and 99
	bool dispy = false; // true: information in terminal while running
	bool disp = false; // true: extremely detailed information while running, avoid with large graphs
	bool unique_imputation = false;
	bool core_dist = true;
	bool credit_adjusted = false;
	double M = 0;
	for (N = 4; N < 11; ++N) {
		cycle_dis.clear();
		cycle_dis_d.clear();
		cycle_dis_t_c.clear();
		cycle_dis_arbitrary.clear();
		double relative_d1 = 0;
		double relative_d1_c = 0;
		double relative_d1_arbitrary = 0;
		double relative_d1_initial_allocation = 0;
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
			unsigned short int years = 6;
			unsigned short int periods_per_year = 4;
			// read the data
			string line;
			ifstream inp;
			unsigned short int graph_size = 2000;
			inp.open("/home/kllg37/data/genxml-" + to_string(inst) + ".xml"); // 1 out of the 100 instances generated by William Pettersson's web tool: https://wpettersson.github.io/kidney-webapp/#/
			getline(inp, line);
			inp.close();

			unsigned short int Vp = 4 * (unsigned short int)((graph_size / 4) / N);
			unsigned short int no_of_nodes = N * Vp;
			unsigned short int initialSize = Vp / 4;
			vector<unsigned int> arc_out(0, 0);
			vector<unsigned int> arc_in(0, 0);
			double M_total = 0;//changed by XY: maximum size of cycle packing
			unsigned int m = 0;
			unsigned short int k = 0;
			vector<unsigned short int> node_labels(no_of_nodes, 0);
			vector<unsigned short int> label_positions(graph_size, graph_size + 1);
			ListGraph g;
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
			vector<unsigned short int> no_of_active_nodes(N, Vp / 4);
			ListGraph::NodeMap<bool> active_nodes(g);
			ListDigraph::NodeMap<bool> active_nodes_original(g_original);
			for (unsigned short int i = 0; i < N; i++) {
				for (unsigned short int j = 0; j < Vp; j++) {
					active_nodes[c[i * Vp + j]] = false;
					active_nodes[c_b[i * Vp + j]] = false;
					active_nodes_original[c_original[i * Vp + j]] = false;
				}
			}

			// read the seed
			string line_seed;
			ifstream seed_doc;
			cout << "start reading the seed doc" << endl;
			seed_doc.open("/home/kllg37/seeds/n" + to_string(N) + "inst" + to_string(inst) + ".txt");
			getline(seed_doc, line_seed);
			seed_doc.close();
			unsigned int seed = 0;
			seed = stoi(line_seed); 
			srand(seed);

			// determining starting pairs and arrival times of others
			ofstream res;
			initial_pairs(Vp, N, active_nodes, c, disp, active_nodes_original, c_b, c_original, initialSize);
			vector<unsigned short int> node_arrives(no_of_nodes, 0);
			arrival_times(node_arrives, Vp, N, active_nodes, c, periods, disp);
			if (disp) {
				cout << endl;
				vector<unsigned short int> bla(periods, 0);
				for (unsigned short int i = 0; i < no_of_nodes; i++)
					bla[node_arrives[i]]++;
				cout << "no of arrivals: ";
				for (unsigned short int i = 0; i < periods; i++)
					cout << bla[i] << " ";
				cout << endl << endl;
			}


			ListGraph::EdgeMap<double> edge_card_weight(g, 0);
			ListDigraph::ArcMap<unsigned short int> arc_card_weight(g_original, 0);
			t0 = cpuTime();
			undi_lemon(m, arc_in, arc_out, label_positions, g, g_original, c, c_b, c_original, edge_card_weight, arc_card_weight, no_of_nodes);
			t1 = cpuTime();
			graph_building += t1 - t0;

			unsigned int S = pow(2, N) - 2;
			vector<double> v(S + 1, 0);
			vector<double> v_bar(S + 1, 0);
			vector<double> v_accum(S + 1, 0);
			vector<unsigned short int> s(N, 0);
			vector<double> s_accum(N, 0);
			double prec = pow(10, -7);
			vector<double> init_alloc(N, 0);
			vector<double> init_alloc_accum(N, 0);
			double game_time = 0;
			double init_alloc_time = 0;
			double matching_time = 0;
			vector<double> credit(N, 0);
			vector<double> deviation(N, 0);
			vector<bool> pos(N, false);
			vector<unsigned short int> w(N, 0);
			unsigned short int p;
			vector<unsigned short int> lb(N, 0);
			vector<unsigned short int> ub(N, 0);
			vector<pair<int, int>> arc_pair;
			vector<int> nodeset(graph_size, 0);
			vector<pair<int, int>> cycle_distri;
			vector<double> d(N, 0);
			double d_total = 0;
			double d_c_total = 0;
			double opt = 0;
			vector<bool> leaving(no_of_nodes, false);
			unsigned short int Q = 0;
			double max_d = 0;
			period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit);

			bool unique_impu = false;
			double y_core_dist = 0;
			double s_core_dist = 0;
			//-------------------arbitrary---------------------------
			if (arbitrary) {
				cout << "start arbitrary" << N << "countries" << " " << "instance_" << inst << endl;
				y_core_dist = 0;
				s_core_dist = 0;
				arbitray_maximum = true;
				period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit);
				t0 = cpuTime();
				arbitraryMaximum(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, credit, edge_card_weight, t0, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_arbitrary, numofMaxSolution, arbitray_maximum, initialSize, S, core_100, negative_core, prec, init_alloc, inst, unique_impu, unique_imputation, core_dist, init_alloc_accum, s_accum, v_accum, y_core_dist, s_core_dist, credit_adjusted, M, max_d, game_generation_arbitrary, solution_concept_time_arbitrary, cycle_dis_arbitrary_period);
				t1 = cpuTime();
				total_time_arbitrary += t1 - t0;
				max_d1_arbitrary += max_d;
				arbitray_maximum = false;
				M_100_d_arbitrary += M_total;
				relative_d1_arbitrary += (d_total / M_total);
				cout << N << "countries" << " " << "instance_" << inst << "arbitrary done...";
			}
			//------------------d1----------------
			if (d1) {
				cout << "start minimizing d_1" << N << "countries" << " " << "instance_" << inst << endl;
				period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit);
				t0 = cpuTime();
				min_d_1(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, credit, edge_card_weight, t0, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_d, numofMaxSolution, arbitray_maximum, initialSize, S, core_100, negative_core, prec, init_alloc, d_c_total, inst, unique_impu, unique_imputation, core_dist, init_alloc_accum, s_accum, v_accum, y_core_dist, s_core_dist, credit_adjusted, M, lex_min, max_d, game_generation_d1, solution_concept_time_d1, time_d1, cycle_dis_d_period, cycle_dis_arbitrary_period);
				t1 = cpuTime();
				total_time_d1 += t1 - t0;
				max_d1 += max_d;
				cout << "negative core: " << negative_core;
				relative_d1 += (d_total / M_total);
				cout << "relative_d1: " << relative_d1 << endl;
				cout << "the number of countries: " << N << " " << "relative_d1" << " " << inst << " " << relative_d1 / (inst + 1) << endl;
				M_100 += M_total;
				cout << "the number of countries: " << N << " " << "relative_d1" << " " << inst << " " << M_100 / (inst + 1);
				cout << N << "countries" << " " << "instance_" << inst << "d1 done...";
			}
			//--------------------d1+c--------------------------
			if (d_c) {
				cout << "start minimizing d1_c" << N << "countries" << " " << "instance_" << inst << endl;
				y_core_dist = 0;
				s_core_dist = 0;
				c_involved = true;
				period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit);
				t0 = cpuTime();
				min_d_1(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, credit, edge_card_weight, t0, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_t_c, numofMaxSolution, arbitray_maximum, initialSize, S, core_100, negative_core, prec, init_alloc, d_c_total, inst, unique_impu, unique_imputation, core_dist, init_alloc_accum, s_accum, v_accum, y_core_dist, s_core_dist, credit_adjusted, M, lex_min, max_d, game_generation_d1_c, solution_concept_time_d1_c, time_d1_c, cycle_dis_t_c_period, cycle_dis_arbitrary_period);
				t1 = cpuTime();
				total_time_d1_c += t1 - t0;
				max_d1_c += max_d;
				c_involved = false;
				M_100_d_c += M_total;
				//core_d_c += core_100;
				//negative_core_d_c += negative_core;
				relative_d1_c += (d_total / M_total);
				cout << N << "countries" << " " << "instance_" << inst << "d1+c done... ";
			}
			//--------------------lexmin-----------------------
			if (lexmin_call) {
				cout << "start minimizing lexmin" << N << "countries" << " " << "instance_" << inst << endl;
				y_core_dist = 0;
				s_core_dist = 0;
				lex_min = true;
				period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit);
				t0 = cpuTime();
				min_d_1(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, credit, edge_card_weight, t0, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_lexmin, numofMaxSolution, arbitray_maximum, initialSize, S, core_100, negative_core, prec, init_alloc, d_c_total, inst, unique_impu, unique_imputation, core_dist, init_alloc_accum, s_accum, v_accum, y_core_dist, s_core_dist, credit_adjusted, M, lex_min, max_d, game_generation_lexmin, solution_concept_time_lexmin, time_lex_min, cycle_dis_lexmin_period, cycle_dis_arbitrary_period);
				t1 = cpuTime();
				total_time_lex_min += t1 - t0;
				relative_lexmin_0 += (d_total / M_total);
				lex_min = false;
				M_lex_min += M_total;
				max_lexmin_0 += max_d;
				cout << N << "countries" << " " << "instance_" << inst << "lexmin done..." << endl;
			}
			//--------------------lexmin+c-----------------------
			if (lexmin_c_call) {
				cout << "start minimizing lexmin_c" << N << "countries" << " " << "instance_" << inst << endl;
				y_core_dist = 0;
				s_core_dist = 0;
				c_involved = true;
				lex_min = true;
				period_0(Q, no_of_active_nodes, N, s, Vp, node_arrives, active_nodes, active_nodes_original, c, c_b, c_original, credit);
				t0 = cpuTime();
				min_d_1(node_arrives, g, g_original, arc_pair, leaving, active_nodes, active_nodes_original, c, c_b, c_original, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, credit, edge_card_weight, t0, nodeset, cycle_distri, d, M_total, d_total, c_involved, cycle_dis_lexmin_c, numofMaxSolution, arbitray_maximum, initialSize, S, core_100, negative_core, prec, init_alloc, d_c_total, inst, unique_impu, unique_imputation, core_dist, init_alloc_accum, s_accum, v_accum, y_core_dist, s_core_dist, credit_adjusted, M, lex_min, max_d, game_generation_lexmin_c, solution_concept_time_lexmin_c, time_lex_min_c, cycle_dis_lexmin_c_period, cycle_dis_arbitrary_period);
				t1 = cpuTime();
				total_time_lex_min_c += t1 - t0;
				relative_lexmin_c_0 += (d_total / M_total);
				c_involved = false;
				lex_min = false;
				M_lex_min_c += M_total;
				max_lexmin_c_0 += max_d;
				cout << N << "countries" << " " << "instance_" << inst << "lexmin_c done..." << endl;
			}
		}
		relative_d1_N[N - 4] = relative_d1 / 100;
		relative_d1_N_c[N - 4] = relative_d1_c / 100;
		relative_arbitrary_N[N - 4] = relative_d1_arbitrary / 100;
		relative_lexmin[N - 4] = relative_lexmin_0 / 100;
		relative_lexmin_c[N - 4] = relative_lexmin_c_0 / 100;
		M_N[N - 4] = M_100 / 100;
		M_N_d_c[N - 4] = M_100_d_c / 100;
		M_N_d_arbitrary[N - 4] = M_100_d_arbitrary / 100;
		M_N_lex_min[N - 4] = M_lex_min / 100;
		M_N_lex_min_c[N - 4] = M_lex_min_c / 100;
		max_d1_N[N - 4] = max_d1 / 100;
		max_d1_N_c[N - 4] = max_d1_c / 100;
		max_arbitrary_N[N - 4] = max_d1_arbitrary / 100;
		max_lexmin[N - 4] = max_lexmin_0 / 100;
		max_lexmin_c[N - 4] = max_lexmin_c_0 / 100;
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
		res.open(version + "/" + "results_Banzhaf.txt", ofstream::out | ofstream::trunc);
		for (unsigned short int i = 0; i < N - 3; i++) {
			res << "number of countries: " << N << endl;
			res << "minimizing arbitrary: " << relative_arbitrary_N[i] << endl;
			res << "minimizing d1: " << relative_d1_N[i] << endl;
			res << "minimizing d1+c: " << relative_d1_N_c[i] << endl;
			res << "minimizing lexmin: " << relative_lexmin[i] << endl;
			res << "minimizing lexmin+c: " << relative_lexmin_c[i] << endl;
			res << "max arbitrary: " << max_arbitrary_N[i] << endl;
			res << "max d1: " << max_d1_N[i] << endl;
			res << "max d1+c: " << max_d1_N_c[i] << endl;
			res << "max lexmin: " << max_lexmin[i] << endl;
			res << "max lexmin+c: " << max_lexmin_c[i] << endl;
			res << "averange number of transplants arbitrary: " << M_N_d_arbitrary[i] << endl;
			res << "averange number of transplants d1: " << M_N[i] << endl;
			res << "averange number of transplants d1+c: " << M_N_d_c[i] << endl;
			res << "averange number of transplants lexmin: " << M_N_lex_min[i] << endl;
			res << "averange number of transplants lexmin+c: " << M_N_lex_min_c[i] << endl;
			res << "data preparation: " << data_preparation_N[i] << endl;
			res << "build graph: " << graph_building_N[i] << endl;
			res << "total time arbitrary: " << total_time_arbitrary_N[i] << endl;
			res << "total time d1: " << total_time_d1_N[i] << endl;
			res << "total time d1+c: " << total_time_d1_c_N[i] << endl;
			res << "total time lexmin: " << total_time_lex_min_N[i] << endl;
			res << "total time lexmin+c: " << total_time_lex_min_c_N[i] << endl;
			res << "scenario time arbitrary: " << time_arbitrary_N[i] << endl;
			res << "scenario time d1: " << time_d1_N[i] << endl;
			res << "scenario time d1+c: " << time_d1_c_N[i] << endl;
			res << "scenario time lexmin: " << time_lex_min_N[i] << endl;
			res << "scenario time lexmin+c: " << time_lex_min_c_N[i] << endl;
			res << "game generation arbitrary: " << game_generation_arbitrary_N[i] << endl;
			res << "game generation d1: " << game_generation_d1_N[i] << endl;
			res << "game generation d1+c: " << game_generation_d1_c_N[i] << endl;
			res << "game generation lexmin: " << game_generation_lexmin_N[i] << endl;
			res << "game generation lexmin+c: " << game_generation_lexmin_c_N[i] << endl;
			res << "solution concept arbitrary: " << solution_concept_time_arbitrary_N[i] << endl;
			res << "solution concept d1: " << solution_concept_time_d1_N[i] << endl;
			res << "solution concept d1+c: " << solution_concept_time_d1_c_N[i] << endl;
			res << "solution concept lexmin: " << solution_concept_time_lexmin_N[i] << endl;
			res << "solution concept lexmin+c: " << solution_concept_time_lexmin_c_N[i] << endl;
			res << endl;
		}

		res.close();

		vector<long> check(5, 0);
		if (d1) {
			ofstream res_dis;
			res_dis.open(version + "/" + "cycle_dis_d" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
			for (const auto& elem : cycle_dis_d) {
				res_dis << elem.first << ": " << elem.second << endl;
				check[0] += elem.first * elem.second;
			}
			res_dis << endl;
			res_dis.close();
		}


		if (d_c) {
			ofstream res_dis_c;
			res_dis_c.open(version + "/" + "cycle_dis_c" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
			for (const auto& elem : cycle_dis_t_c) {
				res_dis_c << elem.first << ": " << elem.second << endl;
				check[1] += elem.first * elem.second;
			}
			res_dis_c << endl;
			res_dis_c.close();
		}

		if (arbitrary) {
			ofstream res_dis_arbitrary;
			res_dis_arbitrary.open(version + "/" + "cycle_dis_arbitrary" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
			for (const auto& elem : cycle_dis_arbitrary) {
				res_dis_arbitrary << elem.first << ": " << elem.second << endl;
				check[2] += elem.first * elem.second;
			}
			res_dis_arbitrary << endl;
			res_dis_arbitrary.close();
		}

		if (lexmin_call) {
			ofstream res_dis_lexmin;
			res_dis_lexmin.open(version + "/" + "cycle_dis_lexmin" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
			for (const auto& elem : cycle_dis_lexmin) {
				res_dis_lexmin << elem.first << ": " << elem.second << endl;
				check[3] += elem.first * elem.second;
			}
			res_dis_lexmin << endl;
			res_dis_lexmin.close();
		}

		if (lexmin_c_call) {
			ofstream res_dis_lexmin_c;
			res_dis_lexmin_c.open(version + "/" + "cycle_dis_lexmin_c" + to_string(N) + ".txt", ofstream::out | ofstream::trunc);
			for (const auto& elem : cycle_dis_lexmin_c) {
				res_dis_lexmin_c << elem.first << ": " << elem.second << endl;
				check[4] += elem.first * elem.second;
			}
			res_dis_lexmin_c << endl;
			res_dis_lexmin_c.close();
		}



		//cycle distribution based on periods
		vector<long> value(5, 0);
		if (d1) {
			ofstream res_cycle_dis_d_period;
			for (unsigned short int i = 0; i < 24; ++i) {
				res_cycle_dis_d_period.open(version + "/" + "cycle_dis_d_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_d_period[(N) * (i + 1)]) {
					value[0] += elem.first * elem.second;
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
				res_cycle_dis_c_period.open(version + "/" + "cycle_dis_c_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_t_c_period[(N) * (i + 1)]) {
					value[1] += elem.first * elem.second;
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
				res_cycle_dis_arbitrary_period.open(version + "/" + "cycle_dis_arbitrary_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_arbitrary_period[(N) * (i + 1)]) {
					value[2] += elem.first * elem.second;
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
				res_dis_lexmin_period.open(version + "/" + "cycle_dis_lexmin_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_lexmin_period[(N) * (i + 1)]) {
					value[3] += elem.first * elem.second;
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
				res_cycle_dis_lexmin_c_period.open(version + "/" + "cycle_dis_lexmin_c_period" + to_string(N) + "_" + to_string(i) + ".txt", ofstream::out | ofstream::trunc);
				for (const auto& elem : cycle_dis_lexmin_c_period[(N) * (i + 1)]) {
					value[4] += elem.first * elem.second;
					res_cycle_dis_lexmin_c_period << elem.first << ": " << elem.second << endl;
				}
				res_cycle_dis_lexmin_c_period << endl;
				res_cycle_dis_lexmin_c_period.close();
			}
			res_cycle_dis_lexmin_c_period.close();
		}

		for (unsigned short int i = 0; i < 5; ++i) {
			cout << "scenario: " << i << " " << value[i] << " " << check[i] << endl;
			if (value[i] != check[i]) {
				cout << "scenario: " << i << " " << value[i] << " " << check[i] << endl;
				cout << "Error in the number of transplants" << endl;
			}
		}


	}
	return 0;
}



void coop_game(ListGraph& g, vector<double>& v, unsigned int& S, vector<unsigned short int>& s, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, ListGraph::EdgeMap<double>& edge_card_weight, bool& dispy, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, map<int, int>& numofMaxSolution, unsigned short int& Q, bool& arbitray_maximum, vector<pair<int, int>>& cycle_distri, map<int, int>& cycle_dis, bool& unique_impu, double& prec, bool& unique_imputation, bool& credit_adjusted, vector<double>& credit, double& M, double& game_generation, std::map<int, std::map<int, int >>& cycle_dis_arbitrary_period) {
	vector<bool> a(N, false);
	vector<double> v_copy(S + 1, 0);
	double t0 = cpuTime();
	for (unsigned int i = 0; i < S; i++) {
		v_copy[i] = 0;
		de2bi(i, a, N);
		ListGraph::NodeMap<bool> coal3(g, false);
		for (unsigned short int j = 0; j < N; j++) {
			if (a[j]) {
				if (credit_adjusted) {
					v_copy[i] += credit[j];
				}
				for (unsigned short int k = j * Vp; k < (j + 1) * Vp; ++k) {
					if (active_nodes[c[k]]) {
						coal3[c[k]] = true;
						coal3[c_b[k]] = true;
					}
				}
			}

		}
		MaxWeightedPerfectMatching<FilterNodes<ListGraph>, ListGraph::EdgeMap<double>> coal_m1(FilterNodes<ListGraph>(g, coal3), edge_card_weight);
		coal_m1.run();
		v_copy[i] += coal_m1.matchingWeight();
		v[i] = v_copy[i];
		cout << "v[i]: " << v[i] << "v_copy[i]: " << v_copy[i] << endl;
	}
	cout << "finish generating the copy" << endl;

	MaxWeightedPerfectMatching<FilterNodes<ListGraph>, ListGraph::EdgeMap<double>> grand_coal(FilterNodes<ListGraph>(g, active_nodes), edge_card_weight);
	cout << "finish generating v[N]_STEP_1" << endl;
	grand_coal.run();
	cout << "finish generating v[N]_STEP_2" << endl;
	//grand_coal.matchingWeight();
	//cout << "finish generating v[N]" << endl;
	double t1 = cpuTime();
	game_generation += t1 - t0;
	grand_coal.matchingMap();

	FilterNodes<ListGraph> sg(g, active_nodes);
	if (arbitray_maximum) {
		for (FilterNodes<ListGraph>::NodeIt n(FilterNodes<ListGraph>(g, active_nodes)); n != INVALID; ++n) {
			//cout << FilterNodes<ListGraph>(g, active_nodes).id(n) << endl;
			if (!(grand_coal.matching(n) == INVALID) && edge_card_weight[grand_coal.matching(n)] > 0 && FilterNodes<ListGraph>(g, active_nodes).id(n) % 2 == 0) {
				cycle_distri.push_back({ FilterNodes<ListGraph>(g, active_nodes).id(n) / 2, (FilterNodes<ListGraph>(g, active_nodes).id(grand_coal.mate(n)) - 1) / 2 });
				leaving[FilterNodes<ListGraph>(g, active_nodes).id(n) / 2] = true;
				//cout << "FilterNodes<ListGraph>(g, active_nodes).id(n):"<<FilterNodes<ListGraph>(g, active_nodes).id(n) << endl;
				//cout << "FilterNodes<ListGraph>(g, active_nodes).id(n): " << FilterNodes<ListGraph>(g, active_nodes).id(n) << endl;
				for (unsigned short int i = 0; i < N; ++i) {
					if (i * Vp <= FilterNodes<ListGraph>(g, active_nodes).id(n) / 2 && FilterNodes<ListGraph>(g, active_nodes).id(n) / 2 < (i + 1) * Vp) {
						++s[i];
					}
				}

			}
			//cycle_distri.push_back(arc_pair[i - 1]);
		}

		cycle_distribution(cycle_dis_arbitrary_period, cycle_dis, cycle_distri, N, Q);
		for (unsigned short int i = 0; i < N; ++i) {
			v_copy[S] += s[i];
		}
		v[S] = v_copy[S];
		M = v[S];
	}
	else {
		for (FilterNodes<ListGraph>::NodeIt n(FilterNodes<ListGraph>(g, active_nodes)); n != INVALID; ++n) {
			cout << FilterNodes<ListGraph>(g, active_nodes).id(n) << endl;
			if (!(grand_coal.matching(n) == INVALID) && edge_card_weight[grand_coal.matching(n)] > 0 && FilterNodes<ListGraph>(g, active_nodes).id(n) % 2 == 0) {
				++v_copy[S];
			}
		}
		v[S] = v_copy[S];
		M = v[S];
	}
	

	if (unique_imputation) {
		double surplus = v[S];
		for (unsigned short int i = 0; i < N; i++)
			surplus -= v[pow(2, i) - 1];
		if (surplus > prec)
			unique_impu = false;
		else
			unique_impu = true;
	}
	if (dispy)
		cout << "grand coal: " << v[S] << endl;

	if (dispy) {
		cout << "s: ";
		for (unsigned short int i = 0; i < N; i++) {
			cout << s[i] << " ";
		}
		cout << endl;
	}
	return;
}

void norm_banzhaf(vector<double>& banz, vector<double>& v, unsigned short int& n, unsigned int& s) {
	double w = 1 / pow(2, n - 1);
	vector<double> expo(n, 1);
	for (unsigned short int j = 1; j < n; j++)
		expo[j] = pow(2, j);
	vector<bool> a(n, false);
	unsigned short int k = 0;
	for (unsigned short int i = 0; i < n; i++)
		banz[i] = w * (double)(v[expo[i] - 1]);
	for (unsigned int i = 0; i < s; i++) {
		de2bi_card(i, a, n, k);
		for (unsigned short int j = 0; j < n; j++) {
			if (!a[j])
				banz[j] += w * (v[i + expo[j]] - v[i]);
		}
	}
	double norm = banz[0];
	for (unsigned short int i = 1; i < n; i++)
		norm += banz[i];
	for (unsigned short int i = 0; i < n; i++)
		banz[i] = banz[i] * v[s] / norm;
	return;
}

void de2bi_card(unsigned int& k, vector<bool>& a, unsigned short int& n, unsigned short int& card) {
	vector<bool> zero(n, false);
	card = 0;
	a = zero;
	unsigned int i = 2;
	for (unsigned short int c = 0; c < n - 2; c++)
		i += i;
	unsigned int j = k + 1;
	unsigned short int l = n - 1;
	while (j > 0) {
		if (j >= i) {
			a[l] = true;
			card++;
			j -= i;
		}
		i /= 2;
		l--;
	}
	return;
}

void de2bi(unsigned int& k, vector<bool>& a, unsigned short int& n) {
	vector<bool> zero(n, false);
	a = zero;
	unsigned int i = 2;
	for (unsigned short int c = 0; c < n - 2; c++)
		i += i;
	unsigned int j = k + 1;
	unsigned short int l = n - 1;
	while (j > 0) {
		if (j >= i) {
			a[l] = true;
			j -= i;
		}
		i /= 2;
		l--;
	}
	return;
}

double core_distance(vector<double>& x, vector<double>& v, unsigned short int& n, unsigned int& S) {
	double eps = x[0] - v[0];
	if (x[1] - v[1] < eps)
		eps = x[1] - v[1];
	vector<bool> a(n, false);
	double xS = 0;
	for (unsigned int i = 2; i < S; i++) {
		de2bi(i, a, n);
		for (unsigned short int j = 0; j < n; j++)
			if (a[j])
				xS += x[j];
		if (xS - v[i] < eps)
			eps = xS - v[i];
		xS = 0;
	}
	return eps;
}


void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, vector<unsigned short int>& s, vector<double>& d) {
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
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

void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, unsigned short int& initialSize) {
	unsigned short int coal = rand() % Vp;
	unsigned short int count = 0;
	if (disp) {
		cout << "Initial pairs:" << endl;
	}
	for (unsigned short int i = 0; i < N; i++) {
		if (disp) {
			cout << "Country " << i << " : ";
		}
		while (count < initialSize) {
			if (active_nodes[c[i * Vp + coal]]) {
				coal = rand() % Vp;
			}
			else {
				active_nodes[c[i * Vp + coal]] = true;
				count++;
				if (disp) {
					cout << i * Vp + coal << " (" << count << ") ";
				}
				coal = rand() % Vp;
			}
		}
		if (disp) {
			cout << endl;
		}
		count = 0;
	}
	return;
}

void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods, bool& disp) {
	for (unsigned short int i = 0; i < N; i++) {
		if (disp)
			cout << "Country " << i << " arrivals: ";
		for (unsigned short int j = 0; j < Vp; j++) {
			if (!(active_nodes[c[i * Vp + j]])) {
				node_arrives[i * Vp + j] = rand() % (periods - 1) + 1;
			}
			if (disp)
				cout << node_arrives[i * Vp + j] << " ";
		}
		if (disp)
			cout << endl;
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
					cout << "ID ERROR!" << endl;
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					if (node_labels[n - 1] != 10 * char2uint(line[l]) + char2uint(line[l + 1]))
						cout << "ID ERROR!" << endl;
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						if (node_labels[n - 1] != 100 * char2uint(line[l]) + 10 * char2uint(line[l + 1]) + char2uint(line[l + 2]))
							cout << "ID ERROR!" << endl;
						l = l + 2;
					}
					else {
						//if (!is_next_char_digit(line, l + 3)){
						if (node_labels[n - 1] != 1000 * char2uint(line[l]) + 100 * char2uint(line[l + 1]) + 10 * char2uint(line[l + 2]) + char2uint(line[l + 3]))
							cout << "ID ERROR!" << endl;
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

void arbitraryMaximum(vector<unsigned short int>& node_arrives, ListGraph& g, ListDigraph& g_original, vector<pair<int, int>>& arc_pair, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<double>& v, vector<double>& credit, ListGraph::EdgeMap<double>& edge_card_weight, double& t0, vector<int>& nodeset, vector<pair<int, int>>& cycle_distri, vector<double>& d, double& M_total, double& d_total, bool& c_involved, map<int, int>& cycle_dis, map<int, int>& numofMaxSolution, bool& arbitray_maximum, unsigned short int& initialSize, unsigned int& S, double& core_100, long& negative_core, double& prec, vector<double>& init_alloc, unsigned short int inst, bool& unique_impu, bool& unique_imputation, bool& core_dist, vector<double>& init_alloc_accum, vector<double>& s_accum, vector<double>& v_accum, double& y_core_dist, double& s_core_dist, bool& credit_adjusted, double& M, double& max_d, double& game_generation, double& solution_concept_time, std::map<int, std::map<int, int>>& cycle_dis_arbitrary_period)
{
	Q = 0;
	d_total = 0;
	M_total = 0;
	core_100 = 0;
	negative_core = 0;
	y_core_dist = 0;
	s_core_dist = 0;
	max_d = 0;
	if (dispy)
		cout << " --== Without lex min matching == -- " << endl;
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
			cout << "--== PERIOD " << Q + 1 << " ==--" << endl;
		}
		if (dispy) {
			cout << "Number of active nodes: ";
			for (unsigned short int i = 0; i < N; i++)
				cout << no_of_active_nodes[i] << " ";
			cout << endl;
		}
		// cooperative game and target
		cout << "start generating values" << endl;
		coop_game(g, v, S, s, c, c_b, edge_card_weight, dispy, Vp, N, active_nodes, leaving, numofMaxSolution, Q, arbitray_maximum, cycle_distri, cycle_dis, unique_impu, prec, unique_imputation, credit_adjusted, credit, M, game_generation, cycle_dis_arbitrary_period);
		double t0 = cpuTime();
		norm_banzhaf(init_alloc, v, N, S);
		double t1 = cpuTime();
		solution_concept_time += t1 - t0;
		M_total += M;



		for (unsigned short int i = 0; i < N; ++i) {
			if (credit_adjusted) {
				d[i] += init_alloc[i] - credit[i] - s[i];
				credit[i] = init_alloc[i] - s[i];

				cout << "credit[i]" << credit[i] << endl;
			}
			else {
				if (c_involved) {
					credit[i] += init_alloc[i] - s[i];
					d[i] += init_alloc[i] - s[i];
				}
				else {
					credit[i] = 0;
					d[i] += init_alloc[i] - s[i];
				}

			}
			cout << "country" << to_string(i) << "init_alloc[i]: " << init_alloc[i] << '/n' << "s[i]: " << s[i] << "d[i]: " << d[i] << "credit[i]: " << credit[i] << endl;
			//actual_alloc[Q].push_back(s[i]);
		}


		Q++;
		changing_nodes(active_nodes, active_nodes_original, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, c_b, c_original, s, d);
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

void period_0(unsigned short int& Q, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, vector<unsigned short int>& s, unsigned short int& Vp, vector<unsigned short int>& node_arrives, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, vector<double>& credit) {
	Q = 0;
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		credit[i] = 0;
		no_of_active_nodes[i] = Vp / 4;
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
	return;
}

void min_d_1(vector<unsigned short int>& node_arrives, ListGraph& g, ListDigraph& g_original, vector<pair<int, int>>& arc_pair, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<ListGraph::Node>& c, vector<ListGraph::Node>& c_b, vector<ListDigraph::Node>& c_original, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<double>& v, vector<double>& credit, ListGraph::EdgeMap<double>& edge_card_weight, double& t0, vector<int>& nodeset, vector<pair<int, int>>& cycle_distri, vector<double>& d, double& M_total, double& d_total, bool& c_involved, map<int, int>& cycle_dis, map<int, int>& numofMaxSolution, bool& arbitray_maximum, unsigned short int& initialSize, unsigned int& S, double& core_100, long& negative_core, double& prec, vector<double>& init_alloc, double& d_c_total, unsigned short int inst, bool& unique_impu, bool& unique_imputation, bool& core_dist, vector<double>& init_alloc_accum, vector<double>& s_accum, vector<double>& v_accum, double& y_core_dist, double& s_core_dist, bool& credit_adjusted, double& M, bool lex_min, double& max_d, double& game_generation, double& solution_concept_time, double& scenario_time, std::map<int, std::map<int, int>>& cycle_dis_period, std::map<int, std::map<int, int>>& cycle_dis_arbitrary_period) {
	Q = 0;
	d_total = 0;
	d_c_total = 0;
	M_total = 0;
	max_d = 0;
	core_100 = 0;
	negative_core = 0;
	y_core_dist = 0;
	s_core_dist = 0;
	if (dispy)
		cout << " --== Without lex min matching == -- " << endl;
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
			cout << "--== PERIOD " << Q + 1 << " ==--" << endl;
		}
		if (dispy) {
			cout << "Number of active nodes: ";
			for (unsigned short int i = 0; i < N; i++)
				cout << no_of_active_nodes[i] << " ";
			cout << endl;
		}
		// cooperative game and target
		cout << "start generating values" << endl;
		coop_game(g, v, S, s, c, c_b, edge_card_weight, dispy, Vp, N, active_nodes, leaving, numofMaxSolution, Q, arbitray_maximum, cycle_distri, cycle_dis, unique_impu, prec, unique_imputation, credit_adjusted, credit, M, game_generation, cycle_dis_arbitrary_period);
		double t0 = cpuTime();
		norm_banzhaf(init_alloc, v, N, S);
		double t1 = cpuTime();
		solution_concept_time += t1 - t0;

		//init_alloc[Q] = target;
		t0 = cpuTime();
		ILP_d1_gurobi(Q, N, g_original, Vp, node_arrives, active_nodes, active_nodes_original, arc_pair, nodeset, M, M_total, s, cycle_distri, leaving, d, d_total, c_involved, credit, cycle_dis, init_alloc, credit_adjusted, lex_min, inst, cycle_dis_period);
		t1 = cpuTime();
		scenario_time += t1 - t0;
		Q++;
		changing_nodes(active_nodes, active_nodes_original, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, c_b, c_original, s, d);
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

void ILP_d1_gurobi(unsigned short int& Q, unsigned short int& N, ListDigraph& g_original, unsigned short int& Vp, vector<unsigned short int>& node_arrives, ListGraph::NodeMap<bool>& active_nodes, ListDigraph::NodeMap<bool>& active_nodes_original, vector<pair<int, int>>& arc_pair, vector<int>& nodeset, double& M, double& M_total, vector<unsigned short int>& s, vector<pair<int, int>>& cycle_distri, vector<bool>& leaving, vector<double>& d, double& d_total, bool& c_involved, vector<double>& credit, map<int, int>& cycle_dis, vector<double>& init_alloc, bool& credit_adjusted, bool lex_min, unsigned short int inst, std::map<int, std::map<int, int>>& cycle_dis_period) {
	pair_arcs(Q, g_original, node_arrives, active_nodes_original, arc_pair, nodeset);
	FilterNodes<ListDigraph> sg(g_original, active_nodes_original);
	long col_num = countArcs(sg) + 1;
	cout << "column: " << col_num - 1 << endl;
	cout << "countArcs(g_original): " << countArcs(g_original) << endl;
	M_total += M;
	//int node_number = countNodes(g_original);
	const unsigned short int row_num = N + 2 * nodeset.size() + 1;
	// Create an environment
	GRBEnv env = GRBEnv(true);
	env.set("LogFile", "mip1_banzhaf_d_c_star.log");
	env.start();

	// Create an empty model
	GRBModel model = GRBModel(env);
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
		if (credit_adjusted) {
			bound[i - 1] = init_alloc[i - 1];
		}
		else {
			bound[i - 1] = init_alloc[i - 1] + credit[i - 1];
		}
	}
	/*for (unsigned short int i = 1 + N; i < 2 * N + 1; ++i) {
		bound[i - 1] = target[i - 1];
	}*/
	for (unsigned short int i = N + 1; i < N + 1 + nodeset.size(); ++i) {
		bound[i - 1] = 0;
	}
	for (unsigned short int i = N + nodeset.size() + 1; i < N + 2 * (nodeset.size()) + 1; ++i) {
		bound[i - 1] = 1;
	}
	bound[row_num - 1] = M;
	cout << "size of max-weighted perfect matching: " << M << endl;
	long matrix_num = row_num * col_num;
	vector<int> ia(matrix_num + 1, 0);
	//cout << "ia.size(): " << ia.size();
	vector<int> ja(matrix_num + 1, 0);
	vector<double> ar(matrix_num + 1, 0);

	//cout << "Period" << to_string(Q) << endl;
	//cout << "N: " << N << "col_num: " << col_num << "N*col_num: " << N * col_num;
	long cnt_2 = 0;


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
	cout << "finish loading efficiencies" << endl;
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
		epsilon_func(init_alloc, credit, epsilon, N, credit_adjusted);
		sort_d_t(d_t, var_bi, col_num, N, Vp, arc_pair, init_alloc, t_star, credit, epsilon, var_lexmin, N_star, credit_adjusted);
	}
	std::cout << "finish sorting" << "epsilon:" << epsilon << endl;
	std::cout << "start n_star_1" << endl;
	if (lex_min && d_t[0] > 0.5 && abs(epsilon) > pow(10, -4)) {
		lex_min_n_star(d_t, lex_min, t_star, N, col_num, epsilon, n_star, model, ia, ja, ar, row_num, cnt_2, bound, nodeset, N_star, var_lexmin, vector_model, track);
	}
	std::cout << "finish n_star_1" << endl;
	lexmin_searching(d_t, lex_min, t_star, N, col_num, epsilon, n_star, model, ia, ja, ar, row_num, cnt_2, bound, nodeset, N_star, Vp, arc_pair, init_alloc, credit, var_lexmin, inst, vector_model, track, credit_adjusted);

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
	cout << "maximum size: " << t;
	for (unsigned short int i = 0; i < N; ++i) {
		if (credit_adjusted) {
			d[i] += init_alloc[i] - credit[i] - s[i];
			credit[i] = init_alloc[i] - s[i];
		}
		else {
			if (c_involved) {
				credit[i] += init_alloc[i] - s[i];
				d[i] += init_alloc[i] - s[i];
			}
			else {
				credit[i] = 0;
				d[i] += init_alloc[i] - s[i];
			}

		}
		cout << "country" << to_string(i) << "init_alloc[i]: " << init_alloc[i] << '/n' << "s[i]: " << s[i] << "d[i]: " << d[i] << "credit[i]: " << credit[i] << endl;
		//actual_alloc[Q].push_back(s[i]);
	}
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

void sort_d_t(vector<double>& d_t, vector<GRBVar>& var_bi, long& col_num, unsigned short int& N, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, unsigned short int& t, vector<double>& credit, double& epsilon, vector<GRBVar>& var_lexmin, vector<unsigned short int>& N_star, bool& credit_adjusted) {
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
		if (credit_adjusted) {
			d_copy[i] = abs(target[i] - s_copy[i]);
		}
		else {
			d_copy[i] = abs(target[i] + credit[i] - s_copy[i]);
		}

	}
	//sort
	std::sort(d_copy.begin(), d_copy.end());
	for (unsigned short int i = 0; i < N; ++i) {
		std::cout << "d_copy[i]" << d_copy[i] << endl;
		//std::cout << "s_copy" << s_copy[i] << endl;
	}
	//epsilon_func(target, credit, epsilon, N);
	if (epsilon > 3 * pow(10, -4)) {
		d_t[t] = d_copy[N - 1 - t] + pow(10, -4);
	}
	else {
		d_t[t] = d_copy[N - 1 - t];
	}
	if (t < N - 1) {
		d_t[t + 1] = d_copy[N - 2 - t];
	}
	std::cout << "d" << to_string(t - 1) << d_t[t - 1] << endl;
	std::cout << "epsilon: " << epsilon << endl;
}

void epsilon_func(vector<double>& target, vector<double>& credit, double& epsilon, unsigned short int N, bool& credit_adjusted) {
	vector<double> target_credit(N, 0);
	vector<double> epsilon_sort(N * (N - 1), 0);
	for (unsigned short int i = 0; i < N; ++i) {
		//cout << "target[i]: " << target[i] << endl;
		if (credit_adjusted) {
			target_credit[i] = target[i];
		}
		else {
			target_credit[i] = target[i] + credit[i];
		}
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

	/*for (unsigned short int i = 0; i < epsilon_sort.size(); ++i) {
		cout << "epsilon_sort: " << epsilon_sort[i] << endl;
	}*/

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

void lex_min_d_star(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, vector<double>& credit, vector<GRBVar>& var_lexmin, vector<GRBModel>& vector_model, unsigned short int& track, bool& credit_adjusted) {
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
	sort_d_t(d_t, var_bi, col_num, N, Vp, arc_pair, target, t_star, credit, epsilon, var_lexmin, N_star, credit_adjusted);
}

void lexmin_searching(vector<double>& d_t, bool& lex_min, unsigned short int& t_star, unsigned short int& N, long& col_num, double& epsilon, unsigned short int& n_star, GRBModel& model, vector<int>& ia, vector<int>& ja, vector<double>& ar, const unsigned short int& row_num, long& cnt_2, vector<double>& bound, vector<int>& nodeset, vector<unsigned short int>& N_star, unsigned short int& Vp, vector<pair<int, int>>& arc_pair, vector<double>& target, vector<double>& credit, vector<GRBVar>& var_lexmin, unsigned short int inst, vector<GRBModel>& vector_model, unsigned short int& track, bool& credit_adjusted) {
	while (lex_min && abs(d_t[t_star - 1]) > 0.5 && n_star < N) {
		std::cout << "begin search d_t" << endl;
		lex_min_d_star(d_t, lex_min, t_star, N, col_num, epsilon, n_star, model, ia, ja, ar, row_num, cnt_2, bound, nodeset, N_star, Vp, arc_pair, target, credit, var_lexmin, vector_model, track, credit_adjusted);
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
				lexmin_searching(d_t, lex_min, t_star, N, col_num, epsilon, n_star, model, ia, ja, ar, row_num, cnt_2, bound, nodeset, N_star, Vp, arc_pair, target, credit, var_lexmin, inst, vector_model, track, credit_adjusted);
			}
		}
		else {
			std::cout << "congratulations break" << endl;
			//break
			n_star = N;
		}
	}
}

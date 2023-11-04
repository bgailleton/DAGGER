#pragma once
#include "dodcontexts.hpp"
#include "graphflood_enums.hpp"
#include "graphflood_parts.hpp"
#include "parambag.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class i_t,
				 class f_t,
				 class CONNECTOR_T,
				 class GRAPH_T,
				 class DBAG_T,
				 class PARAM_T>
class Graphflood2
{
public:
	Graphflood2(){};

	Graphflood2(CONNECTOR_T& con, GRAPH_T& gra, DBAG_T& data, PARAM_T& param)
	{
		this->con = &con;
		this->gra = &gra;
		this->data = &data;
		this->param = &param;
	}

	CONNECTOR_T* con;
	GRAPH_T* gra;
	DBAG_T* data;
	PARAM_T* param;

	WATER_INPUT water_input_mode = WATER_INPUT::PRECIPITATIONS_CONSTANT;
	f_t Prate = 1e-5; // mm.s^-1. Yarr.
	void set_uniform_P(f_t val) { this->Prate = val; }

	std::vector<i_t> entry_node_PQ;
	std::vector<i_t> input_node_Qw;
	std::vector<f_t> input_Qw;
	std::vector<f_t> input_Qs;
	std::vector<std::uint8_t> isInQ;

	// Water transfer params
	f_t min_part_Qw = 0.1;

	// time management thingies
	// std::vector<f_t> data->_timetracker; // MOVED TO THE DATA BAG
	f_t time = 0.;
	f_t dt = 1e-3;
	void set_dt(f_t val) { this->dt = val; }
	f_t get_dt() { return this->dt; }

	// Graphflood params
	f_t mannings = 0.033;
	f_t hw_increment_LM = 1e-3;

	bool prefill_lakes = true;

	int debugyolo = 0;

	// std::vector<i_t> active_nodes;

	void init()
	{

		if (this->data == nullptr)
			throw std::runtime_error(
				"Cannot init Graphflood2 -> no data bag connected");
		if (this->con == nullptr)
			throw std::runtime_error(
				"Cannot init Graphflood2 -> no connector connected");
		// if(this->gra == nullptr) throw std::runtime_error("Cannot init
		// Graphflood2 -> no graph connected");

		if (this->data->_surface.size() == 0)
			throw std::runtime_error(
				"databag has no surface data, cannot init GraphFlood2");
		if (this->data->_hw.size() == 0)
			this->data->_hw = std::vector<f_t>(this->con->nxy(), 0);
		if (this->data->_Qwin.size() == 0)
			this->data->_Qwin = std::vector<f_t>(this->con->nxy(), 0);
		if (this->data->_Qwout.size() == 0)
			this->data->_Qwout = std::vector<f_t>(this->con->nxy(), 0);
		this->data->_vmot_hw = std::vector<f_t>(this->con->nxy(), 0);
		this->isInQ = std::vector<std::uint8_t>(this->con->nxy(), false);
	}

	void initial_fill()
	{
		// copying the bedrock elevation
		std::vector<f_t> bedrock(this->data->_surface);

		this->con->reinit();

		// filling and computing the graph from it
		this->con->PFcompute_all();
		if (this->prefill_lakes) {
			for (int i = 0; i < this->con->nxy(); ++i)
				this->data->_hw[i] += this->data->_surface[i] - bedrock[i];
		}

		// backpropagating the bedrock to topo
		this->data->_surface = bedrock;
	}

	void compute_entry_points_from_P(f_t Qw_threshold)
	{
		this->_computeEntryPoints_prec(Qw_threshold);
	}

	void _computeEntryPoints_prec(f_t Qw_threshold)
	{

		this->entry_node_PQ.clear();
		this->input_node_Qw.clear();
		this->input_Qw.clear();

		// copying the bedrock elevation
		std::vector<f_t> bedrock(this->data->_surface);

		this->con->reinit();

		// filling and computing the graph from it
		this->con->PFcompute_all();
		if (this->prefill_lakes) {
			for (int i = 0; i < this->con->nxy(); ++i)
				this->data->_hw[i] += this->data->_surface[i] - bedrock[i];
		}

		// backpropagating the bedrock to topo
		this->data->_surface = bedrock;

		// temporary QWin
		std::vector<f_t> QW(this->con->_nxy, 0.);

		// Calculating QWin
		for (int i = this->con->nxy() - 1; i >= 0; --i) {
			int node = this->data->_Sstack[i];
			if (nodata(this->data->_boundaries[node]))
				;
			int rec = this->con->Sreceivers(node);
			QW[node] += Prate * this->con->area(node);
			if (node != rec) {
				QW[rec] += QW[node];
			}
		}

		std::vector<std::uint8_t> isdone(this->con->nxy(), false);

		fillvec(this->data->_Qwin, 0.);

		for (int i = this->con->nxy() - 1; i >= 0; --i) {
			int node = this->data->_Sstack[i];
			if (nodata(this->data->_boundaries[node]))
				continue;

			int rec = this->con->Sreceivers(node);
			if (QW[rec] >= Qw_threshold && QW[node] < Qw_threshold &&
					isdone[rec] == false) {
				this->input_node_Qw.emplace_back(node);
				this->input_Qw.emplace_back(QW[node]);
				this->entry_node_PQ.emplace_back(node);
				isdone[node] = true;
				while (rec != node) {
					isdone[rec] = true;
					node = rec;
					rec = this->con->Sreceivers(node);
				}
			} else if (QW[node] >= Qw_threshold) {
				this->data->_Qwin[node] += Prate * this->con->area(node);

			} else if (QW[rec] >= Qw_threshold) {
				this->data->_Qwin[rec] += QW[node];
			}
		}
		for (int i = this->con->nxy() - 1; i >= 0; --i) {
			if (this->data->_Qwin[i] > 0) {
				this->input_node_Qw.emplace_back(i);
				this->input_Qw.emplace_back(this->data->_Qwin[i]);
			}
		}
	}

	template<class arrin_i_t, class arrin_f_t>
	void set_Qw_input_points(arrin_i_t& tarri, arrin_f_t& tarrf)
	{
		auto arri = format_input<arrin_i_t>(tarri);
		this->input_node_Qw = to_vec(arri);
		this->entry_node_PQ = to_vec(arri);
		auto arrf = format_input<arrin_f_t>(tarrf);
		this->input_Qw = to_vec(arrf);
		this->input_Qs = std::vector<f_t>(this->input_Qw.size(), 0.);
	}

	template<class arrin_i_t, class arrin_f_t>
	void set_QwQs_input_points(arrin_i_t& tarri,
														 arrin_f_t& tarrf,
														 arrin_f_t& tarrfqs)
	{
		auto arri = format_input<arrin_i_t>(tarri);
		this->input_node_Qw = to_vec(arri);
		this->entry_node_PQ = to_vec(arri);
		auto arrf = format_input<arrin_f_t>(tarrf);
		this->input_Qw = to_vec(arrf);
		auto arrfqs = format_input<arrin_f_t>(tarrfqs);
		this->input_Qs = to_vec(arrfqs);
	}

	template<class CELL, class DSTACK>
	void init_dstack_Water(DSTACK& dstack)
	{

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			this->data->_Qwin[this->input_node_Qw[i]] = this->input_Qw[i];
		}
		for (size_t i = 0; i < this->entry_node_PQ.size(); ++i) {
			dstack.emplace(CELL(this->entry_node_PQ[i],
													this->data->_surface[this->entry_node_PQ[i]]));
			this->isInQ[this->entry_node_PQ[i]] = true;
		}
	}

	template<class CELL, class DSTACK>
	void init_dstack_WaterSed(DSTACK& dstack)
	{

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			dstack.emplace(CELL(this->input_node_Qw[i],
													this->data->_surface[this->input_node_Qw[i]],
													this->input_Qw[i],
													this->input_Qs[i]));
		}
	}

	void run_subgraphflood()
	{

		if (this->data->_timetracker.size() == 0) {
			this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_debug = std::vector<f_t>(this->con->nxy(), 0.);
		}

		std::priority_queue<WaCell<i_t, f_t>,
												std::vector<WaCell<i_t, f_t>>,
												std::less<WaCell<i_t, f_t>>>
			dynastack;

		fillvec(this->data->_Qwin, 0.);
		fillvec(this->isInQ, false);

		this->init_dstack_Water<WaCell<i_t, f_t>, decltype(dynastack)>(dynastack);

		// CT_neighbourer_WaCell<i_t, f_t> ctx;
		CT_neighbours<i_t, f_t> ctx;

		// int nndt = 0;
		// for(auto v:this->data->_boundaries){
		// 	if(nodata(v)) ++nndt;
		// }
		// std::cout << "I have " << nndt << "no data " << std::endl;

		// fillvec(this->data->_vmot_hw,0.);
		this->time += this->dt;
		// std::vector<bool> isdone(this->con->nxy(), false);
		// int ndone = 0;
		// int nredone = 0;

		std::array<i_t, 8> receivers;
		std::array<f_t, 8> receiversWeights;

		// std::cout << "Starting the process" << std::endl;

		while (dynastack.empty() == false) {

			// Getting the next node
			auto next = this->_dstack_next<WaCell<i_t, f_t>>(dynastack);
			// auto next = dynastack.top();
			// dynastack.pop();
			// std::cout << next.node << "|" << this->data->_surface[next.node] <<
			// "||";

			this->isInQ[next.node] = false;

			bool ispast = this->data->_timetracker[next.node] != this->time;

			// if(isdone[next.node] == false){
			// 	ndone++;
			// 	isdone[next.node] = true;
			// 	if(ndone % 100 == 0)
			// 		std::cout << ndone << " vs " << nredone << " PQsizzla: " <<
			// dynastack.size() << " this->debugyolo " << this->debugyolo <<
			// std::endl; }else{ 	nredone++;
			// }

			// Updating the timer
			this->data->_timetracker[next.node] = this->time;

			// std::cout << next.node << std::endl;
			ctx.update(next.node, *this->con);
			;

			// std::cout << BC2str(ctx.boundary) << std::endl;

			if (ispast) {
				// this->data->_Qwin[next.node] = std::max(this->data->_Qwin[next.node],
				// next.Qw);
				this->data->_Qwin[next.node] += next.Qw;
			}

			if (std::isfinite(this->data->_Qwin[next.node]) == false)
				throw std::runtime_error("allah");
			// std::cout << ispast << std::endl;

			if (can_out(ctx.boundary)) {
				// std::cout << "out" << std::endl;
				continue;
			}

			if (nodata(ctx.boundary)) {
				// std::cout << "nodata reached" << std::endl;
				continue;
			}

			int nr = 0;
			f_t SS = 0;
			f_t SSdy = 1.;

			// std::cout << "A1" << std::endl;
			this->update_receivers(ctx, receivers, receiversWeights, nr, SS, SSdy);
			// std::cout << "A2" << std::endl;

			f_t& thw = this->data->_hw[next.node];
			f_t& tsurf = this->data->_surface[next.node];
			// Actual flux calculation

			f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
								std::sqrt(std::max(1e-6, SS));
			f_t tQwout = thw * u_w * SSdy;
			this->data->_Qwout[ctx.node] = tQwout;

			if (!ispast) {
				thw += this->hw_increment_LM;
				tsurf += this->hw_increment_LM;
			}

			f_t baseQw = (ispast) ? this->data->_Qwin[next.node] : next.Qw;

			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				int rec = ctx.neighbours[i];
				bool tizdone = this->data->_timetracker[rec] == this->time;
				if (tizdone) {
					dynastack.emplace(WaCell<i_t, f_t>(
						rec, this->data->_surface[rec], baseQw * receiversWeights[j]));
				} else {
					if (this->isInQ[rec] == false) {
						dynastack.emplace(WaCell<i_t, f_t>(rec, this->data->_surface[rec]));
						this->isInQ[rec] = true;
					}
					this->data->_Qwin[rec] += baseQw * receiversWeights[j];
				}
			}
		}
		// std::cout << "done" << std::endl;

		CT_neighbourer_WaCell<i_t, f_t> ctx2;

		for (int i = 0; i < this->con->nxy(); ++i) {

			if (nodata(this->data->_boundaries[i]) ||
					can_out(this->data->_boundaries[i]))
				continue;
			if (this->data->_timetracker[i]<this->time&& this->data->_hw[i]> 0) {
				this->data->_timetracker[i] = this->time;
				this->data->_Qwin[i] = 0.;
				ctx.update(i, *this->con);
				this->_calculate_Qwout_for_disconnected_nodes(ctx2);
			}
		}

		for (int i = 0; i < this->con->nxy(); ++i) {

			if (nodata(this->data->_boundaries[i]) ||
					can_out(this->data->_boundaries[i]))
				continue;

			f_t& thw = this->data->_hw[i];
			f_t& tsurf = this->data->_surface[i];

			f_t dhw = this->data->_Qwin[i] - this->data->_Qwout[i];
			dhw *= this->dt;
			dhw /= this->con->area(i);

			thw += dhw;
			tsurf += dhw;

			if (std::isfinite(tsurf) == false) {
				std::cout << this->data->_Qwout[i] << "|" << this->data->_Qwin[i]
									<< std::endl;
				;
				throw std::runtime_error("blug");
			}

			if (thw < 0) {
				tsurf -= thw;
				thw = 0;
			}
		}
	}

	template<class CTX>
	void update_receivers(CTX& ctx,
												std::array<i_t, 8>& receivers,
												std::array<f_t, 8>& receiversWeights,
												int& nr,
												f_t& SS,
												f_t& SSdy)
	{
		nr = 0;
		bool recout;
		bool recdone;

		bool allout;
		bool alldone;

		f_t sumSdw = 0.;
		SS = 0.;

		// pass 1: check the receivers, if no receivers fill up
		while (nr == 0) {
			recout = false;
			recdone = false;
			allout = true;
			alldone = true;

			for (int i = 0; i < ctx.nn; ++i) {
				if (this->data->_surface[ctx.node] >
						this->data->_surface[ctx.neighbours[i]]) {
					if (can_receive(ctx.neighboursCode[i])) {
						receivers[nr] = i;
						if (can_out(ctx.neighboursCode[i])) {
							recout = true;
						} else {
							allout = false;
						}
						if (this->time == this->data->_timetracker[ctx.neighbours[i]]) {
							recdone = true;
						} else {
							alldone = false;
						}
						++nr;
					}
				}
			}
			if (nr == 0) {
				f_t tinc = this->hw_increment_LM +
									 (this->data->randu->get() * this->hw_increment_LM) / 2;
				this->data->_surface[ctx.node] += tinc;
				this->data->_hw[ctx.node] += tinc;
			}
			// std::cout << this->data->_surface[ctx.node] << "|";
		}

		// pass2: recast to only the that can out
		if (recout) {
			int nnr = 0;
			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				if (can_out(ctx.neighboursCode[i])) {
					receivers[nnr] = receivers[j];
					++nnr;
				}
			}
			nr = nnr;
		}
		// Pass 3: if no out, some done but not all done
		else if (recdone && !alldone) {
			int nnr = 0;
			for (int j = 0; j < nr; ++j) {
				int i = receivers[j];
				if (this->time != this->data->_timetracker[ctx.neighbours[i]]) {
					receivers[nnr] = receivers[j];
					++nnr;
				}
			}
			nr = nnr;
		}
		// pass 4: all recs are done, I then choose a single one randomly
		else if (alldone) {
			int ti = static_cast<int>(this->data->randu->get() * nr);
			nr = 1;
			receivers[0] = receivers[ti];
		}

		// pass 5 precalculate SS and weights
		for (int j = 0; j < nr; ++j) {
			int i = receivers[j];
			f_t tS = this->get_Sw(ctx.node, ctx.neighbours[i], ctx.neighboursDx[i]);
			f_t tSdwinc = tS * ctx.neighboursDy[i];
			receiversWeights[j] = tSdwinc;
			sumSdw += tSdwinc;
			if (tS > SS) {
				SS = tS;
				SSdy = ctx.neighboursDy[i];
			}
		}

		// pass 6: calculate weights
		if (sumSdw > 0) {
			for (int j = 0; j < nr; ++j) {
				receiversWeights[j] /= sumSdw;
			}
		} else {
			for (int j = 0; j < nr; ++j) {
				receiversWeights[j] /= nr;
			}
		}

		// f_t sumsum = 0;
		// for(int i=0; i<nr; ++i){
		// 	if(receiversWeights[i] > 1)
		// 		throw std::runtime_error("false");
		// 	sumsum += receiversWeights[i];
		// }
		// if(std::abs(sumsum - 1) > 1e-4 && sumsum != 0)
		// 	throw std::runtime_error("false");

		// done
	}

	f_t get_Sw(int d, int r, f_t dx)
	{
		return (this->data->_surface[d] - this->data->_surface[r]) / dx;
	}

	template<class CELL>
	CELL _dstack_next(
		std::priority_queue<CELL, std::vector<CELL>, std::less<CELL>>& dynastack)
	{
		auto next = dynastack.top();
		// std::cout << next.Qw;

		dynastack.pop();

		if (dynastack.empty() == false) {
			while (dynastack.top().node == next.node) {
				next.ingest(dynastack.top());
				dynastack.pop();
				if (dynastack.empty())
					break;
			}
		}
		// std::cout << " | " << next.Qw << std::endl;;

		return next;
	}

	template<class CELL, class CTX, class Q_t>
	void _subGF_process_node(CELL& next, CTX& ctx, Q_t& dynastack)
	{

		f_t& thw = this->data->_hw[next.node];
		f_t& tsurf = this->data->_surface[next.node];

		// Regulating weights to avoid useless spreading
		this->regulate_weights(ctx);

		// Actual flux calculation

		f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
							std::sqrt(ctx.receiversSlopes[ctx.SSj]);
		f_t tQwout = thw * u_w * ctx.receiversDy[ctx.SSj];

		// std::cout << ctx.receiversSlopes[ctx.SSj] << "|";

		// f_t dhw = next.Qw - tQwout;
		// dhw *= this->dt;
		// dhw /= this->con->area(next.node);

		// thw += dhw;
		// tsurf += dhw;

		// if(dhw>1){
		// 	std::cout << dhw << std::endl;

		// // if(std::isfinite(thw) == false){
		// 	std::cout <<  next.Qw <<  "|" << tQwout << std::endl;
		// 	throw std::runtime_error("bite");
		// }

		// this->data->_Qwin[ctx.node] =
		// 	next.Qw; // std::max( this->data->_Qwin[ctx.node], next.Qw);
		this->data->_Qwout[ctx.node] = tQwout;

		// if (thw < 0) {
		// 	tsurf -= thw;
		// 	thw = 0;
		// }

		this->emplace_transfer(next, ctx, dynastack);
	}

	template<class CTX>
	void _calculate_Qwout_for_disconnected_nodes(CTX& ctx)
	{

		f_t& thw = this->data->_hw[ctx.node];
		f_t& tsurf = this->data->_surface[ctx.node];

		// Actual flux calculation

		if (ctx.receiversSlopes[ctx.SSj] <= 0)
			return;

		f_t u_w = std::pow(thw, (2. / 3.)) / this->mannings *
							std::sqrt(ctx.receiversSlopes[ctx.SSj]);
		f_t tQwout = thw * u_w * ctx.receiversDy[ctx.SSj];
		this->data->_Qwout[ctx.node] = tQwout;
	}

	template<class CELL, class CTX, class Q_t>
	void _subGF_reprocess_node(CELL& next, CTX& ctx, Q_t& dynastack)
	{

		f_t& thw = this->data->_hw[next.node];
		f_t& tsurf = this->data->_surface[next.node];

		// Regulating weights to avoid useless spreading
		this->regulate_weights(ctx);
		f_t tinc = this->hw_increment_LM * (this->data->randu->get() + 0.5);
		thw += tinc;
		tsurf += tinc;

		// this->data->_Qwin[ctx.node] =
		// 	next.Qw; // std::max( this->data->_Qwin[ctx.node], next.Qw);
		// this->data->_Qwin[ctx.node] = 0.;
		// this->data->_Qwout[ctx.node] = 0.;

		this->emplace_transfer(next, ctx, dynastack);
	}

	template<class CTX>
	void regulate_weights(CTX& ctx)
	{
		if (ctx.nr <= 1)
			return;

		if (this->min_part_Qw > 0) {
			f_t nSumSlopesDw = 0.;
			for (int i = 0; i < ctx.nr; ++i) {
				if (ctx.receiversWeights[i] < this->min_part_Qw)
					continue;
				nSumSlopesDw += ctx.receiversSlopes[i] * ctx.receiversDy[i];
			}

			if (nSumSlopesDw != 0. && nSumSlopesDw != ctx.sumslopesdw) {
				for (int i = 0; i < ctx.nr; ++i) {
					if (ctx.receiversWeights[i] < this->min_part_Qw) {
						ctx.receiversWeights[i] = 0.;
						ctx.receiversSlopes[i] = 0.;
					} else
						ctx.receiversWeights[i] =
							ctx.receiversSlopes[i] * ctx.receiversDy[i] / nSumSlopesDw;
				}
			}
			ctx.sumslopesdw = nSumSlopesDw;
		}
	}

	template<class CELL, class CTX, class Q_t>
	void emplace_transfer(CELL& next, CTX& ctx, Q_t& dynastack)
	{
		if (ctx.canout)
			return;

		if (ctx.nr == 0) {
			next.topo = this->data->_surface[next.node];
			dynastack.emplace(next);
			++this->debugyolo;
		} else {
			// std::cout << ctx.nr << "|";
			f_t sumW = 0;
			for (int i = 0; i < ctx.nr; ++i) {
				if (ctx.receiversWeights[i] > 0) {
					sumW += ctx.receiversWeights[i];
					WaCell<i_t, f_t> tnext;
					if (this->time != this->data->_timetracker[ctx.receivers[i]]) {
						tnext = WaCell<i_t, f_t>(ctx.receivers[i],
																		 this->data->_surface[ctx.receivers[i]]);

						this->data->_Qwin[ctx.receivers[i]] +=
							this->data->_Qwin[next.node] * ctx.receiversWeights[i];
					} else {
						tnext = WaCell<i_t, f_t>(ctx.receivers[i],
																		 this->data->_surface[ctx.receivers[i]],
																		 this->data->_Qwin[next.node] *
																			 ctx.receiversWeights[i]);
					}
					// if(this->data->_hw[tnext.node] == 0) std::cout << tnext.Qw <<
					// std::endl;
					// this->data->_debug[tnext.node] = 1;
					dynastack.emplace(tnext);
				}
			}

			if (abs(sumW - 1) > 1e-5)
				std::cout << sumW << std::endl;
		}
	}

	template<class CTX, class Q_t>
	void emplace_transfer(ExpCell<i_t, f_t>& next, CTX& ctx, Q_t& dynastack)
	{
		if (ctx.canout)
			return;

		if (ctx.nr == 0) {
			next.topo = this->data->_surface[next.node];
			dynastack.emplace(next);
		} else {
			// std::cout << ctx.nr << "|";
			f_t sumW = 0;
			for (int i = 0; i < ctx.nr; ++i) {
				if (ctx.receiversWeights[i] > 0) {
					sumW += ctx.receiversWeights[i];
					auto tnext = ExpCell<i_t, f_t>(ctx.receivers[i],
																				 this->data->_surface[ctx.receivers[i]],
																				 next.Qw * ctx.receiversWeights[i],
																				 next.Qs * ctx.receiversWeights[i]);
					dynastack.emplace(tnext);
				}
			}

			if (abs(sumW - 1) > 1e-5)
				std::cout << sumW << std::endl;
		}
	}

	void run_subgraphflood_expA(f_t propsed, f_t erate, f_t th_Qw, f_t stochED)
	{

		if (this->data->_timetracker.size() == 0) {
			this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_debug = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_Qsin = std::vector<f_t>(this->con->nxy(), 0.);
			this->data->_Qsout = std::vector<f_t>(this->con->nxy(), 0.);
		}

		std::priority_queue<ExpCell<i_t, f_t>,
												std::vector<ExpCell<i_t, f_t>>,
												std::less<ExpCell<i_t, f_t>>>
			dynastack;
		this->init_dstack_WaterSed<ExpCell<i_t, f_t>, decltype(dynastack)>(
			dynastack);
		CT_neighbourer_WaCell<i_t, f_t> ctx;

		// fillvec(this->data->_vmot_hw,0.);
		this->time += this->dt;

		while (dynastack.empty() == false) {

			// std::cout << dynastack.size() << "|";

			// Getting the next node
			auto next = this->_dstack_next<ExpCell<i_t, f_t>>(dynastack);

			bool ispast = this->data->_timetracker[next.node] != this->time;

			// Updating the timer
			this->data->_timetracker[next.node] = this->time;

			// std::cout << "A" << std::endl;
			ctx.update(next.node, *this->con);
			// std::cout << "B" << std::endl;

			if (can_out(ctx.boundary))
				continue;

			if (next.Qw < th_Qw) {
				// this->data->_Qsin[next.node] += next.Qs;
				continue;
			}

			//
			if (ctx.nr > 0) {
				// std::cout << "A" << std::endl;
				this->_subGF_process_node(next, ctx, dynastack);

				if (ispast) {
					f_t tprop = next.Qs / next.Qw;
					if (tprop < propsed) {
						f_t te = std::max(erate,
															(propsed - tprop) /
																(this->con->area(next.node) * this->dt)) *
										 (1 + (this->data->randu->get()) * 2 * stochED - stochED);
						this->data->_Qsout[next.node] += te * this->con->area(next.node);
						next.Qs += te * this->con->area(next.node);
					} else {
						f_t td = (tprop - propsed) /
										 (this->con->area(next.node) * this->dt) *
										 (1 + (this->data->randu->get()) * 2 * stochED - stochED);
						this->data->_Qsin[next.node] += td * this->con->area(next.node);
						next.Qs -= td * this->con->area(next.node);
					}
				}

				// std::cout << "Ad" << std::endl;
			} else {
				// std::cout << "B " << next.node << " " << this->data->_hw[next.node]
				// << " " << this->data->_surface[next.node]  << " " << ctx.nr << "|" <<
				// ctx.nn << std::endl;
				this->_subGF_reprocess_node(next, ctx, dynastack);
				// std::cout << "Bd:: " << this->data->_surface[next.node] << std::endl;
				// for(int j=0; j< ctx.nn; ++j){
				// 	std::cout << "|" << this->data->_surface[ctx.neighbours[j]];
				// }
			}
		}

		for (int i = 0; i < this->con->nxy(); ++i) {
			if (this->data->_timetracker[i]<this->time&& this->data->_hw[i]> 0) {
				this->data->_timetracker[i] = this->time;
				this->data->_Qwin[i] = 0.;
				ctx.update(i, *this->con);
				this->_calculate_Qwout_for_disconnected_nodes(ctx);
			}
		}

		for (int i = 0; i < this->con->nxy(); ++i) {

			f_t& thw = this->data->_hw[i];
			f_t& tsurf = this->data->_surface[i];

			f_t dhw = this->data->_Qwin[i] - this->data->_Qwout[i];
			dhw *= this->dt;
			dhw /= this->con->area(i);

			thw += dhw;
			tsurf += dhw;

			tsurf += (this->data->_Qsin[i] - this->data->_Qsout[i]) /
							 this->con->area(i) * this->dt;

			if (thw < 0) {
				tsurf -= thw;
				thw = 0;
			}

			this->data->_Qsin[i] = 0.;
			this->data->_Qsout[i] = 0.;
		}
	}

	f_t bedrockZatI(int i)
	{
		return this->data->_surface[i] - this->data->_hw[i];
	}

	// Access to steepest rec of node i but reprocessed on the spot (as opposed to
	// preprocessed from a compute() operation)
	int Sreceivers_bedrock_raw(i_t i,
														 std::array<i_t, 8>& arr,
														 std::array<f_t, 8>& arrdx,
														 i_t& trec,
														 f_t& tdx)
	{
		i_t nn = this->con->Neighbours(i, arr);
		nn = this->con->NeighboursDx(i, arrdx);
		trec = -1;
		tdx = this->con->_dx;
		f_t tSS = 0.;
		for (int j = 0; j < nn; ++j) {
			// arr[j] += i;

			f_t ttSS = (this->bedrockZatI(i) - this->bedrockZatI(arr[j])) / arrdx[j];

			if (ttSS > tSS) {
				tSS = ttSS;
				trec = arr[j];
				tdx = arrdx[j];
			}
		}

		return nn;
	}

	f_t _quick_slipos_from_point(i_t node,
															 f_t friction_slope,
															 f_t S_c,
															 f_t Kdep,
															 f_t Hpart,
															 f_t minH)
	{

		std::unordered_set<int> alreadyIn;

		std::queue<int> tQ;
		tQ.emplace(node);
		alreadyIn.insert(node);

		std::array<int, 8> tneigh;
		std::array<f_t, 8> tneighdx;

		f_t HLS = 0;

		// determine rupture_slope
		f_t rupture_slope = 0.;

		i_t nn = this->con->Neighbours(node, tneigh);
		this->con->NeighboursDx(node, tneighdx);

		f_t tsurf = this->bedrockZatI(node);

		for (int j = 0; j < nn; ++j) {
			i_t tn = tneigh[j];
			f_t tdx = tneighdx[j];
			std::cout << this->bedrockZatI(tn) << "vs" << tsurf << std::endl;
			if ((this->bedrockZatI(tn) - tsurf) / tdx > rupture_slope) {
				rupture_slope = (this->bedrockZatI(tn) - tsurf) / tdx;
			}
		}

		if (rupture_slope < friction_slope) {
			std::cout << "No ldsl" << std::endl;
			return 0.;
		}

		rupture_slope += friction_slope;
		rupture_slope /= 2;

		// Gather operation
		while (tQ.empty() == false) {
			i_t next = tQ.front();
			tQ.pop();
			// get donors
			i_t nn = this->con->Neighbours(next, tneigh);
			this->con->NeighboursDx(next, tneighdx);

			tsurf = this->data->_surface[next];

			for (int j = 0; j < nn; ++j) {
				i_t tn = tneigh[j];
				if (alreadyIn.find(tn) != alreadyIn.end())
					continue;
				f_t tdx = tneighdx[j];

				f_t tS = (this->bedrockZatI(tn) - tsurf) / tdx;
				if (tS > rupture_slope) {
					alreadyIn.insert(tn);
					tQ.emplace(tn);
					f_t tgz = tsurf + rupture_slope * tdx;
					HLS += (this->bedrockZatI(tn) - tgz);
					this->data->_surface[tn] = tgz + this->data->_hw[tn];
				}
			}
		}
		// return 0;

		i_t trec;
		f_t tdx;

		i_t Npart = std::floor(HLS / Hpart);
		std::cout << "n parts: " << Npart << std::endl;
		for (int i = 0; i < Npart; ++i) {
			f_t tH = Hpart;
			int nexnode = node;
			// this->Sreceivers_bedrock_raw(nexnode, tneigh, tneighdx, trec, tdx);
			// nexnode = trec;

			while (true) {
				this->Sreceivers_bedrock_raw(nexnode, tneigh, tneighdx, trec, tdx);

				if (can_out(this->data->_boundaries[nexnode]))
					break;

				if (trec == -1) {
					f_t delta_hk = tH * Kdep;
					this->data->_surface[nexnode] += delta_hk;
					tH -= delta_hk;

					if (tH < minH) {
						this->data->_surface[nexnode] += tH;
						break;
					}
					continue;
				}
				if (node == nexnode) {
					nexnode = trec;

					continue;
				}

				// std::cout << "YOLO:" << trec << std::endl;

				f_t tS = (this->bedrockZatI(nexnode) - this->bedrockZatI(trec)) / tdx;

				if (tS < S_c) {
					f_t delta_h = std::min((S_c - tS) * tdx, tH);
					this->data->_surface[nexnode] += delta_h;
					tH -= delta_h;
					// std::cout << delta_h << std::endl;
				}

				f_t delta_hk = tH * Kdep;
				this->data->_surface[nexnode] += delta_hk;
				tH -= delta_hk;
				// std::cout << delta_hk << "::" << std::endl;

				if (tH < minH) {
					this->data->_surface[nexnode] += tH;
					break;
				}

				nexnode = trec;
			}
		}

		return HLS;

		// Sreceivers_raw(i_t i, std::array<i_t, 8>& arr, std::array<f_t, 8>& arrdx,
		// i_t& trec, f_t& tdx)
	}

}; // end of class

} // end of namespace DAGGER

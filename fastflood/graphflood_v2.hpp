#include "dodcontexts.hpp"
#include "graphflood_enums.hpp"
#include "utils.hpp"

namespace DAGGER {

template<class T, class U>
class WaCell
{
public:
	// empty constructor
	WaCell(){};
	// Constructor by default
	WaCell(T node, U score, U Qw)
	{
		this->node = node;
		this->topo = score;
		this->Qw = Qw;
	};

	// Node index
	T node;
	// Score data
	U topo;
	U Qw;
	U Qs = 0.;

	// void ingest(WaCell<T,U>& other){this->Qw += other.Qw;}
	void ingest(WaCell<T, U> other) { this->Qw += other.Qw; }
};
;

// Custom operator sorting the nodes by scores
template<class T, class U>
inline bool
operator>(const WaCell<T, U>& lhs, const WaCell<T, U>& rhs)
{
	if (lhs.topo != rhs.topo)
		return lhs.topo > rhs.topo;
	else
		return lhs.node > rhs.node;
}

// Custom operator sorting the nodes by topos
template<class T, class U>
inline bool
operator<(const WaCell<T, U>& lhs, const WaCell<T, U>& rhs)
{
	if (lhs.topo != rhs.topo)
		return lhs.topo < rhs.topo;
	else
		return lhs.node < rhs.node;
}

template<class i_t, class f_t, class CONNECTOR_T, class GRAPH_T, class DBAG_T>
class Graphflood2
{
public:
	Graphflood2(){};

	Graphflood2(CONNECTOR_T& con, GRAPH_T& gra, DBAG_T& data)
	{
		this->con = &con;
		this->gra = &gra;
		this->data = &data;
	}

	CONNECTOR_T* con;
	GRAPH_T* gra;
	DBAG_T* data;

	WATER_INPUT water_input_mode = WATER_INPUT::PRECIPITATIONS_CONSTANT;
	f_t Prate = 1e-5; // mm.s^-1
	std::vector<i_t> input_node_Qw;
	std::vector<f_t> input_Qw;

	// Water transfer params
	f_t min_part_Qw = 0.01;

	// time management thingies
	// std::vector<f_t> data->_timetracker; // MOVED TO THE DATA BAG
	f_t time = 0.;
	f_t dt = 1e-3;

	// Graphflood params
	f_t mannings = 0.033;
	f_t hw_increment_LM = 1e-1;

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
		this->data->_vmot_hw = std::vector<f_t>(this->con->nxy(), 0);
	}

	void _computeEntryPoints_prec(f_t Qw_threshold)
	{

		// copying the bedrock elevation
		std::vector<f_t> bedrock(this->data->_surface);

		this->con->reinit();

		// filling and computing the graph from it
		this->con->PFcompute_all();

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

		for (int i = this->con->nxy() - 1; i >= 0; --i) {
			int node = this->data->_Sstack[i];
			if (nodata(this->data->_boundaries[node]))
				;
			int rec = this->con->Sreceivers(node);
			if (QW[rec] >= Qw_threshold && QW[node] < Qw_threshold) {
				input_node_Qw.emplace_back(node);
				input_Qw.emplace_back(QW[node]);
			} else if (QW[node] >= Qw_threshold) {
				input_node_Qw.emplace_back(node);
				input_Qw.emplace_back(Prate * this->con->area(node));
			}
		}
	}

	template<class DSTACK>
	void init_dstack_Water(DSTACK& dstack)
	{

		for (size_t i = 0; i < this->input_node_Qw.size(); ++i) {
			dstack.emplace(
				WaCell<i_t, f_t>(this->input_node_Qw[i],
												 this->data->_surface[this->input_node_Qw[i]],
												 this->input_Qw[i]));
		}
	}

	void run_subgraphflood()
	{

		if (this->data->_timetracker.size() == 0) {
			this->data->_timetracker = std::vector<f_t>(this->con->nxy(), 0.);
		}

		std::priority_queue<WaCell<i_t, f_t>,
												std::vector<WaCell<i_t, f_t>>,
												std::less<WaCell<i_t, f_t>>>
			dynastack;
		this->init_dstack_Water(dynastack);
		CT_neighbourer_WaCell<i_t, f_t> ctx;

		// fillvec(this->data->_vmot_hw,0.);
		this->time += this->dt;

		while (dynastack.empty() == false) {
			// Getting the next node
			auto next = this->_dstack_next(dynastack);

			bool ispast = this->data->_timetracker[next.node] != this->time;

			// Updating the timer
			this->data->_timetracker[next.node] = this->time;

			// std::cout << "A" << std::endl;
			ctx.update(next.node, *this->con);
			// std::cout << "B" << std::endl;

			if (can_out(ctx.boundary))
				continue;

			//
			if (ispast && ctx.nr > 0) {
				// std::cout << "A" << std::endl;
				this->_subGF_process_node(next, ctx, dynastack);
				// std::cout << "Ad" << std::endl;
			} else {
				continue;
				// std::cout << "B " << next.node << " " << this->data->_hw[next.node]
				// << " " << this->data->_surface[next.node]  << " " << ctx.nr << "|" <<
				// ctx.nn << std::endl;
				this->_subGF_reprocess_node(next, ctx, dynastack);
				// std::cout << "Bd" << std::endl;
			}
		}
	}

	WaCell<i_t, f_t> _dstack_next(
		std::priority_queue<WaCell<i_t, f_t>,
												std::vector<WaCell<i_t, f_t>>,
												std::less<WaCell<i_t, f_t>>>& dynastack)
	{
		auto next = dynastack.top();

		dynastack.pop();

		if (dynastack.empty() == false) {
			while (dynastack.top().node == next.node) {
				// std::cout << next.Qw;
				next.ingest(dynastack.top());
				// std::cout << " | " << next.Qw << std::endl;;
				dynastack.pop();
				if (dynastack.empty())
					break;
			}
		}

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
		f_t dhw = next.Qw - tQwout;
		dhw *= this->dt;
		dhw /= this->con->area(next.node);

		thw += dhw;
		tsurf += dhw;

		// if(dhw>1){
		// 	std::cout << dhw << std::endl;

		// // if(std::isfinite(thw) == false){
		// 	std::cout <<  next.Qw <<  "|" << tQwout << std::endl;
		// 	throw std::runtime_error("bite");
		// }

		this->data->_Qwin[ctx.node] = next.Qw;

		if (thw < 0) {
			tsurf -= thw;
			thw = 0;
		}

		this->emplace_transfer(next, ctx, dynastack);
	}

	template<class CELL, class CTX, class Q_t>
	void _subGF_reprocess_node(CELL& next, CTX& ctx, Q_t& dynastack)
	{

		f_t& thw = this->data->_hw[next.node];
		f_t& tsurf = this->data->_surface[next.node];

		// Regulating weights to avoid useless spreading
		this->regulate_weights(ctx);
		thw += this->hw_increment_LM;
		tsurf += this->hw_increment_LM;

		this->emplace_transfer(next, ctx, dynastack);
	}

	template<class CTX>
	void regulate_weights(CTX& ctx)
	{
		// if(ctx.nr<=1)
		return;

		if (this->min_part_Qw > 0) {
			f_t nSumSlopesDw = 0.;
			for (int i = 0; i < ctx.nr; ++i) {
				if (ctx.receiversWeights[i] < this->min_part_Qw)
					continue;
				nSumSlopesDw = ctx.receiversSlopes[i] * ctx.receiversDy[i];
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

		if (ctx.nr == 0 && can_out(ctx.boundary) == false) {
			next.topo = this->data->_surface[next.node];
			dynastack.emplace(next);
		} else {
			f_t sumW = 0;
			for (int i = 0; i < ctx.nr; ++i) {
				if (ctx.receiversWeights[i] > 0) {
					sumW += ctx.receiversWeights[i];
					dynastack.emplace(
						WaCell<i_t, f_t>(ctx.receivers[i],
														 this->data->_surface[ctx.receivers[i]],
														 next.Qw * ctx.receiversWeights[i]));
				}
			}
			if (abs(sumW - 1) > 1e-5)
				std::cout << sumW << std::endl;
		}
	}
};

} // end of namespace DAGGER

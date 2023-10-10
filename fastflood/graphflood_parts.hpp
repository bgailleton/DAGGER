#pragma once
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

template<class T, class U>
class ExpCell // iarmus ... lol
{
public:
	// empty constructor
	ExpCell(){};
	// Constructor by default
	ExpCell(T node, U score, U Qw, U Qs)
	{
		this->node = node;
		this->topo = score;
		this->Qw = Qw;
		this->Qs = Qs;
	};

	// Node index
	T node;
	// Score data
	U topo;
	U Qw;
	U Qs = 0.;

	// void ingest(ExpCell<T,U>& other){this->Qw += other.Qw;}
	void ingest(ExpCell<T, U> other) { this->Qw += other.Qw; }
};
;

// Custom operator sorting the nodes by scores
template<class T, class U>
inline bool
operator>(const ExpCell<T, U>& lhs, const ExpCell<T, U>& rhs)
{
	if (lhs.topo != rhs.topo)
		return lhs.topo > rhs.topo;
	else
		return lhs.node > rhs.node;
}

// Custom operator sorting the nodes by topos
template<class T, class U>
inline bool
operator<(const ExpCell<T, U>& lhs, const ExpCell<T, U>& rhs)
{
	if (lhs.topo != rhs.topo)
		return lhs.topo < rhs.topo;
	else
		return lhs.node < rhs.node;
}

} // End of the namespace DAGGER

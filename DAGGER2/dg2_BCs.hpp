/*
Holds all the info about boundary conditions

Using the Boundary codes of scabbard/DAGGER:


*/

#pragma once

#include <string>
#include <unordered_map>

namespace dagger2 {

/**
 * NodeType: Defines boundary conditions for flux and data flow through nodes
 *
 * This enum class specifies how water, sediment, and other materials behave
 * at each node in the landscape grid. It controls flow behavior for:
 * - Hydraulic processes (water flow, drainage)
 * - Erosion and sediment transport
 * - Heat transfer, chemical diffusion, etc.
 */
enum class NodeType : uint8_t
{
	/**
	 * NO_DATA: Invalid or uninitialized node
	 * - No flux processing occurs
	 * - Node is effectively removed from calculations
	 * - Used for masking areas outside the domain
	 */
	NO_DATA = 0,

	/**
	 * NORMAL: Standard internal node
	 * - Flux can flow in and out freely
	 * - Follows standard flow routing algorithms
	 * - Most common node type for interior domain
	 */
	NORMAL = 1,

	/**
	 * HAS_TO_OUT: Forced outlet node
	 * - Any flux reaching this node must exit the domain
	 * - Cannot redistribute to neighboring nodes
	 * - Used for river outlets, boundary sinks
	 * - Flux is permanently removed from system
	 */
	HAS_TO_OUT = 2,

	/**
	 * CAN_OUT: Optional outlet node
	 * - Flux can either exit domain OR continue to neighbors
	 * - Decision based on flow routing algorithm
	 * - Allows natural flow patterns near boundaries
	 * - More flexible than HAS_TO_OUT
	 */
	CAN_OUT = 3,

	/**
	 * IN: Inlet/source node
	 * - Flux enters the domain at this node
	 * - Cannot exit through this node (flux is trapped)
	 * - Used for springs, rainfall input, sediment sources
	 * - May have constant or time-varying input rates
	 */
	IN = 4,

	/**
	 * PERIODIC: Periodic boundary condition
	 * - Flux exiting this node re-enters at corresponding periodic node
	 * - Maintains conservation across domain boundaries
	 * - Used for simulating infinite or wrapped domains
	 * - Requires proper pairing with opposite boundary
	 */
	PERIODIC = 5,

	/**
	 * REFLECT: Reflective boundary condition
	 * - Flux hitting this node is reflected back
	 * - No flux crosses the boundary
	 * - Equivalent to no-flux or impermeable boundary
	 * - Used for watershed divides, impermeable barriers
	 */
	REFLECT = 6
};

/**
 * Utility functions for NodeType operations
 */
class NodeTypeUtils
{
public:
	/**
	 * Convert NodeType to human-readable string
	 */
	static std::string to_string(NodeType type)
	{
		static const std::unordered_map<NodeType, std::string> type_names = {
			{ NodeType::NO_DATA, "NO_DATA" },
			{ NodeType::NORMAL, "NORMAL" },
			{ NodeType::HAS_TO_OUT, "HAS_TO_OUT" },
			{ NodeType::CAN_OUT, "CAN_OUT" },
			{ NodeType::IN, "IN" },
			{ NodeType::PERIODIC, "PERIODIC" },
			{ NodeType::REFLECT, "REFLECT" }
		};

		auto it = type_names.find(type);
		return (it != type_names.end()) ? it->second : "UNKNOWN";
	}

	/**
	 * Parse string to NodeType
	 */
	static NodeType from_string(const std::string& str)
	{
		static const std::unordered_map<std::string, NodeType> string_to_type = {
			{ "NO_DATA", NodeType::NO_DATA },
			{ "NORMAL", NodeType::NORMAL },
			{ "HAS_TO_OUT", NodeType::HAS_TO_OUT },
			{ "CAN_OUT", NodeType::CAN_OUT },
			{ "IN", NodeType::IN },
			{ "PERIODIC", NodeType::PERIODIC },
			{ "REFLECT", NodeType::REFLECT }
		};

		auto it = string_to_type.find(str);
		if (it != string_to_type.end()) {
			return it->second;
		}
		throw std::invalid_argument("Unknown NodeType string: " + str);
	}

	/**
	 * Check if node allows flux to exit the domain
	 */
	static bool allows_outflow(NodeType type)
	{
		return type == NodeType::HAS_TO_OUT || type == NodeType::CAN_OUT ||
					 type == NodeType::PERIODIC;
	}

	/**
	 * Check if node allows flux to enter the domain
	 */
	static bool allows_inflow(NodeType type)
	{
		return type == NodeType::NORMAL || type == NodeType::CAN_OUT ||
					 type == NodeType::IN || type == NodeType::PERIODIC;
	}

	/**
	 * Check if node processes flux at all
	 */
	static bool is_active(NodeType type) { return type != NodeType::NO_DATA; }

	/**
	 * Check if node is a boundary condition
	 */
	static bool is_boundary(NodeType type)
	{
		return type == NodeType::HAS_TO_OUT || type == NodeType::CAN_OUT ||
					 type == NodeType::IN || type == NodeType::PERIODIC ||
					 type == NodeType::REFLECT;
	}

	/**
	 * Get default flux behavior description
	 */
	static std::string get_description(NodeType type)
	{
		switch (type) {
			case NodeType::NO_DATA:
				return "No flux processing - node is inactive";
			case NodeType::NORMAL:
				return "Standard node - flux flows freely in all directions";
			case NodeType::HAS_TO_OUT:
				return "Forced outlet - all flux must exit domain";
			case NodeType::CAN_OUT:
				return "Optional outlet - flux can exit or continue to neighbors";
			case NodeType::IN:
				return "Inlet/source - flux enters but cannot exit";
			case NodeType::PERIODIC:
				return "Periodic boundary - flux wraps to opposite boundary";
			case NodeType::REFLECT:
				return "Reflective boundary - flux bounces back";
			default:
				return "Unknown node type";
		}
	}
};

} // namespace dagger2

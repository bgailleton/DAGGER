#pragma once

#include "dg2_BCs.hpp"
#include "dg2_array.hpp"
#include "dg2_connector.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace dagger2 {

template<typename T>
class HydraulicErosionEngine
{
public:
	struct ErosionParams
	{
		int numDroplets = 50000;
		float inertia = 0.1f;
		float capacity = 4.0f;
		float deposition = 0.3f;
		float erosion = 0.3f;
		float evaporation = 0.01f;
		float radius = 1.5f;
		float minSlope = 0.001f;
		int maxLifetime = 64;
		float initialWater = 1.0f;
		float initialSpeed = 1.0f;
	};

private:
	struct Droplet
	{
		float x, y;
		float dx, dy;
		float speed;
		float water;
		float sediment;
	};

	std::shared_ptr<Connector<T>> connector_;
	std::mt19937 rng_;
	ErosionParams params_;

public:
	HydraulicErosionEngine(std::shared_ptr<Connector<T>> connector,
												 unsigned int seed = std::random_device{}())
		: connector_(connector)
		, rng_(seed)
	{
	}

	void setParams(const ErosionParams& params) { params_ = params; }
	const ErosionParams& getParams() const { return params_; }

	void preset_light_erosion()
	{
		params_.numDroplets = 25000;
		params_.erosion = 0.1f;
		params_.deposition = 0.2f;
		params_.capacity = 2.0f;
	}

	void preset_moderate_erosion()
	{
		params_.numDroplets = 50000;
		params_.erosion = 0.3f;
		params_.deposition = 0.3f;
		params_.capacity = 4.0f;
	}

	void preset_heavy_erosion()
	{
		params_.numDroplets = 100000;
		params_.erosion = 0.5f;
		params_.deposition = 0.2f;
		params_.capacity = 6.0f;
	}

	void erode(std::shared_ptr<Grid2D<T>> heightmap)
	{
		if (!heightmap || heightmap->rows() != connector_->rows() ||
				heightmap->cols() != connector_->cols()) {
			throw std::invalid_argument("Invalid heightmap dimensions");
		}

		for (int i = 0; i < params_.numDroplets; ++i) {
			simulateDroplet(*heightmap);
		}
	}

	std::shared_ptr<Grid2D<T>> erode_copy(std::shared_ptr<Grid2D<T>> heightmap)
	{
		auto result = std::make_shared<Grid2D<T>>(*heightmap);
		erode(result);
		return result;
	}

private:
	void simulateDroplet(Grid2D<T>& heightmap)
	{
		Droplet drop = createDroplet();

		for (int lifetime = 0; lifetime < params_.maxLifetime; ++lifetime) {
			int xi = static_cast<int>(drop.x);
			int yi = static_cast<int>(drop.y);

			if (!connector_->is_valid_coord(yi, xi))
				break;

			// Check if location can output
			NodeType bc = connector_->get_boundary_type(yi, xi);
			bool can_output = NodeTypeUtils::allows_outflow(bc);
			if (can_output)
				break;

			// Get height and gradient
			float height = sampleHeight(heightmap, drop.x, drop.y);
			auto [gradX, gradY] = calculateGradient(heightmap, xi, yi);

			// Update velocity (blend old direction with gradient)
			drop.dx = drop.dx * params_.inertia - gradX * (1.0f - params_.inertia);
			drop.dy = drop.dy * params_.inertia - gradY * (1.0f - params_.inertia);

			// Check for sufficient flow
			float speed = std::sqrt(drop.dx * drop.dx + drop.dy * drop.dy);
			if (speed < params_.minSlope)
				break;

			// Normalize velocity
			drop.dx /= speed;
			drop.dy /= speed;

			// Calculate new position with boundary handling
			float newX = drop.x + drop.dx;
			float newY = drop.y + drop.dy;

			// Handle boundary conditions
			if (!connector_->is_valid_coord(static_cast<int>(newY),
																			static_cast<int>(newX))) {
				// Determine direction for boundary handling
				Direction dir;
				if (std::abs(drop.dx) > std::abs(drop.dy)) {
					dir = (drop.dx > 0) ? Direction::EAST : Direction::WEST;
				} else {
					dir = (drop.dy > 0) ? Direction::SOUTH : Direction::NORTH;
				}

				auto neighbor = connector_->get_effective_neighbor(yi, xi, dir);
				if (!neighbor.is_valid)
					break;

				auto [nr, nc] = connector_->to_2d(neighbor.index);
				newX = static_cast<float>(nc) + (drop.x - xi);
				newY = static_cast<float>(nr) + (drop.y - yi);
			}

			float newHeight = sampleHeight(heightmap, newX, newY);

			// Height difference (positive = downhill)
			float deltaHeight = height - newHeight;

			// Update speed from gravity
			drop.speed =
				std::sqrt(drop.speed * drop.speed + std::max(0.0f, deltaHeight));

			// Calculate sediment capacity
			float capacity = std::max(0.0f, deltaHeight) * drop.speed * drop.water *
											 params_.capacity;

			// Erosion/deposition at current position
			if (drop.sediment > capacity) {
				// Deposit excess sediment
				float deposit = (drop.sediment - capacity) * params_.deposition;
				drop.sediment -= deposit;
				depositSediment(heightmap, drop.x, drop.y, deposit);
			} else if (deltaHeight > 0) {
				// Erode when moving downhill
				float maxErode = deltaHeight * 0.5f;
				float erode =
					std::min((capacity - drop.sediment) * params_.erosion, maxErode);
				drop.sediment += erode;
				erodeTerrain(heightmap, drop.x, drop.y, erode);
			}

			// Move droplet
			drop.x = newX;
			drop.y = newY;

			// Evaporate water
			drop.water *= (1.0f - params_.evaporation);
			if (drop.water < 0.01f)
				break;
		}
	}

	Droplet createDroplet()
	{
		std::uniform_real_distribution<float> xDist(
			params_.radius, connector_->cols() - params_.radius - 1);
		std::uniform_real_distribution<float> yDist(
			params_.radius, connector_->rows() - params_.radius - 1);
		std::uniform_real_distribution<float> dirDist(-1.0f, 1.0f);

		Droplet drop;
		drop.x = xDist(rng_);
		drop.y = yDist(rng_);
		drop.dx = dirDist(rng_);
		drop.dy = dirDist(rng_);

		// Normalize initial direction
		float len = std::sqrt(drop.dx * drop.dx + drop.dy * drop.dy);
		if (len > 0) {
			drop.dx /= len;
			drop.dy /= len;
		} else {
			drop.dx = 1.0f;
			drop.dy = 0.0f;
		}

		drop.speed = params_.initialSpeed;
		drop.water = params_.initialWater;
		drop.sediment = 0.0f;

		return drop;
	}

	std::pair<float, float> calculateGradient(const Grid2D<T>& heightmap,
																						int x,
																						int y)
	{
		float gx = 0.0f, gy = 0.0f;

		// X gradient - only calculate if both neighbors exist
		auto east = connector_->get_effective_neighbor(y, x, Direction::EAST);
		auto west = connector_->get_effective_neighbor(y, x, Direction::WEST);

		if (east.is_valid && west.is_valid) {
			auto [er, ec] = connector_->to_2d(east.index);
			auto [wr, wc] = connector_->to_2d(west.index);
			float h_east = static_cast<float>(heightmap(er, ec));
			float h_west = static_cast<float>(heightmap(wr, wc));
			gx = (h_west - h_east) * 0.5f;
		}

		// Y gradient - only calculate if both neighbors exist
		auto north = connector_->get_effective_neighbor(y, x, Direction::NORTH);
		auto south = connector_->get_effective_neighbor(y, x, Direction::SOUTH);

		if (north.is_valid && south.is_valid) {
			auto [nr, nc] = connector_->to_2d(north.index);
			auto [sr, sc] = connector_->to_2d(south.index);
			float h_north = static_cast<float>(heightmap(nr, nc));
			float h_south = static_cast<float>(heightmap(sr, sc));
			gy = (h_north - h_south) * 0.5f;
		}

		return { gx, gy };
	}

	float sampleHeight(const Grid2D<T>& heightmap, float x, float y)
	{
		int x0 = static_cast<int>(std::floor(x));
		int y0 = static_cast<int>(std::floor(y));
		int x1 = x0 + 1;
		int y1 = y0 + 1;

		// Clamp to valid bounds
		x0 = std::max(0, std::min(x0, static_cast<int>(connector_->cols()) - 1));
		y0 = std::max(0, std::min(y0, static_cast<int>(connector_->rows()) - 1));
		x1 = std::max(0, std::min(x1, static_cast<int>(connector_->cols()) - 1));
		y1 = std::max(0, std::min(y1, static_cast<int>(connector_->rows()) - 1));

		float fx = x - std::floor(x);
		float fy = y - std::floor(y);

		float h00 = static_cast<float>(heightmap(y0, x0));
		float h10 = static_cast<float>(heightmap(y0, x1));
		float h01 = static_cast<float>(heightmap(y1, x0));
		float h11 = static_cast<float>(heightmap(y1, x1));

		// Bilinear interpolation
		float h0 = h00 * (1.0f - fx) + h10 * fx;
		float h1 = h01 * (1.0f - fx) + h11 * fx;
		return h0 * (1.0f - fy) + h1 * fy;
	}

	void depositSediment(Grid2D<T>& heightmap, float x, float y, float amount)
	{
		if (std::isnan(amount) || std::isinf(amount) || amount <= 0)
			return;

		int cx = static_cast<int>(std::round(x));
		int cy = static_cast<int>(std::round(y));
		int radius = static_cast<int>(params_.radius);

		for (int dy = -radius; dy <= radius; ++dy) {
			for (int dx = -radius; dx <= radius; ++dx) {
				int nx = cx + dx;
				int ny = cy + dy;

				if (!connector_->is_valid_coord(ny, nx) ||
						!connector_->is_active_node(ny, nx))
					continue;

				float dist = std::sqrt(dx * dx + dy * dy);
				if (dist > params_.radius)
					continue;

				float weight = std::max(0.0f, 1.0f - dist / params_.radius);
				float delta = amount * weight / (params_.radius * params_.radius);

				heightmap(ny, nx) += static_cast<T>(delta);
			}
		}
	}

	void erodeTerrain(Grid2D<T>& heightmap, float x, float y, float amount)
	{
		if (std::isnan(amount) || std::isinf(amount) || amount <= 0)
			return;

		int cx = static_cast<int>(std::round(x));
		int cy = static_cast<int>(std::round(y));
		int radius = static_cast<int>(params_.radius);

		for (int dy = -radius; dy <= radius; ++dy) {
			for (int dx = -radius; dx <= radius; ++dx) {
				int nx = cx + dx;
				int ny = cy + dy;

				if (!connector_->is_valid_coord(ny, nx) ||
						!connector_->is_active_node(ny, nx))
					continue;

				float dist = std::sqrt(dx * dx + dy * dy);
				if (dist > params_.radius)
					continue;

				float weight = std::max(0.0f, 1.0f - dist / params_.radius);
				float delta = amount * weight / (params_.radius * params_.radius);

				heightmap(ny, nx) -= static_cast<T>(delta);
			}
		}
	}
};

} // namespace dagger2

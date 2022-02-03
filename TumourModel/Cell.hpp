#pragma once
#include <assert.h>

#define _USE_MATH_DEFINES

#include "Helper.hpp"
#include <math.h>
#include <vector>
#include <map>

#include <iostream>
#include <random>

static std::random_device rd;
static std::mt19937 gen(rd());

using namespace::std;

enum class CellType { Cancer = 1, Dead = 3, Healthy = 4, Empty = 5, PSC = 51 };

class Cell;

struct CellState
{
	CellType type;
	int age;
	double spring_length;
	Cell* sibling;
	double X, Y;
};

class Cell
{
public:
	vector<Cell*> Neighbours;
	CellState currentState, newState;

    Cell(double X, double Y, double spring_length, Cell* sibling, CellType type, int age)
	{
		this->currentState.X = X;
		this->currentState.Y = Y;
		this->currentState.spring_length = spring_length;
		this->currentState.sibling = sibling;
		this->currentState.type = type;
		this->currentState.age = age;
	}

	Cell(Cell* cell)
	{
		currentState = cell->currentState;
	}

	void Renew()
	{
		newState = currentState;
	}

	void UpdateState()
	{
		currentState = newState;
		currentState.age++;
	}

	void Infect()
	{
		currentState.type = CellType::Cancer;
	}

	void updateSibling(map<Cell*, Cell*> &map)
	{
		if (currentState.sibling != NULL)
			newState.sibling = map[currentState.sibling];
	}

	void clearNeighbours()
	{
		Neighbours.clear();
	}

	bool OnBoundary()
	{
		if (currentState.type == CellType::Healthy)
			return false;
		else
			for (Cell* neighbour : Neighbours)
				if (neighbour->currentState.type == CellType::Healthy)
					return true;
		return false;
	}

	bool TooCrowded()
	{
		// if there is more than one neighbour within distance threshold
		int closeNeighbours = 0;
		for (Cell* neighbour : Neighbours)
			if (DistanceSquaredTo(neighbour) < Sqr(Params::rmin))
			{
				closeNeighbours++;
				if (closeNeighbours > 1)
					return true;
			}
		return false;
	}

	double DistanceSquaredTo(Cell* cell)
	{
		return DistanceSquared(currentState.X, currentState.Y, cell->currentState.X, cell->currentState.Y);
	}

	double DistanceSquaredFromCentre()
	{
		return DistanceSquared(currentState.X, currentState.Y, 0, 0);
	}

	double DistanceFromBoundary(vector<Cell*> &boundaryCells)
	{
		double minDistance = MAX_DBL;
		for (Cell* boundaryCell : boundaryCells)
		{
			double distance = DistanceSquaredTo(boundaryCell);
			if (distance < minDistance)
				minDistance = distance;
		}
		return sqrt(minDistance);
	}

	bool Necrotic(double distanceToBoundary, Params* parameters)
	{
		// a long way from any boundary cell
		return distanceToBoundary >= parameters->dmax;
	}

	bool TooYoung(Params* parameters)
	{
		return currentState.age < parameters->gage;
	}

	void LengthenSpring(Params* parameters)
	{
		if (currentState.spring_length < Params::s)
			if (currentState.sibling != NULL)
				newState.spring_length = currentState.spring_length + Params::s / parameters->page;
			else
			{
				newState.spring_length = Params::s;
				newState.sibling = NULL;
			}
		else
			newState.sibling = NULL;
	}

	bool DrugInducedDeath(Params* parameters, double* drugConcentration, int gridRadius)
	{
		assert(-gridRadius <= currentState.X && currentState.X <= gridRadius);
		assert(-gridRadius <= currentState.Y && currentState.Y <= gridRadius);

		int nearestX = round(currentState.X) + gridRadius;
		int nearestY = round(currentState.Y) + gridRadius;

		double concentration = drugConcentration[nearestY * (gridRadius * 2 + 1) + nearestX];
		if (concentration > 0)
			return parameters->WithProbability(concentration / (concentration + parameters->EC50));
		else
			return false;
	}

	void Die()
	{
		newState.type = CellType::Dead;
		newState.age = 0;
	}

	void Disintegrate()
	{
		int dage = 3; //number of hours it takes for cell to completely disintegrate
		if (currentState.age > dage)
			newState.type = CellType::Empty;
		else
			newState.spring_length /= 8;
	}

	void PossiblyPSCInfectNeighbour(Params* parameters)
	{
		if (parameters->WithProbability(parameters->p_psc))
		{
			vector<Cell*> healthyNeighbours;
			for (Cell* neighbour : Neighbours)
				if (neighbour->currentState.type == CellType::Healthy)
					healthyNeighbours.push_back(neighbour);
			Cell* randomNeighbour = healthyNeighbours[(int)(parameters->RandomDouble() * healthyNeighbours.size())];
			randomNeighbour->newState.type = CellType::PSC;
		}
	}

	void Move()
	{
		double Fx = 0, Fy = 0;
		for (Cell* neighbour : Neighbours)
		{
			double dist = sqrt(DistanceSquaredTo(neighbour));
			double dx = (neighbour->currentState.X - currentState.X) / dist;
			double dy = (neighbour->currentState.Y - currentState.Y) / dist;

			double s = Params::s;
			if (neighbour->currentState.spring_length < Params::s && neighbour == currentState.sibling)
				s = neighbour->currentState.spring_length;

			double diff = dist - s;
			if (diff > 0.05)
				diff = 0;

			Fx += Params::mu * diff * dx;
			Fy += Params::mu * diff * dy;
		}

		newState.X = currentState.X + Fx * Params::Delta_t;
		newState.Y = currentState.Y + Fy * Params::Delta_t;
	}

	Cell* PossiblyPoliferate(vector<Cell*> &boundaryCells, Params* parameters)
	{
		if (!TooYoung(parameters) && !TooCrowded())
		{
			double distanceToBoundry = DistanceFromBoundary(boundaryCells);
			if (!Necrotic(distanceToBoundry, parameters) && parameters->WithProbability(parameters->p_0 * (1 - distanceToBoundry / parameters->dmax)))
			    return Proliferate(parameters);
		}

		return NULL;
	}

	Cell* Proliferate(Params* parameters)
	{
		double theta = 2 * M_PI * parameters->RandomDouble();
		double scale = Params::s / parameters->page / 2;
		newState.sibling = new Cell(currentState.X + scale * cos(theta), currentState.Y + scale * sin(theta), Params::s / parameters->page, this, CellType::Cancer, 1);
		newState.X = currentState.X + scale * cos(M_PI + theta);
		newState.Y = currentState.Y + scale * sin(M_PI + theta);
		newState.spring_length = Params::s / parameters->page;
		
		/*double age_prob = parameters->RandomDouble();
		if(age_prob<0.1)
		{newState.age = 1;}
		else if(age_prob<0.2)
		{newState.age = 2;}
		else if(age_prob<0.3)
		{newState.age = 3;}
		else if(age_prob<0.4)
		{newState.age = 4;}
		else if(age_prob<0.5)
		{newState.age = 5;}
		else if(age_prob<0.6)
		{newState.age = 6;}
		else if(age_prob<0.7)
		{newState.age = 7;}
		else if(age_prob<0.8)
		{newState.age = 8;}
		else if(age_prob<0.9)
		{newState.age = 9;}
		else 
		{newState.age = 10;}*/
	

		std::poisson_distribution<int> pd(50);
		
		newState.age = pd(gen);
			
		return newState.sibling;
	}
};

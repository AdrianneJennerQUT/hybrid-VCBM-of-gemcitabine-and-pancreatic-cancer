#pragma once

#include "Params.hpp"
#include "Cell.hpp"
#include "Diffuse.hpp"
#include "delaunator.hpp"
#include <vector>
#include <set>
#include <assert.h>
#include <chrono>


#include <iostream>
#include <random>


using namespace::std;

double totalTime = 0;

class Pancreas
{
public:
	Pancreas(Pancreas *existing, Params* parameters)
	{
		this->parameters = parameters;

		memcpy(drugConcentration, existing->drugConcentration, sizeof(drugConcentration));
		this->timeSinceInjection = existing->timeSinceInjection;
		this->fibreX = existing->fibreX;
		this->fibreY = existing->fibreY;

		// Clone the input cells in case those cells are used in more than one particle
		map<Cell*, Cell*> map;

		for (Cell* cell : existing->cells)
		{
			Cell* clone = new Cell(cell);
			cells.push_back(clone);
			map.insert(std::map<Cell*, Cell*>::value_type(cell, clone));
		}

		for (Cell* cell : cells)
			cell->updateSibling(map);

		src = dst = NULL;
	}

	Pancreas(Params* parameters)
	{
		this->parameters = parameters;
		memset(drugConcentration, 0, sizeof(drugConcentration));
		timeSinceInjection = -1;
		fibreX = fibreY = 0;
		src = dst = NULL;
	}

	~Pancreas()
	{
		//delete parameters; // is this right???

		for (Cell* cell : cells)
			delete cell;
	}

	Pancreas* CreateNewParticle(Params* parameters)
	{
		return new Pancreas(this, parameters);
	}

#define gridRadius	70
#define gridWidth	(gridRadius*2+1)

	void InjectPoint(int x, int y, double amount)
	{
		timeSinceInjection = 0;
		drugConcentration[(gridRadius + y) * gridWidth + (gridRadius + x)] += amount;
	}

	void InjectFibre(int x, int y, double amount)
	{
		// fibre starts at location (x,y) and extends horizontally to the right
		fibreX = x;
		fibreY = y;
		timeSinceInjection = 0;
		for (int i = 0; i < Params::N-1; i++)
			fibreConcentration[i] = amount;
		fibreConcentration[Params::N - 1] = 0;
	}

private:

	vector<Cell*> cells;
	Params* parameters;
	vector<Cell*> new_cells; // accumulated during the latest iteration
	vector<Cell*> boundaryCells;
	Cell* src, *dst; // the two cells that are furthest apart - using for computing dimensions of tumour

public:
	int fibreX, fibreY;
	long timeSinceInjection;
	double drugConcentration[gridWidth * gridWidth + Params::N];
	double* fibreConcentration = drugConcentration + gridWidth * gridWidth;

#define nextHalfedge(e) ((e % 3 == 2) ? e - 2 : e + 1)
#define triangleOfEdge(e) (e/3)

	double* LoadCellsCoordinates()
	{
		double* coords = new double[cells.size() * 2];
		int i = 0;
		for (Cell* cell : cells)
		{
			cell->clearNeighbours();
			coords[i++] = cell->currentState.X;
			coords[i++] = cell->currentState.Y;
		}
		return coords;
	}

	void circumCentre42(double a0, double a1, double b0, double b1, double c0, double c1, void (*vertex)(double, double))
	{
		double ad = a0 * a0 + a1 * a1;
		double bd = b0 * b0 + b1 * b1;
		double cd = c0 * c0 + c1 * c1;
		double D = 2 * (a0 * (b1 - c1) + b0 * (c1 - a1) + c0 * (a1 - b1));
		double XX = 1 / D * (ad * (b1 - c1) + bd * (c1 - a1) + cd * (a1 - b1));
		double YY = 1 / D * (ad * (c0 - b0) + bd * (a0 - c0) + cd * (b0 - a0));
		vertex(XX, YY);
	}

	void circumCentre(double a0, double a1, double b0, double b1, double c0, double c1, void (*vertex)(double, double))
	{
		double accumulatedArea = 0.0f;
		double centerX = 0.0f;
		double centerY = 0.0f;

		{
			double temp = a0 * c1 - c0 * a1;
			accumulatedArea += temp;
			centerX += (a0 + c0) * temp;
			centerY += (a1 + c1) * temp;
		}
		{
			double temp = b0 * a1 - a0 * b1;
			accumulatedArea += temp;
			centerX += (b0 + a0) * temp;
			centerY += (b1 + a1) * temp;
		}
		{
			double temp = c0 * b1 - b0 * c1;
			accumulatedArea += temp;
			centerX += (c0 + b0) * temp;
			centerY += (c1 + b1) * temp;
		}

		if (abs(accumulatedArea) < 1E-7)
		{
			//vertex(0, 0);
		}
		else
		{
			accumulatedArea *= 3;
			vertex(centerX / accumulatedArea, centerY / accumulatedArea);
		}	
	}

	void EnumerateVoronoiCells(void (*startCell)(CellType), void (*vertex)(double,double), void (*endCell)())
	{
		double* coords = LoadCellsCoordinates();
		int points = cells.size();
		Delaunator TRI(coords, points);
		delete[] coords;

		int seen[10000];
		memset(seen, 0, sizeof(seen));

		for (int e = 0; e < TRI.trianglesLen; e++) 
		{
			int id = TRI.triangles[nextHalfedge(e)];

			if (!seen[id] && cells[id]->currentState.type != CellType::Healthy)
			{
				seen[id] = 1;

				int incoming = e;
				startCell(cells[id]->currentState.type);
				do
				{
					int t = triangleOfEdge(incoming);
					Cell* a = cells[TRI.triangles[3 * t]];
					Cell* b = cells[TRI.triangles[3 * t + 1]];
					Cell* c = cells[TRI.triangles[3 * t + 2]];
					
					double a0 = a->currentState.X;
					double a1 = a->currentState.Y;
					double b0 = b->currentState.X;
					double b1 = b->currentState.Y;
					double c0 = c->currentState.X;
					double c1 = c->currentState.Y;
					circumCentre(a0, a1, b0, b1, c0, c1, vertex);

					int outgoing = nextHalfedge(incoming);
					incoming = TRI.halfedges[outgoing];
				} 			
				while (incoming != -1 && incoming != e);
				endCell();
			}
		}
	}

	void EnumerateNeighbours(void (*line)(Cell*, Cell*))
	{
		double* coords = LoadCellsCoordinates();
		Delaunator TRI(coords, cells.size());
		delete[] coords;

		for (size_t e = 0; e < TRI.trianglesLen; e++)
			if (TRI.halfedges[e] == -1 || e > TRI.halfedges[e])
			{
				Cell* P = cells[TRI.triangles[e]];
				Cell* Q = cells[TRI.triangles[nextHalfedge(e)]];
				line(P, Q);
			}
	}

	void DetermineNeighbours()
	{
		double* coords = LoadCellsCoordinates();
		Delaunator TRI(coords, cells.size());
		delete[] coords;

		for (size_t e = 0; e < TRI.trianglesLen; e++)
			if (TRI.halfedges[e] == -1 || e > TRI.halfedges[e])
			{
				Cell* P = cells[TRI.triangles[e]];
				Cell* Q = cells[TRI.triangles[nextHalfedge(e)]];
				if (P->currentState.type != CellType::Empty && Q->currentState.type != CellType::Empty)
				{
					P->Neighbours.push_back(Q);
					Q->Neighbours.push_back(P);
				}
			}
	}

	bool HealthyCellsBeyondRadius(double radius)
	{
		for (Cell* cell : cells)
			if (cell->currentState.type == CellType::Healthy && cell->DistanceSquaredFromCentre() >= Sqr(radius))
				return true;
		return false;
	}

	void AddNewCell(Cell* new_cell)
	{
		if (new_cell != NULL)
			new_cells.push_back(new_cell);
	}

	void AddMoreTissue(double moving_rim, double max_tumour_radius)
	{
		const double Y_spacing = 3;
		double X_spacing = Y_spacing * sqrt(3) / 2;

		double new_radius = max_tumour_radius + moving_rim + 10;
		int i = 0;
		for (double X = -new_radius; X <= new_radius; X += X_spacing)
		{
			// should offset odd columns by 1 or Y_spacing/2 ???
			for (double Y = -new_radius + (i % 2) * Y_spacing / 2; Y <= new_radius; Y += Y_spacing)
			{
				double distance_squared = DistanceSquared(X, Y, 0, 0);
				if (distance_squared > Sqr(max_tumour_radius) && distance_squared <= Sqr(max_tumour_radius + moving_rim)) // why +10 vs +20 above?
					AddNewCell(new Cell(X, Y, Params::s, NULL, CellType::Healthy, parameters->gage));
			}
			i++;
		}
	}
	void MoreTissueAddedIfNecessary()
	{
		const int moving_rim = 10;
		double tumour_radius = TumourRadius();
		if (!HealthyCellsBeyondRadius(tumour_radius + moving_rim))
			AddMoreTissue(moving_rim, tumour_radius);
	}

	void DetermineBoundaryCells()
	{
		if (cells[0]->Neighbours.size() == 0)
			DetermineNeighbours();

		boundaryCells.clear();
		for (Cell* cell : cells)
			if (cell->currentState.type != CellType::Healthy)
			{
				for (Cell* neighbour: cell->Neighbours)
					if (neighbour->currentState.type == CellType::Healthy)
					{
						boundaryCells.push_back(cell);
						break;
					}
			}
	}

	// Still used for determining when to add more cells
	double TumourRadius()
	{
		double maxDistance = 0;
		for (Cell* cell : cells)
			if (cell->currentState.type != CellType::Healthy)
			{
				double distance = cell->DistanceSquaredFromCentre();
				if (distance > maxDistance)
					maxDistance = distance;
			}
		return sqrt(maxDistance);
	}

	double DistanceToLine(Cell* cell)
	{
		return (dst->currentState.X - src->currentState.X)  * (src->currentState.Y - cell->currentState.Y) - 
			   (src->currentState.X - cell->currentState.X) * (dst->currentState.Y - src->currentState.Y);
	}


public:
	double TumourVolume()
	{
		if (boundaryCells.size() == 0)
			DetermineBoundaryCells();

		if (boundaryCells.size() <= 1)
			return 0;

		double longest = -1;

		// Find the two cells that are furthest apart from each other
		for (Cell* cell1 : boundaryCells)
		{
			for (Cell* cell2 : boundaryCells)
			{
				double distSquared = cell1->DistanceSquaredTo(cell2);
				if (distSquared > longest)
				{
					longest = distSquared;
					src = cell1;
					dst = cell2;
				}
			}
		}

		double length = sqrt(longest);

		// Find width from cells furthest from the centre line
		double minDist = 0, maxDist = 0;
		for (Cell* cell : boundaryCells)
		{
			double distance = DistanceToLine(cell);
			if (distance < minDist) minDist = distance;
			if (distance > maxDist) maxDist = distance;
		}

		double scale = sqrt(Sqr(dst->currentState.X - src->currentState.X) + Sqr(dst->currentState.Y - src->currentState.Y));
		minDist /= scale;
		maxDist /= scale;

		double width = maxDist - minDist;

		assert(width <= length);

		width *= 0.1728;
		length *= 0.1728;

		double volume = width * width * length / 2;
		//assert(!isnan(volume));
		return volume;
	}

	double getPscRatio()
	{
		int cancer = 0;
		int psc = 0;
		for (Cell* cell : cells)
		{
			if (cell->currentState.type == CellType::PSC)
				psc++;
			if (cell->currentState.type != CellType::Healthy)
				cancer++;
		}
		return (double)psc / cancer;
	}

	void CreateInitialTumour()
	{
		AddMoreTissue(5,0);
		// append the new cells to the current cells ...
		cells.insert(cells.end(), new_cells.begin(), new_cells.end());
		new_cells.clear();

		Cell* closestDistanceToCentre = NULL;
		double closestDistance = MAX_DBL;
		for (Cell* cell : cells)
		{
			double distance = cell->DistanceSquaredFromCentre();
			if (distance < closestDistance)
			{
				closestDistanceToCentre = cell;
				closestDistance = distance;
			}
		}
		closestDistanceToCentre->Infect();
	}

	void SimulateOneHour()
	{
		DetermineBoundaryCells();

		if (timeSinceInjection >= 0)
		{
			DiffuseSimple(timeSinceInjection, timeSinceInjection + 60, drugConcentration, gridWidth, fibreX, fibreY);
			timeSinceInjection += 60;
		}

		for (Cell* cell : cells)
		{
			cell->Renew(); // set new state to current state

			if (cell->OnBoundary())
				cell->PossiblyPSCInfectNeighbour(parameters);

			if (cell->currentState.type == CellType::Healthy)
				cell->Move();
			else if (cell->currentState.type == CellType::Dead)
				cell->Disintegrate();
			else if (cell->currentState.type == CellType::Empty)
			{
			}
			else
			{
				if (cell->DrugInducedDeath(parameters, drugConcentration, gridRadius))
					cell->Die();
				else
				{
					if (cell->currentState.age < parameters->gage)
						cell->LengthenSpring(parameters);

					// cancer cells either proliferate or move ...
					Cell* newCell = cell->PossiblyPoliferate(boundaryCells, parameters);
					if (newCell != NULL)
						AddNewCell(newCell);
					else
						cell->Move();
				}
			}
		}

		for (Cell* cell : cells)
			cell->UpdateState(); // set current state to new state

		MoreTissueAddedIfNecessary();

		// Only add new cells after all processing for this hour is complete.
		cells.insert(cells.end(), new_cells.begin(), new_cells.end());
		new_cells.clear();

		DetermineNeighbours();
	}

	double SimulateOneDay(int day, void (*render)(int, int, Pancreas*, int))
	{
		DetermineNeighbours();
		
		for (int hour = 0; hour < Params::tinterval; hour++)
		{
			if (render) render(day, hour, this, gridRadius);
			SimulateOneHour();
		}
		if (render) render(day, Params::tinterval, this, gridRadius);
		return TumourVolume();
	}
    
    double SimulateOneDay(int day)
	{
        return SimulateOneDay(day, NULL);
	}
	
	int ReturnTotalNumberTumourCells()
	{
		int num_of_cancer_cells = 0;
		for (Cell* cell : cells)
			if (cell->currentState.type == CellType::Cancer)
				num_of_cancer_cells++;
		return num_of_cancer_cells;
	}
	
	int ReturnTotalNumberDeadCells()
	{
		int num_of_dead_cells = 0;	
		for (Cell* cell : cells)
			if (cell->currentState.type == CellType::Dead)
				num_of_dead_cells++;
		return num_of_dead_cells;
	}
	
	int ReturnTotalNumberPSCCells()
	{
		int num_of_PSC_cells = 0;	
		for (Cell* cell : cells)
			if (cell->currentState.type == CellType::PSC)
				num_of_PSC_cells++;		
		return num_of_PSC_cells;
	}
	
	int ReturnTotalNumberHealthyCells()
	{
		int num_of_Healthy_cells = 0;	
		for (Cell* cell : cells)
			if (cell->currentState.type == CellType::Healthy)
				num_of_Healthy_cells++;
		
		return num_of_Healthy_cells;
	}
	
	int TestingPoissonDist()
	{
		std::poisson_distribution<int> pd(parameters->gage);
		
		int newage = pd(gen);
		return newage;
	}
	
	double ReturnDrugConcentrationDomain()
	{
		double drug_conc_total = 0;
		for(int ii = 0; ii <gridWidth * gridWidth; ii++)
			drug_conc_total += drugConcentration[ii];
		return drug_conc_total;
		
	}
	
	double ReturnDrugConcentrationinFibre()
	{
		double drug_conc_fibre_total = 0;	
		for(int ii = 0; ii <Params::N; ii++)
			drug_conc_fibre_total += fibreConcentration[ii];
		return drug_conc_fibre_total;
		
	}
	
	double ReturnDrugConcentrationAout()
	{
		double drug_conc_Aout = fibreConcentration[Params::N-1];
		return drug_conc_Aout;
		
	}
	
	double ReturnCellPositions(int index)
	{
		double* coords = LoadCellsCoordinates();
		return coords[index];
	}
	int ReturnNumberCells()
	{
		return cells.size();
	}
	int ReturnCellType(int index)
	{
		//Cancer = 1, Dead = 3, Healthy = 4, Empty = 5, PSC = 51 
		if (cells[index]->currentState.type == CellType::Healthy)
			return 4;
		else if (cells[index]->currentState.type == CellType::Cancer)
			return 1;
		else if (cells[index]->currentState.type == CellType::Dead)
			return 3;
		else if (cells[index]->currentState.type == CellType::Empty)
			return 5;
		else if (cells[index]->currentState.type == CellType::PSC)
			return 51;
		else 
				return 6;
	}
	
};

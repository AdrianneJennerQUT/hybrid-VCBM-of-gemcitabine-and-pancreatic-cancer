#include "Model.h"

Pancreas* SeedAndGrowToStartVolume(double p0, double psc, int dmax, int gage, int page, double EC50, double startVolume, void(*render)(int, int, Pancreas*, int))
{
	Params* parameters = new Params(p0, psc, dmax, gage, page, EC50);
	Pancreas* pancreas = new Pancreas(parameters);
	// start with just one infected cancer cell nearest to (0, 0)
	pancreas->CreateInitialTumour();

	// pre-observation phase - run until tumour reaches start volume
	double volume = 0;
	int days = 0;
	while (volume < startVolume && days < 200)
		volume = pancreas->SimulateOneDay(days++, render);

	// who disposes parameters???
	return pancreas;
}

Pancreas* SeedAndGrowToStartVolumeM(double p0, double psc, int dmax, int gage, int page, double EC50, double startVolume)
{
    return SeedAndGrowToStartVolume(p0, psc, dmax, gage, page, EC50, startVolume, NULL);
}

void SimulateWholeExperiment(double p0, double psc, int dmax, int gage, int page, double EC50, double startVolume, int timeSteps, double volumes[], void(*render)(int, int, Pancreas*, int))
{
	Pancreas* pancreas = SeedAndGrowToStartVolume(p0, psc, dmax, gage, page, EC50, startVolume, render);

	volumes[0] = pancreas->TumourVolume();

	//pancreas->InjectPoint(0, 0, 10);

	pancreas->InjectFibre(10, 10, 2000 * Params::C0 / (Params::N + 1));

	for (int day = 1; day < timeSteps; day++)
	{
		volumes[day] = pancreas->SimulateOneDay(day, render);

		if (volumes[day] > 10000) // terminate early if growth is out of control
		{
			for (int i = day + 1; i < timeSteps; i++)
				volumes[i] = 1000000;
			break;
		}
	}

	delete pancreas;
}

void SimulateWholeExperimentM(double p0, double psc, int dmax, int gage, int page, double EC50, double startVolume, int timeSteps, double volumes[])
{
    SimulateWholeExperiment(p0, psc, dmax, gage, page, EC50, startVolume, timeSteps, volumes, NULL);
}

void PerformMultipleRuns(double p0, double psc, int dmax, int gage, int page, double EC50, double startVolume, int timeSteps, int iterations, double volumes[], void(*render)(int, int, Pancreas*, int))
{
	for (int j = 0; j < timeSteps; j++)
		volumes[j] = 0;

	for (int i = 0; i < iterations; i++)
	{
		double* v = new double[timeSteps];
		SimulateWholeExperiment(p0, psc, dmax, gage, page, EC50, startVolume, timeSteps, v, render);
		for (int j = 0; j < timeSteps; j++)
			volumes[j] += v[j];

		delete[] v;
	}

	for (int j = 0; j < timeSteps; j++)
		volumes[j] /= iterations;
}

void PerformMultipleRunsM(double p0, double psc, int dmax, int gage, int page, double EC50, double startVolume, int timeSteps, int iterations, double volumes[])
{
    PerformMultipleRuns(p0, psc, dmax, gage, page, EC50, startVolume, timeSteps, iterations, volumes, NULL);
}

Pancreas* CreateNewParticle(double p0, double psc, int dmax, int gage, int page, double EC50, Pancreas* pancreas)
{
	return pancreas->CreateNewParticle(new Params(p0, psc, dmax, gage, page, EC50));
}
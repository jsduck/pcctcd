#pragma once

//-- VTK Includes
#include <vtkCommand.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>

class vtkUpdateCallback : public vtkCommand
{
public:
	static vtkUpdateCallback *New() { return new vtkUpdateCallback; }
	virtual void Execute(vtkObject* caller, unsigned long, void*);
	
	void* p;
	vtkPolyData* poly_data;
	std::vector<std::pair<vtkSphereSource*, double>> weighted_spheres;
};
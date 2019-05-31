#pragma once

//-- VTK Includes
#include <vtkCommand.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyle.h>

class vtkAntTweakBar : public vtkCommand {
public:	
	static vtkAntTweakBar* New();
	void Initialize(vtkRenderWindowInteractor* interactor);
	void Execute(vtkObject* caller, unsigned long event, void *) VTK_OVERRIDE;
	void Terminate();
	void Delete() VTK_OVERRIDE;

	void EnableAlwaysRenderOnEvent();
	void DisableAlwaysRenderOnEvent();
	bool AlwaysRenderOnEvent() const;
private:
	vtkAntTweakBar();

	bool initialized_;
	bool alwaysRender_;
	int mouseWheelPosition_;
	vtkRenderWindowInteractor* interactor_;
	vtkInteractorStyle* interactorStyle_;
};
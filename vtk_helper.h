#pragma once

#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL2)
#include "vtkAutoInit.h" 

// Initialise graphical modules before any object construction
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkVertex.h>
//#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkPolyLine.h>
#include <vtkPolyVertex.h>
#include <vtkSmartPointer.h>

#include <vtkSphereSource.h>
#include <vtkTextProperty.h>
#include <vtkProperty2D.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
//#include <vtkSliderWidget.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkCommand.h>
//#include <vtkWidgetEvent.h>
#include <vtkCallbackCommand.h>
//#include <vtkWidgetEventTranslator.h>
#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkSliderWidget.h>
//#include <vtkSliderRepresentation2D.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkOpenGLPolyDataMapper.h>
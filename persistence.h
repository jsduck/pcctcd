#pragma once

//-- AntTweakBar Includes
#include <AntTweakBar.h>

//-- Local Includes
#include "persistence_helper.h"
#include "crystal_str.h"
#include "vtk_helper.h"
#include "vtkAntTweakBar.hpp"

class Persistence
{
public:
	//Persistence();
	/*
	 * Construction taking in a Crystal object
	 */
	Persistence(Crystal& c);

	/*
	 * Set the crystal object used for the computation of the alpha shape
	 */
	void setCrystal(Crystal c) { c_ = c; }

	/*
	 * Calculate persistent cohomology
	 * param1: string of the .off points file -- DEPRECATED --
	 * param2: string of the output file -- DEPRECATED --
	 * param3: int value of the coefficient
	 * param4: double value for the minimum persistence kept in the pair vector
	 */
	void calculate(std::string off_file_points, std::string output_file_diag, int coeff_field_characteristic, Filtration_value min_persistence);
	/*
	 * Calculate clusters of persistent pairs -- DEPRECATED --
	 * param1: double value of the tolerance for cluster pairs
	 */
	void calculate_clusters(double tolerance);

	/*
	 * Output the complex polygon data to a .vtk file
	 * param1: string of the output filepath
	 */
	void output_complex(const std::string& path) const;
	/*
	 * Output the persistent pairs to a file
	 * param1: string of the output filepath
	 */
	void output_pairs(const std::string& path);
	//void output_intervals(const std::string& path) const;
	//void output_heatmap(const std::string& path) const;
	/*
	 * Update the master output file with the computed most persistent pair
	 * param1: string of the output filepath
	 */
	void output(const std::string& path);

	//void TW_CALL set_alpha(const void* alpha, void* data);
	//void TW_CALL get_alpha(void* alpha, void* data);

private:
	/*
	 * Initialise window objects
	 */
	void init_win();
	/*
	 * Initialise renderer and interactor objects
	 */
	void init_ren();
	/*
	 * Initialise and create camera object
	 */
	void init_camera();
	/*
	 * Initialise the complex actor object and its polygonal data
	 */
	void init_complex_actor();
	/*
	 * Initialise the spheres used to display points
	 */
	void init_point_sphere_actors();
	/*
	 * Initialise the spheres used for displaying alpha balls
	 */
	void init_radii_sphere_actors();
	/*
	 * Master method that calls other actor initialisation methods
	 */
	void init_actors();
	/*
	 * Initialise the VTK slider widget -- DEPRECATED --
	 */
	void init_slider();
	/*
	 * Master method that calls other widget initialisation methods
	 */
	void init_widgets();
	/*
	 * Initialise VTK slider callback -- DEPRECATED --
	 */
	void init_slider_callback();
	/*
	 * Initialise the actor data update callback
	 */
	void init_update_callback();
	/*
	 * Initialise AntTweakBar callback
	 */
	void init_atb_callback();
	/*
	 * Master method that calls other callback initialisation methods
	 */
	void init_callbacks();
	/*
	 * Initialise AntTweakbar and create options menu
	 */
	void init_abt();

public:
	/* 
	 * Visualise the computed alpha shape 
	 */
	void display_complex();

	/*
	 * Retrieve the persitent pairs
	 * return: vecotr containing pairs with Birth, Death double values
	 */
	std::vector<std::pair<double, double>> get_pairs() const;
	/*
	 * Retrieve present alpha values
	 * return: vector of double containg the alpha spectrum
	 */
	std::vector<double> get_alpha_values() const;

	/*
	 * Update the point map from an alpha shape
	 */
	void update_point_map(Alpha_shape_3* alpha_shape);
	/*
	 * Retrieve the current weighted point_map
	 * return: map of points and their indices
	 */
	std::map<Weighted_point_3, int> get_point_map() const;

	/*
	 * Retrieve all points present in the alpha shape
	 * param1: pointer to the alpha shape object
	 * return: VTK pointer of points obtained from the alpha shape object
	 */
	vtkSmartPointer<vtkPoints>		get_alpha_points(Alpha_shape_3* alpha_shape);
	/*
	 * Retrieve all edges present in the alpha shape
	 * param1: pointer to the alpha shape object
	 * return: VTK pointer of cells obtained from the alpha shape object
	 */
	vtkSmartPointer<vtkCellArray>	get_alpha_edges(Alpha_shape_3* alpha_shape, Alpha_shape_3::Classification_type type = Alpha_shape_3::SINGULAR);
	/*
	 * Retrieve all faces present in the alpha shape
	 * param1: pointer to the alpha shape object
	 * return: VTK pointer of cells obtained from the alpha shape object
	 */
	vtkSmartPointer<vtkCellArray>	get_alpha_faces(Alpha_shape_3* alpha_shape);
	/*
	 * Update the current polygonal data
	 * param1: pointer to the alpha shape object
	 * param2: alpha value for the extraction of the PolyData
	 * return: VTK pointer of the new PolyData object
	 */
	vtkSmartPointer<vtkPolyData>	update_poly_data(Alpha_shape_3* alpha_shape, double alpha = 0);

	//void vtkUpdateCallback(vtkObject* caller, long unsigned int id, void* clientData, void*);

	/*
	 * Clear all array data
	 */
	void clean();

	/*
	 * Retrieve Weighted points vector
	 * return: vector of Weighted_points_3
	 */
	std::vector<Weighted_point_3> get_weighted_points() const;

	/*
	 * Set alpha for the alpha shape
	 * param1: double value of the desired alpha
	 */
	void set_alpha(double alpha);
	/*
	 * Retrieve the alpha of the current alpha shape
	 * return: double value of alpha
	 */
	double get_alpha() const;

	/*
	 * Retrieve a pointer to the complex actor
	 * return: pointer to the complex actor
	 */
	vtkActor* GetComplexActor() const { return complex_actor_; }
	/*
	 * Retrieve the renderer
	 * return: pointer to the current renderer
	 */
	vtkRenderer* GetRenderer() const { return ren_; }
	/*
	 * Retrieve the window object
	 * return: pointer to the current window
	 */
	vtkRenderWindow* GetRenderWindow() const { return win_; }

	/*
	 * Retrieve the actors of the spheres used to display alpha balls
	 * return: vector of pointers to each sphere
	 */
	std::vector<vtkActor*> GetRadiiSpheresActors() const { return radii_sphere_actors_; }

private:
	bool update_flag = false;
	double max_alpha_ = 1;
	double min_alpha_ = 0;
	Crystal c_;
	Alpha_shape_3* as_;
	double alpha_;
	bool render_flag = false;
	std::vector<std::pair<double, double>> pairs_;
	std::map<Weighted_point_3, int> point_map_;
	std::vector<Weighted_point_3> points_;
	std::vector<double> pairs_vals_;
	std::vector<double> alpha_vals_;
	std::vector<int> clusters_;

	vtkSmartPointer<vtkRenderer> ren_;
	vtkSmartPointer<vtkRenderWindow> win_;
	vtkSmartPointer<vtkRenderWindowInteractor> inter_;

	vtkSmartPointer<vtkActor> complex_actor_;
	vtkSmartPointer<vtkPolyData> complex_source_;

	std::vector<vtkActor*> point_sphere_actors_;
	std::vector<vtkActor*> radii_sphere_actors_;

	std::vector<vtkSphereSource*> point_sphere_sources_;
	std::vector<std::pair<vtkSphereSource*, double>> radii_sphere_sources_; 

	//vtkSmartPointer<vtkSliderWidget> slider_widget_;

	//friend class vtkSliderCallback;
	friend class vtkUpdateCallback;

	//vtkSmartPointer<vtkUpdateCallback> update_callback_;
	vtkSmartPointer<vtkAntTweakBar> atb_callback_;

	//friend void TW_CALL set_alpha(const void * alpha, void * data);
	//friend void TW_CALL get_alpha(void * alpha, void * data);
};


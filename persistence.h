#pragma once

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
	Persistence(Crystal& c);

	void setCrystal(Crystal c) { c_ = c; }

	void calculate(std::string off_file_points, std::string output_file_diag, int coeff_field_characteristic, Filtration_value min_persistence);
	void calculate_clusters(double tolerance);

	void output_complex(const std::string& path);
	void output_pairs(const std::string& path);
	void output(const std::string& path);

	//void TW_CALL set_alpha(const void* alpha, void* data);
	//void TW_CALL get_alpha(void* alpha, void* data);

	void init_win();
	void init_ren();
	void init_complex_actor();
	void init_point_sphere_actors();
	void init_radii_sphere_actors();
	void init_actors();
	void init_slider();
	void init_widgets();
	void init_slider_callback();
	void init_atb_callback();
	void init_callbacks();
	void display_complex();

	std::vector<std::pair<double, double>> get_pairs() const;
	std::vector<double> get_alpha_values() const;

	void update_point_map(Alpha_shape_3 *alpha_shape);
	std::map<Weighted_point_3, int> get_point_map() const;

	vtkSmartPointer<vtkPoints>		get_alpha_points(Alpha_shape_3 *alpha_shape);
	vtkSmartPointer<vtkCellArray>	get_alpha_edges(Alpha_shape_3 *alpha_shape);
	vtkSmartPointer<vtkCellArray>	get_alpha_faces(Alpha_shape_3 *alpha_shape);
	vtkSmartPointer<vtkPolyData>	update_poly_data(Alpha_shape_3 *alpha_shape, double alpha = 0);

	//void vtkUpdateCallback(vtkObject* caller, long unsigned int id, void* clientData, void*);

	std::vector<Weighted_point_3> get_weighted_points() const;

	void set_alpha(const double alpha) {
		alpha_ = alpha;
		//std::cout << "Attempting to reconstruct complex with alpha as : " << alpha << std::endl;
		update_flag = true;
	}

	double get_alpha() const { return alpha_; }

private:
	bool update_flag = false;
	Crystal c_;
	Alpha_shape_3* as_;
	double alpha_;
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

	std::vector<vtkSphereSource*> point_sphere_sources_;
	std::vector<std::pair<vtkSphereSource*, double>> radii_sphere_sources_; 

	vtkSmartPointer<vtkSliderWidget> slider_widget_;

	friend class vtkSliderCallback;
	friend class vtkUpdateCallback;

	vtkSmartPointer<vtkAntTweakBar> atb_callback_;

	//friend void TW_CALL set_alpha(const void * alpha, void * data);
	//friend void TW_CALL get_alpha(void * alpha, void * data);
};


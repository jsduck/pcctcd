//-- Base Includes
#include "persistence.h"
#include "atom_data.h"
#include "file_helper.h"

#include <boost/assign/list_of.hpp>

//#include "vtkAntTweakBar.hpp"

#include "vtkAutoInit.h" 
//#include "atom_data.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

class vtkSliderCallback : public vtkCommand
{
public:
	static vtkSliderCallback *New() { return new vtkSliderCallback; }
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkSliderWidget *slider_widget =
			reinterpret_cast<vtkSliderWidget*>(caller);

		double val = pow(static_cast<vtkSliderRepresentation *>(slider_widget->GetRepresentation())->GetValue(), 2);
		p->alpha_ = val;
		//sliderWidget->GetRepresentation()->

		p->update_flag = true;

	}
	//vtkPolyData* poly_data;
	Persistence* p;
	////std::vector<vtkSmartPointer<vtkSphereSource>> spheres;
	//std::vector<std::pair<vtkSphereSource*, double>> weighted_spheres;
};

class vtkUpdateCallback : public vtkCommand
{
public:
	static vtkUpdateCallback *New() { return new vtkUpdateCallback; }
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		//vtkRenderWindowInteractor* interactor =
		//	reinterpret_cast<vtkRenderWindowInteractor*>(caller);
		//std::cout << p->update_flag << std::endl;
		if (p->update_flag) {
			p->update_poly_data(p->as_, p->alpha_);
			auto pd = p->complex_source_;
			poly_data->SetPoints(pd->GetPoints());
			poly_data->SetPolys(pd->GetPolys());
			poly_data->SetLines(pd->GetLines());

			//auto points = p_->get_points();
			for (auto &sphere : weighted_spheres) {
				//sphere->SetRadius(alpha); // + points[i].weight()
				sphere.first->SetRadius(sqrt(p->alpha_ + sphere.second));
			}

			p->update_flag = false;
		}

		//sliderWidget->GetRepresentation()->

		/*p->update_poly_data(p->as_, val);
		auto pd = p->complex_source_;
		poly_data->SetPoints(pd->GetPoints());
		poly_data->SetPolys(pd->GetPolys());
		poly_data->SetLines(pd->GetLines());

		//auto points = p_->get_points();
		for (auto &sphere : weighted_spheres) {
			//sphere->SetRadius(alpha); // + points[i].weight()
			sphere.first->SetRadius(sqrt(val + sphere.second));
		}
		*/
	}

	vtkPolyData* poly_data;
	Persistence* p;
	//std::vector<vtkSmartPointer<vtkSphereSource>> spheres;
	std::vector<std::pair<vtkSphereSource*, double>> weighted_spheres;
};

Persistence::Persistence(Crystal& c) {
	c_ = c;
	as_ = new Alpha_shape_3;
	//alpha_ = 0;
}

struct cmp_intervals_by_dim_then_length {
	explicit cmp_intervals_by_dim_then_length(ST * sc)
		: sc_(sc) { }
	template<typename Persistent_interval>
	bool operator()(const Persistent_interval & p1, const Persistent_interval & p2) {
		if (sc_->dimension(get < 0 >(p1)) == sc_->dimension(get < 0 >(p2)))
			return (sc_->filtration(get < 1 >(p1)) - sc_->filtration(get < 0 >(p1))
> sc_->filtration(get < 1 >(p2)) - sc_->filtration(get < 0 >(p2)));
		else
			return (sc_->dimension(get < 0 >(p1)) > sc_->dimension(get < 0 >(p2)));
	}
	ST* sc_;
};

void Persistence::calculate(std::string off_file_points, std::string output_file_diag, int coeff_field_characteristic, Filtration_value min_persistence) {

	//Gudhi::Points_3D_off_reader<Point_3> off_reader(off_file_points);
	// Check the read operation was correct
	//if (!off_reader.is_valid()) {
	//	std::cerr << "Unable to read OFF file " << off_file_points << std::endl;
	//	exit(-1);
	//}
	// Retrieve the points
	//std::vector<Point_3> lp = off_reader.get_point_cloud();

	//std::ifstream weights_ifstr(weight_file);
	std::vector<Weighted_point_3> wp;
	//auto points = c_.get_points();
	auto points = c_.get_orthogonal_points();
	std::vector<double> weights = c_.get_vdW_weights();
	for (auto point : points) {
		wp.emplace_back(point.first, point.second);
	}
	//for (int i = 0; i < lp.size(); i++) {
		//wp.emplace_back(Point_3(points[i].x, points[i].y, points[i].z), weights[i]);
	//}
	// alpha shape construction from points. CGAL has a strange behavior in REGULARIZED mode.
	//Alpha_shape_3 as(wp.begin(), wp.end(), 0, Alpha_shape_3::GENERAL);
	as_->set_mode(Alpha_shape_3::GENERAL);
	as_->set_alpha(0);
	as_->make_alpha_shape(wp.begin(), wp.end());
	
	for (auto i = as_->alpha_begin(); i != as_->alpha_end(); ++i) {
		alpha_vals_.push_back(*i);
	}

	// filtration with alpha values from alpha shape
	std::vector<Object> the_objects;
	std::vector<Alpha_value_type> the_alpha_values;
	Dispatch disp = CGAL::dispatch_output<Object, Alpha_value_type>(std::back_inserter(the_objects),
		std::back_inserter(the_alpha_values));
	as_->filtration_with_alpha_values(disp);
	Alpha_shape_3::size_type count_vertices = 0;
	Alpha_shape_3::size_type count_edges = 0;
	Alpha_shape_3::size_type count_facets = 0;
	Alpha_shape_3::size_type count_cells = 0;
	// Loop on objects vector
	Vertex_list vertex_list;
	ST simplex_tree;
	Alpha_shape_simplex_tree_map map_cgal_simplex_tree;
	std::vector<Alpha_value_type>::iterator the_alpha_value_iterator = the_alpha_values.begin();
	for (auto object_iterator : the_objects) {
		// Retrieve Alpha shape vertex list from object
		if (const Cell_handle *cell = CGAL::object_cast<Cell_handle>(&object_iterator)) {
			vertex_list = from_cell<Vertex_list, Cell_handle>(*cell);
			count_cells++;
		}
		else if (const Facet *facet = CGAL::object_cast<Facet>(&object_iterator)) {
			vertex_list = from_facet<Vertex_list, Facet>(*facet);
			count_facets++;
		}
		else if (const Edge_3 *edge = CGAL::object_cast<Edge_3>(&object_iterator)) {
			vertex_list = from_edge<Vertex_list, Edge_3>(*edge);
			count_edges++;
		}
		else if (const Vertex_handle *vertex = CGAL::object_cast<Vertex_handle>(&object_iterator)) {
			count_vertices++;
			vertex_list = from_vertex<Vertex_list, Vertex_handle>(*vertex);
		}
		// Construction of the vector of simplex_tree vertex from list of alpha_shapes vertex
		Simplex_tree_vector_vertex the_simplex;
		for (auto the_alpha_shape_vertex : vertex_list) {
			Alpha_shape_simplex_tree_map::iterator the_map_iterator = map_cgal_simplex_tree.find(the_alpha_shape_vertex);
			if (the_map_iterator == map_cgal_simplex_tree.end()) {
				// alpha shape not found
				Simplex_tree_vertex vertex = map_cgal_simplex_tree.size();
				the_simplex.push_back(vertex);
				map_cgal_simplex_tree.emplace(the_alpha_shape_vertex, vertex);
			}
			else {
				// alpha shape found
				Simplex_tree_vertex vertex = the_map_iterator->second;
				the_simplex.push_back(vertex);
			}
		}
		// Construction of the simplex_tree
		Filtration_value filtr = /*std::sqrt*/ (*the_alpha_value_iterator);
		simplex_tree.insert_simplex(the_simplex, filtr);
		if (the_alpha_value_iterator != the_alpha_values.end())
			++the_alpha_value_iterator;
		else
			std::cout << "This shall not happen" << std::endl;
	}
	// Sort the simplices in the order of the filtration
	simplex_tree.initialize_filtration();
	// std::cout << "Simplex_tree dim: " << simplex_tree.dimension() << std::endl;
	// Compute the persistence diagram of the complex
	Persistent_cohomology pcoh(simplex_tree, true);
	// initializes the coefficient field for homology
	pcoh.init_coefficients(coeff_field_characteristic);
	pcoh.compute_persistent_cohomology(min_persistence);
	// Output the diagram in filediag

	cmp_intervals_by_dim_then_length cmp(&simplex_tree);
	auto pairs = pcoh.get_persistent_pairs();
	//pcoh.

	std::sort(std::begin(pairs), std::end(pairs), cmp);

	// Save pairs for later
	for (auto pair : pairs) {
		if (simplex_tree.dimension(get<0>(pair)) == 1) {
			pairs_.emplace_back(simplex_tree.filtration(get<0>(pair)), simplex_tree.filtration(get<1>(pair)));
		}
	}
}

void Persistence::calculate_clusters(double tolerance)
{
	pairs_vals_.clear();
	clusters_.clear();

	for (auto pair : pairs_) {
		pairs_vals_.push_back(pair.second - pair.first);
	}

	std::vector<double> pers_diff;
	int pers_diag_gap = -1;
	int pers_diag_gap_index = -1;

	pers_diff.push_back(0);
	for (int i = 0; i < pairs_vals_.size() - 1; i++) {
		pers_diff.push_back(pairs_vals_[i] - pairs_vals_[i + 1]);
		if (pers_diff[i] > pers_diag_gap) {
			pers_diag_gap = pers_diff[i];
			pers_diag_gap_index = i;
		}
	}

	for (int i = 0; i <= pers_diag_gap_index; i++) {
		int num_points = 0;
		while (pers_diff[i] < tolerance) {
			i++;
			num_points++;
		}

		if (num_points > 0) {
			clusters_.push_back(num_points);
		}
	}
}

void Persistence::output_complex(const std::string& path)
{
	// Write the file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(path.c_str());
	writer->SetInputData(complex_source_);
	writer->Write();
}

void Persistence::output_pairs(const std::string& path)
{
	std::ofstream out(path);
	std::cout << "Result in file: " << path << std::endl;

	for (auto pair : pairs_) {
		auto fc = pair.first;
		auto sc = pair.second;

		out << std::fixed << std::setprecision(6) << "B: "
			<< fc << "\t\tD: "
			<< sc << "\t\tD-B: " << sc - fc << std::endl;
	}

	out.close();
}

void Persistence::output(const std::string& path)
{
	//PAIR - PAIR_BD - PAIR HOLE IN CLUSTER - PAIR CLUSTERS
	std::ofstream out(path, std::ofstream::out | std::ofstream::app);

	out << pairs_[0].first << "," << pairs_[0].second << "\t\t" <<
		pairs_vals_[0] << "\t\t" <<
		clusters_.size() << "\t\t";
	for (auto &cluster_point : clusters_)
		out << cluster_point << "\t\t";
	out << std::endl;

	for (auto &pair : pairs_) {
		out << pair.first << " " << pair.second << std::endl;
	}

	out.close();
}

void Persistence::update_point_map(Alpha_shape_3 * alpha_shape) {
	point_map_.clear();

	//std::vector<vtkIdType> point_ids;

	vtkIdType id = 0;
	for (auto i = alpha_shape->points_begin(); i != alpha_shape->points_end(); ++i) {
		//auto p = i->point();

		if (point_map_.find(*i) == point_map_.end()) {
			point_map_.emplace(*i, id++);
		}
	}
}

vtkSmartPointer<vtkPoints> Persistence::get_alpha_points(Alpha_shape_3 * alpha_shape)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkIdType id = 0;
	for (auto i = alpha_shape->points_begin(); i != alpha_shape->points_end(); ++i) {
		//auto p = *i;

		if (point_map_.find(*i) != point_map_.end()) {
			points->InsertNextPoint(i->x(), i->y(), i->z());
		}
	}

	return points;
}

vtkSmartPointer<vtkCellArray> Persistence::get_alpha_edges(Alpha_shape_3 * alpha_shape)
{
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	std::vector<vtkSmartPointer<vtkLine>> poly_lines;

	std::vector<Edge_3> edges;
	alpha_shape->get_alpha_shape_edges(std::back_inserter(edges), Alpha_shape_3::SINGULAR, alpha_shape->get_alpha());
	for (auto &edge : edges) {
		if (alpha_shape->is_infinite(edge))
			continue;

		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
		Vertex_list vertex_list = from_edge<Vertex_list, Edge_3>(edge);
		std::vector<vtkIdType> line_ids;

		for (auto v : vertex_list) {
			auto r = point_map_.find(v->point());
			if (r != point_map_.end()) {
				line_ids.push_back(r->second);
			}
		}

		line->GetPointIds()->SetNumberOfIds(line_ids.size());

		for (int i = 0; i < line_ids.size(); i++) {
			line->GetPointIds()->SetId(i, line_ids[i]);
		}

		//edges_ids.push_back(line_ids);
		poly_lines.push_back(line);
	}

	for (auto &line : poly_lines) {
		lines->InsertNextCell(line);
	}

	return lines;
}

vtkSmartPointer<vtkCellArray> Persistence::get_alpha_faces(Alpha_shape_3* alpha_shape) {
	vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
	std::vector<vtkSmartPointer<vtkPolygon>> poly_faces;

	for (auto i = alpha_shape->alpha_shape_facets_begin(); i != alpha_shape->alpha_shape_facets_end(); ++i) {
		if (alpha_shape->is_infinite(*i))
			continue;

		vtkSmartPointer<vtkPolygon> face = vtkSmartPointer<vtkPolygon>::New();
		Vertex_list vertex_list = from_facet<Vertex_list, Facet>(*i);
		std::vector<vtkIdType> face_ids;

		for (auto v : vertex_list) {
			auto r = point_map_.find(v->point());
			if (r != point_map_.end()) {
				face_ids.push_back(r->second);
			}
		}

		face->GetPointIds()->SetNumberOfIds(face_ids.size());

		for (int i = 0; i < face_ids.size(); i++) {
			face->GetPointIds()->SetId(i, face_ids[i]);
		}

		//faces_ids.push_back(face_ids);
		poly_faces.push_back(face);
	}

	for (auto &poly_face : poly_faces) {
		faces->InsertNextCell(poly_face);
	}

	return faces;
}

vtkSmartPointer<vtkPolyData> Persistence::update_poly_data(Alpha_shape_3 *alpha_shape, double alpha)
{
	//double sr = alpha;// sqrt(0);
	alpha_shape->set_alpha(alpha);

	const auto points = get_alpha_points(alpha_shape);
	const auto lines = get_alpha_edges(alpha_shape);
	const auto faces = get_alpha_faces(alpha_shape);

	// Create a polydata object and add the points to it.
	complex_source_ = vtkSmartPointer<vtkPolyData>::New();
	complex_source_->SetPoints(points);
	complex_source_->SetLines(lines);
	complex_source_->SetPolys(faces);

	return complex_source_;
}

void Persistence::init_win() {
	win_ = vtkSmartPointer<vtkRenderWindow>::New();
	win_->SetWindowName("complex");
	win_->SetSize(1280, 720);
	//renwin->AddRenderer(renderer);
}

void Persistence::init_ren()
{
	ren_ = vtkSmartPointer<vtkRenderer>::New();

	win_->AddRenderer(ren_);
	
	inter_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	inter_->SetRenderWindow(win_);
}

void Persistence::init_complex_actor() {
	vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(complex_source_);

	complex_actor_ = vtkSmartPointer<vtkActor>::New();
	complex_actor_->SetMapper(mapper);
	complex_actor_->GetProperty()->SetColor(colors->GetColor3d("Silver").GetData());

	ren_->AddActor(complex_actor_);
}

void Persistence::init_point_sphere_actors() {
	//std::vector<vtkSmartPointer<vtkSphereSource>> sphere_sources;
	std::vector<vtkSmartPointer<vtkPolyDataMapper>> sphere_mappers;
	std::vector<vtkSmartPointer<vtkActor>> sphere_actors;

	//point_sphere_actors_.clear();
	//point_sphere_sources_.clear();

	auto points = c_.get_orthogonal_points();
	for (auto& point : points) {
		vtkSmartPointer<vtkSphereSource> sphere_source = vtkSmartPointer<vtkSphereSource>::New();
		vtkSmartPointer<vtkPolyDataMapper> sphere_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkActor> sphere_actor = vtkSmartPointer<vtkActor>::New();

		double x = point.first[0];
		double y = point.first[1];
		double z = point.first[2];
		double w = point.second;

		sphere_mapper->SetInputConnection(sphere_source->GetOutputPort());
		sphere_actor->SetMapper(sphere_mapper);

		sphere_source->SetCenter(x, y, z);
		sphere_source->SetRadius(0.1);
		
		sphere_source->SetPhiResolution(12);
		sphere_source->SetThetaResolution(12);

		sphere_actor->GetProperty()->SetColor(255, 0, 0);
		sphere_actor->GetProperty()->SetOpacity(1);

		//sphere_sources.push_back(sphere_source);

		point_sphere_sources_.push_back(sphere_source);
		sphere_mappers.push_back(sphere_mapper);
		sphere_actors.push_back(sphere_actor);
		point_sphere_actors_.push_back(sphere_actor);
	}

	for (auto &sphere_actor : sphere_actors) {
		ren_->AddActor(sphere_actor);
	}
}

void Persistence::init_radii_sphere_actors() {
	//std::vector<vtkSmartPointer<vtkSphereSource>> sphere_sources;
	//std::vector<std::pair<vtkSmartPointer<vtkSphereSource>, double>> weighted_sphere_sources;
	std::vector<vtkSmartPointer<vtkPolyDataMapper>> sphere_mappers;
	std::vector<vtkSmartPointer<vtkActor>> sphere_actors;

	//radii_sphere_actors_.clear();
	//radii_sphere_sources_.clear();

	auto points = c_.get_orthogonal_points();
	for (auto& point : points) {
		vtkSmartPointer<vtkSphereSource> sphere_source = vtkSmartPointer<vtkSphereSource>::New();
		vtkSmartPointer<vtkPolyDataMapper> sphere_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkActor> sphere_actor = vtkSmartPointer<vtkActor>::New();

		double x = point.first[0];
		double y = point.first[1];
		double z = point.first[2];
		double w = point.second;

		sphere_mapper->SetInputConnection(sphere_source->GetOutputPort());
		sphere_actor->SetMapper(sphere_mapper);

		sphere_source->SetCenter(x, y, z);
		sphere_source->SetRadius(sqrt(w));

		sphere_source->SetPhiResolution(12);
		sphere_source->SetThetaResolution(12);

		sphere_actor->GetProperty()->SetColor(0, 255, 0);
		sphere_actor->GetProperty()->SetOpacity(0.05);

		//sphere_sources.push_back(sphere_source);
		radii_sphere_sources_.emplace_back(sphere_source, w);

		sphere_mappers.push_back(sphere_mapper);
		sphere_actors.push_back(sphere_actor);
		radii_sphere_actors_.push_back(sphere_actor);
	}

	for (auto &sphere_actor : sphere_actors) {
		ren_->AddActor(sphere_actor);
	}
}

void Persistence::init_actors() {
	init_complex_actor();
	init_point_sphere_actors();
	init_radii_sphere_actors();
}

void Persistence::init_slider() {
	vtkSmartPointer<vtkSliderRepresentation2D> slider_rep = vtkSmartPointer<vtkSliderRepresentation2D>::New();

	double max_alpha = 0;
	double min_alpha = std::numeric_limits<double>::max();
	for (auto i = as_->alpha_begin(); i != as_->alpha_end(); ++i) {
		if (*i > max_alpha) {
			max_alpha = *i;
		}
		if (*i < min_alpha) {
			min_alpha = *i;
		}
	}

	double max_death = 0;
	for (auto pair : pairs_) {
		if (pair.second > max_death)
			max_death = pair.second;
	}

	slider_rep->SetMinimumValue(min_alpha);
	slider_rep->SetMaximumValue(sqrt(max_death + 1));
	slider_rep->SetValue(0);
	//sliderRep->SetValue(alpha_vals[0]);

	//sliderRep->SetTitleText("Alpha value");
	slider_rep->GetSliderProperty()->SetColor(1, 0, 0);//red
	slider_rep->GetTitleProperty()->SetColor(1, 0, 0);//red
	slider_rep->GetLabelProperty()->SetColor(1, 0, 0);//red
	slider_rep->GetSelectedProperty()->SetColor(0, 1, 0);//green
	slider_rep->GetTubeProperty()->SetColor(1, 1, 0);//yellow

	// Change the color of the ends of the bar
	slider_rep->GetCapProperty()->SetColor(1, 1, 0);//yellow
	slider_rep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
	slider_rep->GetPoint1Coordinate()->SetValue(40, 40);
	slider_rep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
	slider_rep->GetPoint2Coordinate()->SetValue(740, 40);

	slider_widget_ = vtkSmartPointer<vtkSliderWidget>::New();
	slider_widget_->SetInteractor(inter_);
	slider_widget_->SetRepresentation(slider_rep);
	slider_widget_->SetAnimationModeToAnimate();
	slider_widget_->EnabledOn();
}

void Persistence::init_widgets() {
	init_slider();
}

void Persistence::init_slider_callback() {
	vtkSmartPointer<vtkSliderCallback> callback = vtkSmartPointer<vtkSliderCallback>::New();

	callback->p = this;
	//callback->poly_data = complex_source_;
	//callback->spheres = sphere_sources;
	//callback->weighted_spheres = radii_sphere_sources_;
	//callback->text = text_ptr;

	slider_widget_->AddObserver(vtkCommand::EndInteractionEvent, callback);
}

void Persistence::init_update_callback()
{
	vtkSmartPointer<vtkUpdateCallback> update_callback = vtkSmartPointer<vtkUpdateCallback>::New();
	update_callback->p = this;
	update_callback->poly_data = complex_source_;
	update_callback->weighted_spheres = radii_sphere_sources_;
	static const vtkCommand::EventIds interactorEvents[] = {
		vtkCommand::MouseMoveEvent,
		vtkCommand::LeftButtonPressEvent,
		vtkCommand::LeftButtonReleaseEvent,
		vtkCommand::RightButtonPressEvent,
		vtkCommand::RightButtonReleaseEvent,
		vtkCommand::MiddleButtonPressEvent,
		vtkCommand::MiddleButtonReleaseEvent,
		vtkCommand::MouseWheelForwardEvent,
		vtkCommand::MouseWheelBackwardEvent,
		vtkCommand::KeyPressEvent
	};
	for (auto e : interactorEvents) {
		inter_->AddObserver(e, update_callback);
	}
}

void Persistence::init_atb_callback() {
	atb_callback_ = vtkSmartPointer<vtkAntTweakBar>::New();
	atb_callback_->Initialize(inter_);
	atb_callback_->EnableAlwaysRenderOnEvent();
}

void Persistence::init_callbacks() {
	init_slider_callback();
	init_update_callback();
	init_atb_callback();
}

void TW_CALL SetBgColor(const void *value, void *clientData) {
	auto actor = static_cast<Persistence*>(clientData)->GetRenderer();
	auto color = static_cast<const float*>(value);
	actor->SetBackground(color[0], color[1], color[2]);
}

void TW_CALL GetBgColor(void *value, void *clientData) {
	auto actor = static_cast<Persistence*>(clientData)->GetRenderer();
	auto color = static_cast<float*>(value);
	auto actorColor = actor->GetBackground();
	for (int i = 0; i < 3; i++)
		color[i] = (float)actorColor[i];
}

void TW_CALL SetColor(const void *value, void *clientData) {
	auto actor = static_cast<Persistence*>(clientData)->GetComplexActor();
	auto color = static_cast<const float*>(value);
	actor->GetProperty()->SetColor(color[0], color[1], color[2]);
}

void TW_CALL GetColor(void *value, void *clientData) {
	auto actor = static_cast<Persistence*>(clientData)->GetComplexActor();
	auto color = static_cast<float*>(value);
	auto actorColor = actor->GetProperty()->GetColor();
	for (int i = 0; i < 3; i++)
		color[i] = (float)actorColor[i];
}

void TW_CALL SetAlpha(const void *value, void *clientData) {
	auto p = static_cast<Persistence*>(clientData);
	auto alpha = static_cast<const double*>(value);
	//actor->GetProperty()->SetColor(color[0], color[1], color[2]);
	p->set_alpha(*alpha);
}

void TW_CALL GetAlpha(void *value, void *clientData) {
	auto p = static_cast<Persistence*>(clientData);
	*static_cast<double*>(value) = p->get_alpha();
	//auto actorColor = actor->GetProperty()->GetColor();
	//for (int i = 0; i < 3; i++)
	//	color[i] = (float)actorColor[i];
}

void TW_CALL SetDisplayAlpha(const void *value, void *clientData) {
	auto p = static_cast<Persistence*>(clientData);
	auto alpha = static_cast<const double*>(value);
	//actor->GetProperty()->SetColor(color[0], color[1], color[2]);
	p->set_alpha(*alpha);
}

void TW_CALL GetDisplayAlpha(void *value, void *clientData) {
	auto p = static_cast<Persistence*>(clientData);
	*static_cast<double*>(value) = sqrt(p->get_alpha());
	//auto actorColor = actor->GetProperty()->GetColor();
	//for (int i = 0; i < 3; i++)
	//	color[i] = (float)actorColor[i];
}

void TW_CALL SetSphColor(const void *value, void *clientData) {
	auto actor = static_cast<Persistence*>(clientData)->GetRadiiSpheresActors();
	auto color = static_cast<const float*>(value);
	for (auto& a : actor) {
		a->GetProperty()->SetColor(color[0], color[1], color[2]);
	}
}

void TW_CALL GetSphColor(void *value, void *clientData) {
	auto actor = static_cast<Persistence*>(clientData)->GetRadiiSpheresActors();
	auto color = static_cast<float*>(value);
	auto actorColor = actor.at(0)->GetProperty()->GetColor();
	for (int i = 0; i < 3; i++)
		color[i] = (float)actorColor[i];
}

void TW_CALL SetSphOpacity(const void *value, void *clientData) {
	auto actor = static_cast<Persistence*>(clientData)->GetRadiiSpheresActors();
	auto color = *static_cast<const double*>(value);
	for (auto& a : actor) {
		a->GetProperty()->SetOpacity(color);
	}
}

void TW_CALL GetSphOpacity(void *value, void *clientData) {
	auto actor = static_cast<Persistence*>(clientData)->GetRadiiSpheresActors();
	*static_cast<double*>(value) = actor.at(0)->GetProperty()->GetOpacity();
}

void TW_CALL SetAlwaysRenderOnEvent(const void *value, void *clientData) {
	auto callback = static_cast<vtkAntTweakBar*>(clientData);
	bool alwaysRender = *static_cast<const bool*>(value);
	if (alwaysRender)
		callback->EnableAlwaysRenderOnEvent();
	else
		callback->DisableAlwaysRenderOnEvent();
}

void TW_CALL GetAlwaysRenderOnEvent(void *value, void *clientData) {
	auto callback = static_cast<vtkAntTweakBar*>(clientData);
	*static_cast<bool*>(value) = callback->AlwaysRenderOnEvent();
}

void Persistence::display_complex() {
	update_point_map(as_);
	update_poly_data(as_);

	init_win();
	init_ren();
	init_actors();
	init_widgets();

	auto style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	inter_->SetInteractorStyle(style);

	ren_->ResetCamera();

	init_callbacks();

	vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

	ren_->SetBackground(colors->GetColor3d("Salmon").GetData());

	//inter_->Initialize();
	
	TwBar* bar = TwNewBar("Options");
	TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with VTK.' ");

	TwAddVarCB(bar, "background color", TW_TYPE_COLOR3F, SetBgColor, GetBgColor, this,
		" help='Change the color of the complex' ");
	TwAddVarCB(bar, "complex color", TW_TYPE_COLOR3F, SetColor, GetColor, this,
		" help='Change the color of the complex' ");
	TwAddSeparator(bar, "", "");
	TwAddVarCB(bar, "alpha", TW_TYPE_DOUBLE, SetAlpha, GetAlpha, this,
		" help='Change the value of alpha' ");
	TwAddVarCB(bar, "display alpha", TW_TYPE_DOUBLE, SetDisplayAlpha, GetDisplayAlpha, this,
		" help='Change the value of display alpha' ");
	TwAddSeparator(bar, "", "");
	TwAddVarCB(bar, "sphere color", TW_TYPE_COLOR3F, SetSphColor, GetSphColor, this,
		" help='Change the color of the complex' ");
	TwAddVarCB(bar, "sphere opacity", TW_TYPE_DOUBLE, SetSphOpacity, GetSphOpacity, this,
		" help='Change the value of alpha' min=0 max=1 step=0.05");
	TwAddSeparator(bar, "", "");
	TwAddVarCB(bar, "render on event", TW_TYPE_BOOLCPP, SetAlwaysRenderOnEvent, GetAlwaysRenderOnEvent, atb_callback_,
		" true='Always' false='If handled' help=`Specify vtkAntTweakBar's behaviour.` ");

	//inter_->AddObserver(vtkCommand::AnyEvent, vtkUpdateCallback);

	//win_->Render();
	inter_->Start();

	//std::cout << "really?";
}

std::vector<std::pair<double, double>> Persistence::get_pairs() const{ return pairs_; }

std::vector<double> Persistence::get_alpha_values() const { return alpha_vals_; }

std::map<Weighted_point_3, int> Persistence::get_point_map() const { return point_map_; }

std::vector<Weighted_point_3> Persistence::get_weighted_points() const { return points_; }

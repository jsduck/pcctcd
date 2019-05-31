#include "vtkUpdateCallback.h"

#include "persistence.h"

void vtkUpdateCallback::Execute(vtkObject* caller, unsigned long, void*) {
	auto ptr = static_cast<Persistence*>(p);
	if (ptr->update_flag) {
		ptr->update_poly_data(ptr->as_, ptr->alpha_);
		auto pd = ptr->complex_source_;
		poly_data->SetPoints(pd->GetPoints());
		poly_data->SetPolys(pd->GetPolys());
		poly_data->SetLines(pd->GetLines());

		for (auto& sphere : weighted_spheres) {
			sphere.first->SetRadius(sqrt(ptr->alpha_ + sphere.second));
		}

		ptr->update_flag = false;
	}
}

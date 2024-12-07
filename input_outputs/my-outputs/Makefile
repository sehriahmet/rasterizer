all:
	g++ ../../code_template/*.cpp -O3 -o rasterizer -std=c++17 -w

	./rasterizer ../clipping_example/empty_box_clipped.xml

	./rasterizer ../culling_disabled_inputs/empty_box.xml 
	./rasterizer ../culling_disabled_inputs/filled_box.xml
	./rasterizer ../culling_disabled_inputs/sample.xml
	./rasterizer ../culling_disabled_inputs/flag_czechia_alternative.xml # not correct
	./rasterizer ../culling_disabled_inputs/flag_czechia.xml
	./rasterizer ../culling_disabled_inputs/flag_germany.xml
	./rasterizer ../culling_disabled_inputs/flag_turkey_alternative.xml # not correct 
	./rasterizer ../culling_disabled_inputs/flag_turkey.xml
	./rasterizer ../culling_disabled_inputs/horse_and_mug.xml

	./rasterizer ../culling_enabled_inputs/empty_box.xml
	./rasterizer ../culling_enabled_inputs/filled_box.xml
	./rasterizer ../culling_enabled_inputs/sample.xml
	./rasterizer ../culling_enabled_inputs/flag_czechia_alternative.xml
	./rasterizer ../culling_enabled_inputs/flag_czechia.xml
	./rasterizer ../culling_enabled_inputs/flag_iceland.xml
	./rasterizer ../culling_enabled_inputs/flag_germany.xml
	./rasterizer ../culling_enabled_inputs/flag_turkey_alternative.xml
	./rasterizer ../culling_enabled_inputs/flag_turkey.xml
	./rasterizer ../culling_enabled_inputs/horse_and_mug.xml

	./rasterizer ../different_projection_type/flag_turkey/flag_turkey_orthographic.xml
	./rasterizer ../different_projection_type/flag_turkey/flag_turkey_perspective.xml

	./rasterizer ../different_projection_type/horse_and_mug/horse_and_mug_orthographic.xml
	./rasterizer ../different_projection_type/horse_and_mug/horse_and_mug_perspective.xml
	
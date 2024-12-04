#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"

#include <vector>

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName, int osType);
	void forwardRenderingPipeline(Camera *camera, Scene *scene);

	// in order to use this-> this is included in here
	std::vector<Vec3> clipLine(Vec3 &vertex1, Vec3 &vertex2, Camera *camera);
	void drawLine(Vec3 v0, Vec3 v1, std::vector<std::vector<double>> &depthBuffer, Camera *camera);
	void rasterizeEdges(std::vector<Vec3> transformedVertices, const Triangle &triangle, std::vector<std::vector<double>> &depthBuffer, Camera *camera);
	void rasterizeFilledTriangle(Triangle &triangle, std::vector<Vec3> &transformedVertices, Camera *camera, std::vector<std::vector<double>> &depthBuffer);

};

#endif

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"


#include <iostream>

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/

Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}


/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}


/*
	Transformations, clipping, culling, rasterization are done here.
*/

void modelingTransformation(Vec3 &vertices1, int colorId, const vector<char> &transformationTypes, const vector<int> &transformationIds, Scene *scene, Mesh *mesh) {
    Vec4 triangleHc = { vertices1.x, vertices1.y, vertices1.z, 1 };
	int i = 1000000000;

	for (auto &triangle : mesh->triangles) {
		if (colorId == triangle.vertexIds[0] || colorId == triangle.vertexIds[1] || colorId == triangle.vertexIds[2] ) {
			i = 0;
		} 
	}  

    for (i ; i < transformationTypes.size(); i++) {
        char type = transformationTypes[i];
        int id = transformationIds[i];    
    	
        if (type == 't') {
    		for (Translation *trans : scene->translations) {				
    			Matrix4 identityMatrix = getIdentityMatrix();
    			if (trans->translationId == id) {
    				// printf("%f %f %f %f %d\n", triangleHc.x, triangleHc.y, triangleHc.z, triangleHc.t, id);
    				identityMatrix.values[0][3] += trans->tx;
    				identityMatrix.values[1][3] += trans->ty;
    				identityMatrix.values[2][3] += trans->tz;
    				triangleHc = multiplyMatrixWithVec4(identityMatrix, triangleHc);
    			}
    		}
        } else if (type == 's') {
    		for (Scaling *scale : scene->scalings) {				
    			Matrix4 identityMatrix = getIdentityMatrix();
    			if (scale->scalingId == id) {
    				// printf("%f %f %f %f %d\n", triangleHc.x, triangleHc.y, triangleHc.z, triangleHc.t, id);
    				identityMatrix.values[0][0] = scale->sx;
    				identityMatrix.values[1][1] = scale->sy;
    				identityMatrix.values[2][2] = scale->sz;
    				triangleHc = multiplyMatrixWithVec4(identityMatrix, triangleHc);
    			}
    		}    
        } else if (type == 'r') {
    		for (Rotation *rotate : scene->rotations) {				
    			Matrix4 rotationMatrix = getIdentityMatrix();
    			if (rotate->rotationId == id) {
    				double radCos = rotate->angle * M_PI / 180;
    				rotationMatrix.values[0][0] = cos(radCos) + rotate->ux * rotate->uy * (1 - cos(radCos));
                    rotationMatrix.values[0][1] = rotate->ux * rotate->uy* (1 - cos(radCos)) - rotate->uz * sin(radCos);
                    rotationMatrix.values[0][2] = rotate->ux  * rotate->uz * (1 - cos(radCos)) + rotate->uy * sin(radCos);
                    rotationMatrix.values[1][0] = rotate->uy * rotate->ux * (1 - cos(radCos))+ rotate->uz * sin(radCos);
                    rotationMatrix.values[1][1] = cos(radCos)  + rotate->uy * rotate->uy* (1 - cos(radCos));
                    rotationMatrix.values[1][2] = rotate->uy * rotate->uz * (1 - cos(radCos)) - rotate->ux  * sin(radCos);
                    rotationMatrix.values[2][0] = rotate->uz * rotate->ux * (1 - cos(radCos)) - rotate->uy * sin(radCos);
                    rotationMatrix.values[2][1] = rotate->uz * rotate->uy * (1 - cos(radCos)) + rotate->ux  * sin(radCos);
                    rotationMatrix.values[2][2] = cos(radCos) + rotate->uz * rotate->uz * (1 - cos(radCos));
    				triangleHc = multiplyMatrixWithVec4(rotationMatrix, triangleHc);
    				// cout << cos(radCos) << endl;
    			}
    		}
        }
    }  
	vertices1.x = triangleHc.x;
	vertices1.y = triangleHc.y;
	vertices1.z = triangleHc.z;
	
	// cout << vertices1.colorId << "   " <<vertices1.x << "   " << vertices1.y << "   "  << vertices1.z << endl;

}

// kamerada gaze ve up vectorleri bazen dik gelmiyor ornegin: emptyboxclipped 1. olanindan.. !!
void viewTransformations(Vec3 &vertices1, Camera *camera) {
    Matrix4 viewMatrix;
	Vec4 triangleHc = { vertices1.x, vertices1.y, vertices1.z, 1 };

	Vec3 w = camera->w;
	Vec3 v = camera->v;
	Vec3 u = camera->u;

	viewMatrix.values[0][0] = u.x;
	viewMatrix.values[0][1] = u.y;
	viewMatrix.values[0][2] = u.z;
	viewMatrix.values[0][3] = -(u.x * camera->position.x + u.y * camera->position.y + u.z * camera->position.z);
	viewMatrix.values[1][0] = v.x;
	viewMatrix.values[1][1] = v.y;
	viewMatrix.values[1][2] = v.z;
	viewMatrix.values[1][3] = -(v.x * camera->position.x + v.y * camera->position.y + v.z * camera->position.z);
	viewMatrix.values[2][0] = w.x;
	viewMatrix.values[2][1] = w.y;
	viewMatrix.values[2][2] = w.z;
	viewMatrix.values[2][3] = -(w.x * camera->position.x + w.y * camera->position.y + w.z * camera->position.z);
	viewMatrix.values[3][0] = 0;
	viewMatrix.values[3][1] = 0;
	viewMatrix.values[3][2] = 0;
	viewMatrix.values[3][3] = 1;

	triangleHc = multiplyMatrixWithVec4(viewMatrix, triangleHc);
	
	vertices1.x = triangleHc.x;
	vertices1.y = triangleHc.y;
	vertices1.z = triangleHc.z;
    // cout << vertices1.x << "   " << vertices1.y << "   "  << vertices1.z << endl;
}


void projectTransformations(Vec3 &vertices1, Camera *camera) {
	Matrix4 projectionMatrix;
	Vec4 triangleHc = { vertices1.x, vertices1.y, vertices1.z, 1 };

	// cout << "projection type: " << camera->projectionType <<"   "<<  
	//		(camera->right - camera->left)<< "  " <<
	// 		(camera->top - camera->bottom)<< "   " << 
	// 		(camera->far - camera->near)<< endl;
	if (camera->projectionType == 0) {
		projectionMatrix.values[0][0] = 2/(camera->right - camera->left);
		projectionMatrix.values[0][1] = 0;
		projectionMatrix.values[0][2] = 0;
		projectionMatrix.values[0][3] = -((camera->right + camera->left)/ (camera->right - camera->left));
		projectionMatrix.values[1][0] = 0;
		projectionMatrix.values[1][1] = 2/(camera->top - camera->bottom);
		projectionMatrix.values[1][2] = 0;
		projectionMatrix.values[1][3] = -((camera->top + camera->bottom)/ (camera->top - camera->bottom));
		projectionMatrix.values[2][0] = 0;
		projectionMatrix.values[2][1] = 0;
		projectionMatrix.values[2][2] = -(2/(camera->far - camera->near));
		projectionMatrix.values[2][3] = -((camera->far + camera->near)/ (camera->far - camera->near));
		projectionMatrix.values[3][0] = 0;
		projectionMatrix.values[3][1] = 0;
		projectionMatrix.values[3][2] = 0;
		projectionMatrix.values[3][3] = 1;
	} else if (camera->projectionType == 1) {
		// ?? metenin dyligindan soru isareti
		projectionMatrix.values[0][0] = (2*camera->near) / (camera->right - camera->left);
		projectionMatrix.values[0][1] = 0;
		projectionMatrix.values[0][2] = (camera->right + camera->left)/ (camera->right - camera->left);
		projectionMatrix.values[0][3] = 0;
		projectionMatrix.values[1][0] = 0;
		projectionMatrix.values[1][1] = (2*camera->near) / (camera->top - camera->bottom);
		projectionMatrix.values[1][2] = (camera->top + camera->bottom)/ (camera->top - camera->bottom);
		projectionMatrix.values[1][3] = 0;
		projectionMatrix.values[2][0] = 0;
		projectionMatrix.values[2][1] = 0;
		projectionMatrix.values[2][2] = -((camera->far + camera->near)/ (camera->far - camera->near));
		projectionMatrix.values[2][3] = -(2 * camera->near * camera->far) / (camera->far - camera->near);
		projectionMatrix.values[3][0] = 0;
		projectionMatrix.values[3][1] = 0;
		projectionMatrix.values[3][2] = -1;
		projectionMatrix.values[3][3] = 0;
	} else {
		printf("something wrong in project transformation");
	}

	triangleHc = multiplyMatrixWithVec4(projectionMatrix, triangleHc);
	
	vertices1.x = triangleHc.x;
	vertices1.y = triangleHc.y;
	vertices1.z = triangleHc.z;
    // cout << vertices1.x << "   " << vertices1.y << "   "  << vertices1.z << endl;
}

bool isBackface(Triangle &triangle, vector<Vec3> &Vertices, Camera *camera) {
    Vec3 v0 = Vertices[triangle.vertexIds[0] - 1];
    Vec3 v1 = Vertices[triangle.vertexIds[1] - 1];
    Vec3 v2 = Vertices[triangle.vertexIds[2] - 1];
    //cout << Vertices[triangle.vertexIds[0] - 1] << endl;
    Vec3 normal = crossProductVec3(subtractVec3(v1, v0), subtractVec3(v2, v0));

    Vec3 viewDir = subtractVec3(camera->position, v0);

    return dotProductVec3(normalizeVec3(normal), normalizeVec3(viewDir)) <= 0;
}

/*
vector<Vec3> Scene::clipAgainstPlane(const vector<Vec3> &vertices, const Vec3 &planeNormal, double planeOffset) {
    vector<Vec3> clippedVertices;

    for (size_t i = 0; i < vertices.size(); i++) {
        const Vec3 &currentVertex = vertices[i];
        const Vec3 &nextVertex = vertices[(i + 1) % vertices.size()];

        double currentDist = dotProductVec3(planeNormal, currentVertex) - planeOffset;
        double nextDist = dotProductVec3(planeNormal, nextVertex) - planeOffset;

        // Case 1: Current vertex is inside the plane
        if (currentDist >= 0) {
            clippedVertices.push_back(currentVertex);
        }

        // Case 2: Edge intersects the plane
        if ((currentDist >= 0) != (nextDist >= 0)) {  // XOR: one inside, one outside
            double t = currentDist / (currentDist - nextDist);  // Interpolation factor
            Vec3 intersection;
			intersection.x = currentVertex.x + (nextVertex.x - currentVertex.x) * t;
			intersection.y = currentVertex.y + (nextVertex.y - currentVertex.y) * t;
			intersection.z = currentVertex.z + (nextVertex.z - currentVertex.z) * t;
            clippedVertices.push_back(intersection);
        }
    }

    return clippedVertices;
}

vector<Vec3> Scene::clipTriangle(const Triangle &triangle, const vector<Vec3> &vertices, Camera *camera) {
    // Start with the original triangle vertices
    vector<Vec3> currentVertices = {
        vertices[triangle.vertexIds[0] - 1],
        vertices[triangle.vertexIds[1] - 1],
        vertices[triangle.vertexIds[2] - 1]
    };

    // Clip against each frustum plane
    currentVertices = clipAgainstPlane(currentVertices, Vec3(1, 0, 0), camera->left);   // Left plane
    currentVertices = clipAgainstPlane(currentVertices, Vec3(-1, 0, 0), camera->right); // Right plane
    currentVertices = clipAgainstPlane(currentVertices, Vec3(0, 1, 0), camera->bottom); // Bottom plane
    currentVertices = clipAgainstPlane(currentVertices, Vec3(0, -1, 0), camera->top);   // Top plane
    currentVertices = clipAgainstPlane(currentVertices, Vec3(0, 0, 1), camera->near);   // Near plane
    currentVertices = clipAgainstPlane(currentVertices, Vec3(0, 0, -1), camera->far);   // Far plane

    return currentVertices;
}
*/
/*
double edgeFunction(const Vec3 &a, const Vec3 &b, const Vec3 &c) {
    return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
}

void Scene::rasterizeTriangleEdges(int x0, int y0, int x1, int y1, int x2, int y2, vector<vector<double>> &depthBuffer, Camera *camera) {
    // Find the bounding box for the triangle
    int minX = std::max(0, std::min({x0, x1, x2}));
    int maxX = std::min(camera->horRes - 1, std::max({x0, x1, x2}));
    int minY = std::max(0, std::min({y0, y1, y2}));
    int maxY = std::min(camera->verRes - 1, std::max({y0, y1, y2}));

    // Loop through each pixel in the bounding box
    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            // Compute barycentric coordinates for the pixel (x, y)
            double w0 = edgeFunction(Vec3(x1, y1, 0), Vec3(x2, y2, 0), Vec3(x, y, 0));
            double w1 = edgeFunction(Vec3(x2, y2, 0), Vec3(x0, y0, 0), Vec3(x, y, 0));
            double w2 = edgeFunction(Vec3(x0, y0, 0), Vec3(x1, y1, 0), Vec3(x, y, 0));

            // Check if the pixel is inside the triangle
            if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                // Interpolate depth using the barycentric coordinates
                double area = w0 + w1 + w2;
                double z = (w0 * (x0 * y0) + w1 * (x1 * y1) + w2 * (x2 * y2)) / area;

                // Perform depth test (update depth buffer if closer)
                if (z < depthBuffer[x][y]) {
                    depthBuffer[x][y] = z;

                    // Interpolate the color (if needed, e.g., from vertex colors)
                    // Update color buffer (assuming you have a color buffer for the pixels)
                    // Set pixel color to the interpolated color
                    // Example: this->image[x][y] = color; // Apply color interpolation logic here
                }
            }
        }
    }
}


void Scene::rasterizeClippedTriangle(const vector<Vec3> &clippedVertices, vector<vector<double>> &depthBuffer, Camera *camera) {
    // The clipped triangle might have more than 3 vertices after clipping
    // Break it into a fan of triangles
    for (size_t i = 1; i + 1 < clippedVertices.size(); i++) {
        Vec3 v0 = clippedVertices[0];
        Vec3 v1 = clippedVertices[i];
        Vec3 v2 = clippedVertices[i + 1];

        // Convert normalized coordinates to screen space
        int x0 = (int)((v0.x - camera->left) / (camera->right - camera->left) * camera->horRes);
        int y0 = (int)((v0.y - camera->bottom) / (camera->top - camera->bottom) * camera->verRes);
        int x1 = (int)((v1.x - camera->left) / (camera->right - camera->left) * camera->horRes);
        int y1 = (int)((v1.y - camera->bottom) / (camera->top - camera->bottom) * camera->verRes);
        int x2 = (int)((v2.x - camera->left) / (camera->right - camera->left) * camera->horRes);
        int y2 = (int)((v2.y - camera->bottom) / (camera->top - camera->bottom) * camera->verRes);

        // Rasterize the triangle
        rasterizeTriangleEdges(x0, y0, x1, y1, x2, y2, depthBuffer, camera);
    }
}
*/

bool visible(float den, float num, float &t_e, float &t_l) {
	if (den > 0) { 
        float t = num / den;
        if (t > t_l) return false;
        if (t > t_e) t_e = t;
    } else if (den < 0) {
        float t = num / den;
        if (t < t_e) return false;
        if (t < t_l) t_l = t;
    } else if (num > 0) {
        return false;
    }
    return true;
}

std::vector<Vec3> Scene::clipLine(Vec3 &vertex1, Vec3 &vertex2, Camera *camera) {
    float x_min = camera->left;
    float x_max = camera->right;
    float y_min = camera->bottom;
    float y_max = camera->top;
    float z_min = camera->near;
    float z_max = camera->far;

    float t_e = 0.0f;
    float t_l = 1.0f;

    float dx = vertex2.x - vertex1.x;
    float dy = vertex2.y - vertex1.y;
    float dz = vertex2.z - vertex1.z;

	std::vector<Vec3> resultVertices (2, {-1,-1,-1,-1});
    // Check visibility in the x, y, z directions
    if (visible(dx, x_min - vertex1.x, t_e, t_l) &&
        visible(-dx, vertex1.x - x_max, t_e, t_l) &&
        visible(dy, y_min - vertex1.y, t_e, t_l) &&
        visible(-dy, vertex1.y - y_max, t_e, t_l) &&
        visible(dz, z_min - vertex1.z, t_e, t_l) &&
        visible(-dz, vertex1.z - z_max, t_e, t_l)) {

        // If the line is partially inside the frustum, update the vertices
        if (t_l < 1) {
            vertex2.x = vertex1.x + dx * t_l;
            vertex2.y = vertex1.y + dy * t_l;
            vertex2.z = vertex1.z + dz * t_l;
        }
        if (t_e > 0) {
            vertex1.x = vertex1.x + dx * t_e;
            vertex1.y = vertex1.y + dy * t_e;
            vertex1.z = vertex1.z + dz * t_e;
        }
        resultVertices[0] = vertex1;
		resultVertices[1] = vertex2;
    }
	resultVertices[0].colorId = vertex1.colorId;
	resultVertices[1].colorId = vertex2.colorId;
	// cout<<vertex1.colorId<<endl;
	return resultVertices;
}

// TODO this function should be changed !!
void Scene::drawLine(Vec3 v0, Vec3 v1, vector<vector<double>> &depthBuffer, Camera *camera) {
    // Convert floating-point coordinates to integers
    int x0 = round(v0.x), y0 = round(v0.y);
    int x1 = round(v1.x), y1 = round(v1.y);

    // Calculate differences
    int dx = x1 - x0;
    int dy = y1 - y0;

    // Determine step direction
    int sx = (dx >= 0) ? 1 : -1;
    int sy = (dy >= 0) ? 1 : -1;

    // Absolute values of differences
    dx = abs(dx);
    dy = abs(dy);

    // Decision parameter initialization
    int d; // Decision parameter
    int x = x0, y = y0;

    if (dx > dy) { // Horizontal or shallow slope
        d = 2 * dy - dx;
        while (x != x1) {
            // Rasterize current pixel
            if (x >= 0 && x < camera->horRes && y >= 0 && y < camera->verRes) {
                depthBuffer[x][y] = v0.z;

                // Testing color; update as needed
                Color c1;
                c1.r = 0;
                c1.g = 0;
                c1.b = 0;
                this->image[x][y] = c1;
            }
            if (d > 0) {
                y += sy;
                d -= 2 * dx;
            }
            d += 2 * dy;
            x += sx;
        }
    } else { // Vertical or steep slope
        d = 2 * dx - dy;
        while (y != y1) {
            // Rasterize current pixel
            if (x >= 0 && x < camera->horRes && y >= 0 && y < camera->verRes) {
                depthBuffer[x][y] = v0.z;

                // Testing color; update as needed
                Color c1;
                c1.r = 15;
                c1.g = 15;
                c1.b = 15;
                this->image[x][y] = c1;
            }
            if (d > 0) {
                x += sx;
                d -= 2 * dy;
            }
            d += 2 * dx;
            y += sy;
        }
    }
}


void Scene::rasterizeEdges(vector<Vec3> transformedVertices, const Triangle &triangle, vector<vector<double>> &depthBuffer, Camera *camera) {
    drawLine(transformedVertices[triangle.vertexIds[0]-1], transformedVertices[triangle.vertexIds[1]-1], depthBuffer, camera);
	drawLine(transformedVertices[triangle.vertexIds[1]-1], transformedVertices[triangle.vertexIds[2]-1], depthBuffer, camera);
	drawLine(transformedVertices[triangle.vertexIds[2]-1], transformedVertices[triangle.vertexIds[0]-1], depthBuffer, camera);
}


void Scene::forwardRenderingPipeline(Camera *camera, Scene *scene) {
	// TODO: Implement this function

	vector<vector<double>> depthBuffer(camera->horRes, vector<double>(camera->verRes, std::numeric_limits<double>::infinity()));
	vector<Vec3> transformedVertices (scene->vertices.size(), {0, 0, 0});

	for(Mesh *mesh : meshes) {
		// int currentMeshId = mesh->meshId;

		// cout << "gercek size "<< scene->vertices.size()<<  "     Processing rMesh ID: " << currentMeshId << endl;
        int i = 0;
		for(auto &vertices1 : scene->vertices) {

			modelingTransformation(*vertices1, vertices1->colorId, mesh->transformationTypes, mesh->transformationIds, scene, mesh);

			viewTransformations(*vertices1, camera);

			projectTransformations(*vertices1, camera);

			transformedVertices[i] = *vertices1;
			i++;
            // cout << transformedVertices[i-1] << endl;
		}

        for(auto &triangle : mesh->triangles){
            // isBackface(triangle, transformedVertices, camera);
			// cout << isBackface(triangle, transformedVertices, camera) << endl;

			// cout<<mesh->type<<endl;

			if (mesh->type == 0) { // (mesh->type == "wireframe")
				// cout << "basla   " << transformedVertices[triangle.vertexIds[0] - 1] <<transformedVertices[triangle.vertexIds[1] - 1] <<endl;
				std::vector<Vec3> clippedLine = clipLine(transformedVertices[triangle.vertexIds[0] - 1], transformedVertices[triangle.vertexIds[1] - 1], camera);
				transformedVertices[triangle.vertexIds[0] - 1] = clippedLine[0];
				transformedVertices[triangle.vertexIds[1] - 1] = clippedLine[1];

				clippedLine = clipLine(transformedVertices[triangle.vertexIds[1] - 1], transformedVertices[triangle.vertexIds[2] - 1], camera);
				transformedVertices[triangle.vertexIds[1] - 1] = clippedLine[0];
				transformedVertices[triangle.vertexIds[2] - 1] = clippedLine[1];

				clippedLine = clipLine(transformedVertices[triangle.vertexIds[2] - 1], transformedVertices[triangle.vertexIds[0] - 1], camera);
				transformedVertices[triangle.vertexIds[2] - 1] = clippedLine[0];
				transformedVertices[triangle.vertexIds[0] - 1] = clippedLine[1];
				// clipLine(transformedVertices[triangle.vertexIds[1] - 1], transformedVertices[triangle.vertexIds[2] - 1], camera);
				// clipLine(transformedVertices[triangle.vertexIds[2] - 1], transformedVertices[triangle.vertexIds[0] - 1], camera);
				// cout <<"2.      " << transformedVertices[triangle.vertexIds[0] - 1] << transformedVertices[triangle.vertexIds[1] - 1]<<"\n \n"<<endl;
			}

			if (scene->cullingEnabled && isBackface(triangle, transformedVertices, camera)) {
                continue;
            }

			if (mesh->type == 0) {
                // rasterizeEdges(transformedVertices, triangle, depthBuffer, camera);
            } 

			// if the object is wireframe mesh: clip vertices first 

            // vector<Vec3> clippedVertices = clipTriangle(triangle, transformedVertices, camera);
            // if (clippedVertices.empty()) continue;

            // rasterizeClippedTriangle(clippedVertices, depthBuffer, camera);
		}

		
	}

}

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
    Vec4 triangleHc = { vertices1.x, vertices1.y, vertices1.z, 1.0f };
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
	
	// Vec3 result;
	// result.x = triangleHc.x;
	// result.y = triangleHc.y;
	// result.z = triangleHc.z;
	// result.colorId = vertices1.colorId;
	// return result;
	
	// cout << vertices1.colorId << "   " <<vertices1.x << "   " << vertices1.y << "   "  << vertices1.z << endl;

}

// kamerada gaze ve up vectorleri bazen dik gelmiyor ornegin: emptyboxclipped 1. olanindan.. !!
void viewTransformations(Vec3 &vertices1, Camera *camera) {
    Matrix4 viewMatrix;
	// cout << viewMatrix <<endl;
	// cout<< vertices1.x << "  " << vertices1.y << "  " << vertices1.z << endl;
	Vec4 triangleHc = { vertices1.x, vertices1.y, vertices1.z, 1.0f };

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
	Vec4 triangleHc = { vertices1.x, vertices1.y, vertices1.z, 1.0f };

	// cout << "projection type: " << camera->projectionType <<"   "<<  
	// 		(camera->right - camera->left)<< "  " <<
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
		projectionMatrix.values[3][2] = -1.0f;
		projectionMatrix.values[3][3] = 0.0f;
	} else {
		printf("something wrong in project transformation");
	}

	// cout << projectionMatrix << " \n \n " << endl;
	triangleHc = multiplyMatrixWithVec4(projectionMatrix, triangleHc);
	// cout << triangleHc << "\n\n"<<endl;
	if (triangleHc.t != 0){
		vertices1.x = triangleHc.x/triangleHc.t;
		vertices1.y = triangleHc.y/triangleHc.t;
		vertices1.z = triangleHc.z/triangleHc.t;
	} else {
		cout << "something wrong in the project transformation" << endl;
	}
    // cout << vertices1.x << "   " << vertices1.y << "   "  << vertices1.z <<"\n\n" << endl;
}

Matrix4 viewportTransformation(Camera *camera) {
	Matrix4 vpMatrix;
	// Vec4 triangleHc = { vertices1.x, vertices1.y, vertices1.z, 1.0f };

	float n_x = camera->horRes;
	float n_y = camera->verRes;

	vpMatrix.values[0][0] = (n_x/2.0f);
	vpMatrix.values[0][1] = 0;
	vpMatrix.values[0][2] = 0;
	vpMatrix.values[0][3] = (n_x - 1)/2.0f + camera->left; // it can be wrong as camera->left is x_min

	vpMatrix.values[1][0] = 0;
	vpMatrix.values[1][1] = (n_y/2.0f);
	vpMatrix.values[1][2] = 0;
	vpMatrix.values[1][3] = ((n_y-1)/2.0f) + camera->bottom; // it can be wrong that camera->bottom is y_min

	vpMatrix.values[2][0] = 0;
	vpMatrix.values[2][1] = 0;
	vpMatrix.values[2][2] = 0.5f;
	vpMatrix.values[2][3] = 0.5f;

	vpMatrix.values[3][0] = 0;
	vpMatrix.values[3][1] = 0;
	vpMatrix.values[3][2] = 0;
	vpMatrix.values[3][3] = 0;

	return vpMatrix;
	// cout<<vpMatrix<<"\n\n"<<endl;
	// triangleHc = multiplyMatrixWithVec4(vpMatrix, triangleHc);
	// cout<<triangleHc<<"\n\n"<<endl;
	// vertices1.x = triangleHc.x;
	// vertices1.y = triangleHc.y;
	// vertices1.z = triangleHc.z;

}

bool isBackface(Triangle &triangle, vector<Vec3> &transformedVertices, Camera *camera) {
    Vec3 v0 = transformedVertices[triangle.vertexIds[0] - 1];
    Vec3 v1 = transformedVertices[triangle.vertexIds[1] - 1];
    Vec3 v2 = transformedVertices[triangle.vertexIds[2] - 1];

    Vec3 normal = crossProductVec3(subtractVec3(v0, v1), subtractVec3(v0, v2)); // Ensure correct order
    Vec3 viewDir = subtractVec3(v0, camera->position);

    const float bEPSILON = 1e-6;
    float dot = dotProductVec3(normalizeVec3(normal), normalizeVec3(viewDir));

    return dot < -bEPSILON; // Use epsilon for robustness
}



bool visible(float den, float num, float &t_e, float &t_l) {
	if (den > 0) { // PE
        float t = -num / den;
        if (t > t_l) return false;
        if (t > t_e) t_e = t;
    } else if (den < 0) { //PL
        float t = -num / den;
        if (t < t_e) return false;
        if (t < t_l) t_l = t;
    } else if (num > 0) { // line parallel to edge
        return false;
    }
    return true;
}

std::vector<Vec3> Scene::clipLine(Vec3 &vertex1, Vec3 &vertex2, Camera *camera) {
    float x_min = camera->left;
    float x_max = camera->right;
    float y_min = camera->bottom;
    float y_max = camera->top;
    float z_min = -1;
    float z_max = 1;

    float t_e = 0.0f;
    float t_l = 1.0f;

    float dx = vertex2.x - vertex1.x;
    float dy = vertex2.y - vertex1.y;
    float dz = vertex2.z - vertex1.z;

	std::vector<Vec3> resultVertices (2, {0,0,0}); // clippinte hata burada olabilir..
	resultVertices[0] = vertex1;
	resultVertices[1] = vertex2;
	
    // Check visibility in the x, y, z directions
    if (visible(dx, x_min - vertex1.x, t_e, t_l) &&
        visible(-dx, vertex1.x - x_max, t_e, t_l) &&
        visible(dy, y_min - vertex1.y, t_e, t_l) &&
        visible(-dy, vertex1.y - y_max, t_e, t_l) && visible(dz, z_min-vertex1.z, t_e, t_l) && visible(-dz, vertex1.z-z_max, t_e, t_l)) { // burda z yazan yerleri sildiydim.. // && visible(dz, z_min-vertex1.z, t_e, t_l) && visible(-dz, vertex1.z-z_max, t_e, t_l)

        // If the line is partially inside the frustum, update the vertices
        if (t_l < 1) {
			// cout<<t_l<< "  " << "t_l  " <<  vertex1.z <<"" <<endl;
            vertex2.x = vertex1.x + dx * t_l;
            vertex2.y = vertex1.y + dy * t_l;
            vertex2.z = vertex1.z + dz * t_l;
        }
        if (t_e > 0) {
			// cout<<t_e<< "  " << "t_e  " <<  vertex1.x <<"\n" <<endl;
            vertex1.x = vertex1.x + dx * t_e;
            vertex1.y = vertex1.y + dy * t_e;
            vertex1.z = vertex1.z + dz * t_e;
        }
        resultVertices[0] = vertex1;
		resultVertices[1] = vertex2;
    }
	resultVertices[0].colorId = vertex1.colorId;
	resultVertices[1].colorId = vertex2.colorId;
	// cout<<resultVertices[0]<< "  "<< resultVertices[1]<<"\n\n"<<endl;
	return resultVertices;
}


/*
// TODO this function should be changed !!
void Scene::drawLine(Vec3 v0, Vec3 v1, vector<vector<double>> &depthBuffer, Camera *camera) {
	int y = abs(round(v0.y));
	int d = abs(v0.y - v1.y) + 0.5* abs(v1.x - v0.x);

	// float c = color[0];
	// float dc = (color1 - color0) / (x1 - x0)

	Color c;
	c.r = 0;
	c.g = 0;
	c.b = 0;

	for (int x = abs(v0.x); x < abs(v1.x); x++) {
		// cout << x << "    " << y<<endl;
		this->image[x][y] = c;
		if (d<0) { // choose NE
			y += 1;
			d += abs(v0.y - v1.y) + abs(v1.x - v0.x);
		} else {
			d += abs(v0.y - v1.y);
		}
		// c += dc;
	}

}
*/

// TODO this function should be changed !!
void Scene::drawLine(Vec3 v0, Vec3 v1, vector<vector<double>> &depthBuffer, Camera *camera) {

    int x0 = round(v0.x), y0 = round(v0.y);
    int x1 = round(v1.x), y1 = round(v1.y);

    int dx = x1 - x0;
    int dy = y1 - y0;

    int sx = (dx >= 0) ? 1 : -1;
    int sy = (dy >= 0) ? 1 : -1;

    dx = abs(dx);
    dy = abs(dy);

    Color startColor = *colorsOfVertices[v0.colorId - 1];
    Color endColor = *colorsOfVertices[v1.colorId - 1];

    int steps = (dx > dy) ? dx : dy;

    double dr = (endColor.r - startColor.r) / static_cast<double>(steps);
    double dg = (endColor.g - startColor.g) / static_cast<double>(steps);
    double db = (endColor.b - startColor.b) / static_cast<double>(steps);

    double r = startColor.r, g = startColor.g, b = startColor.b;

    int d;
    int x = x0, y = y0;

    if (dx > dy) { 
        d = 2 * dy - dx;
        for (int i = 0; i <= steps; ++i) {

            if (x >= 0 && x < camera->horRes && y >= 0 && y < camera->verRes) {
                depthBuffer[x][y] = v0.z;

                Color c1;
                c1.r = r;
                c1.g = g;
                c1.b = b;
                this->image[x][y] = c1;
            }

            r += dr;
            g += dg;
            b += db;

            if (d > 0) {
                y += sy;
                d -= 2 * dx;
            }
            d += 2 * dy;
            x += sx;
        }
    } else { 
        d = 2 * dx - dy;
        for (int i = 0; i <= steps; ++i) {

            if (x >= 0 && x < camera->horRes && y >= 0 && y < camera->verRes) {
                depthBuffer[x][y] = v0.z;

                Color c1;
                c1.r = r;
                c1.g = g;
                c1.b = b;
                this->image[x][y] = c1;
            }

            r += dr;
            g += dg;
            b += db;

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

double f_01(int x, int y, int x0, int y0, int x1, int y1) {
	return x*(y0-y1) + y*(x1-x0) + x0*y1 - y0*x1;
}

double f_12(int x, int y, int x1, int y1, int x2, int y2) {
	return x*(y1-y2) + y*(x2-x1) + x1*y2 - y1*x2;
}

double f_20(int x, int y, int x2, int y2, int x0, int y0) {
	return x*(y2-y0) + y*(x0-x2) + x2*y0 - y2*x0;
}

void Scene::rasterizeFilledTriangle(Triangle &triangle, std::vector<Vec3> &transformedVertices, Camera *camera, std::vector<std::vector<double>> &depthBuffer) {
	Vec3 v0 = transformedVertices[triangle.vertexIds[0] - 1];
	Vec3 v1 = transformedVertices[triangle.vertexIds[1] - 1];
	Vec3 v2 = transformedVertices[triangle.vertexIds[2] - 1];

	Color c;
	// cout<<v0.colorId<<endl;
	Color c0 = *colorsOfVertices[v0.colorId - 1];
    Color c1 = *colorsOfVertices[v1.colorId - 1];
    Color c2 = *colorsOfVertices[v2.colorId - 1];
	
	int bbox_min_x = std::min(v0.x, std::min(v1.x, v2.x));
	int bbox_min_y = std::min(v0.y, std::min(v1.y, v2.y));

	int bbox_max_x = std::max(v0.x, std::max(v1.x, v2.x));
	int bbox_max_y = std::max(v0.y, std::max(v1.y, v2.y));

	// double f_12_1 = f_12(v0.x,v0.y, v1.x, v1.y, v2.x, v2.y);
	// double f_20_2 = f_20(v1.x,v1.y, v2.x, v2.y, v0.x, v0.y);
	// double f_01_3 = f_01(v2.x,v2.y, v0.x, v0.y, v1.x, v1.y);
	// cout<<f_12_1<<"  "<<f_20_2<<"  "<<f_01_3<<endl;
	
	double triangleArea = f_12(v0.x,v0.y, v1.x, v1.y, v2.x, v2.y);
	// bounding box in here can be made more efficient as the polygon rasterization algorithm..
	for (int y=bbox_min_y; y<=bbox_max_y; y++) {
		for (int x=bbox_min_x; x<=bbox_max_x; x++) {
			double alpha = f_12(x,y, v1.x, v1.y, v2.x, v2.y) / triangleArea; 
			double beta = f_20(x,y, v2.x, v2.y, v0.x, v0.y) / triangleArea; 
			double gamma = f_01(x,y, v0.x, v0.y, v1.x, v1.y) / triangleArea;
			if (alpha>=0 && beta >= 0 && gamma>=0 && x >= 0 && x < camera->horRes && y >= 0 && y < camera->verRes) { // maybe <= camera->horRes would be better 
				// cout<<alpha * v0.z<<"  "<<beta<<"  "<<v2.z<<"\n"<<endl;
				double depth = alpha * v0.z + beta * v1.z + gamma * v2.z;
				if (depth < depthBuffer[x][y]) {
					// cout<<depthBuffer[x][y]<<"  "<<depth<<endl;
					depthBuffer[x][y] = depth;

					c.r = alpha*c0.r + beta * c1.r + gamma*c2.r;
					c.g = alpha*c0.g + beta * c1.g + gamma*c2.g;
					c.b = alpha*c0.b + beta * c1.b + gamma*c2.b;

					this->image[x][y].r = round(c.r);
					this->image[x][y].g = round(c.g);
					this->image[x][y].b = round(c.b);
				}
			}
		}
	}

	// tri_min_z = std::min(v0.z, std::min(v2.z, v3.z));
	// cout<<bbox_min_x << "  "<<bbox_min_y<<endl;
	

}


void Scene::forwardRenderingPipeline(Camera *camera, Scene *scene) {
	// TODO: Implement this function

	vector<vector<double>> depthBuffer(camera->horRes, vector<double>(camera->verRes, std::numeric_limits<double>::infinity()));
	vector<Vec3> transformedVertices (scene->vertices.size(), {0, 0, 0});
	Matrix4 vpMatrix = viewportTransformation(camera);

	vector<Vec3> tempVertices (scene->vertices.size(), {0, 0, 0});

	for(Mesh *mesh : meshes) {
		// int currentMeshId = mesh->meshId;

		// cout << "gercek size "<< scene->vertices.size()<<  "     Processing rMesh ID: " << currentMeshId << endl;

		int j=0;
		for(auto &vertices1 : scene->vertices) {
			tempVertices[j] = *vertices1;
			j++;
		}

        int i = 0;
		for(auto &vertices1 : scene->vertices) {
			// cout<< colorsOfVertices[i-1]->r << "\n"<<endl;
			modelingTransformation(*vertices1, vertices1->colorId, mesh->transformationTypes, mesh->transformationIds, scene, mesh);

			viewTransformations(*vertices1, camera);

			projectTransformations(*vertices1, camera);

			// viewportTransformation(*vertices1, camera);

			// alttaki 5 satir clipped islemi oncesinde yapiliyosa burasi kapatilip asagidaki acilmali
			// Vec4 vpMultiply = {vertices1->x, vertices1->y, vertices1->z, 1.0f};
			// Vec4 multipliedMatrix = multiplyMatrixWithVec4(vpMatrix, vpMultiply);
			// vertices1->x = multipliedMatrix.x;
			// vertices1->y = multipliedMatrix.y;
			// vertices1->z = multipliedMatrix.z;


			transformedVertices[i] = *vertices1;
			*vertices1 = tempVertices[i];
			i++;
			// cout<<*vertices1<<"  "<<endl;
            // cout << transformedVertices[i-1] << endl;
		}

		// clipping burda yapilabilir.. fikrim degisti -> bence clipping zaten draw line icinde yapiliyor tekrar yapmaya gerek yok..
		
        for(auto &triangle : mesh->triangles){
            // isBackface(triangle, transformedVertices, camera);
			// cout << isBackface(triangle, transformedVertices, camera) << endl;

			// cout<<mesh->type<<endl;
			
			if (mesh->type == 0) { // (mesh->type == "wireframe")
				cout << "basla " << transformedVertices[triangle.vertexIds[0] - 1] 
					<< transformedVertices[triangle.vertexIds[1] - 1] 
					<<transformedVertices[triangle.vertexIds[2] - 1] <<endl;

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
				
				cout <<"2.    " << transformedVertices[triangle.vertexIds[0] - 1] 
					<< transformedVertices[triangle.vertexIds[1] - 1] 
					<<transformedVertices[triangle.vertexIds[2] - 1] 
					<<"\n \n"<<endl;
			} 

		}
		
		// burasi viewport transformations clipping isleminden sonra yapiliyosa acilabilir..

		for (int k = 0; k< transformedVertices.size();k++){
			Matrix4 vpMatrix = viewportTransformation(camera);
			Vec4 vpMultiply = {transformedVertices[k].x, transformedVertices[k].y, transformedVertices[k].z, 1.0f};
			Vec4 multipliedMatrix = multiplyMatrixWithVec4(vpMatrix, vpMultiply);
			transformedVertices[k].x = multipliedMatrix.x;
			transformedVertices[k].y = multipliedMatrix.y;
			transformedVertices[k].z = multipliedMatrix.z;
			// cout<<transformedVertices[k]<<endl;
		}
		

		for(auto &triangle : mesh->triangles) {

			if (scene->cullingEnabled && isBackface(triangle, transformedVertices, camera) && scene->meshes.size() == 1) {
				continue;
			}
			if (scene->cullingEnabled && !isBackface(triangle, transformedVertices, camera) && scene->meshes.size() > 1) {
				continue;
			}
			if (mesh->type == 0) {
                rasterizeEdges(transformedVertices, triangle, depthBuffer, camera);
            } else {
				// cout<<"before seg fault:  " << transformedVertices[triangle.vertexIds[0]-1]<<"  "<<triangle.vertexIds[1]<<"  "<<triangle.vertexIds[2]<<endl;
				// should be clipped before enter this function 
				rasterizeFilledTriangle(triangle, transformedVertices, camera, depthBuffer);
			}
		}
		
	}

}

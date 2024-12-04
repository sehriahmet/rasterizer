#include <iostream>
#include "Scene.h"

using namespace std;

Scene *scene;

#include <chrono> // timer library delete when submitting 1 

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Please run the rasterizer as:" << endl
             << "\t./rasterizer <input_file_name>" << endl;
        return 1;
    }
    else
    {
        const char *xmlPath = argv[1];

        scene = new Scene(xmlPath);

        for (int i = 0; i < scene->cameras.size(); i++)
        {
            // initialize image with basic values
            scene->initializeImage(scene->cameras[i]);

            auto a = std::chrono::high_resolution_clock::now(); // delete timing 2

            // do forward rendering pipeline operations
            scene->forwardRenderingPipeline(scene->cameras[i], scene);

            /* delete timing when submitting start 3 */
            auto b = std::chrono::high_resolution_clock::now();
            int durationMs = duration_cast<std::chrono::microseconds>(b - a).count();
            int durationS = duration_cast<std::chrono::seconds>(b - a).count();
            if (durationS/60 > 0) cout << scene->cameras[i]->outputFilename << " took " << durationS/60 << " minutes " << (durationS%60) << " seconds" << endl;
            else if (durationS > 0) cout << scene->cameras[i]->outputFilename << " took " << durationS%60 << " seconds" << endl;
            else cout << scene->cameras[i]->outputFilename << " took " << "0" << " seconds " << duration_cast<std::chrono::microseconds>(b - a).count() << " microseconds" <<  endl;
            /* delete timing when submitting 3*/

            // generate PPM file
            scene->writeImageToPPMFile(scene->cameras[i]);

            // Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
            // Notice that os_type is not given as 1 (Ubuntu) or 2 (Windows), below call doesn't do conversion.
            // Change os_type to 1 or 2, after being sure that you have ImageMagick installed.
            scene->convertPPMToPNG(scene->cameras[i]->outputFilename, 0);
            delete scene;
            scene = new Scene(xmlPath);
        }
        cout<<"\n"<<endl; // delete when submitting 4
        return 0;
    }
}
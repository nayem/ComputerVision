// B657 assignment 3 skeleton code, D. Crandall
//
// Compile with: "make"
//
// This skeleton code implements nearest-neighbor classification
// using just matching on raw pixel values, but on subsampled "tiny images" of
// e.g. 20x20 pixels.
//
// It defines an abstract Classifier class, so that all you have to do
// :) to write a new algorithm is to derive a new class from
// Classifier containing your algorithm-specific code
// (i.e. load_model(), train(), and classify() methods) -- see
// NearestNeighbor.h for a prototype.  So in theory, you really
// shouldn't have to modify the code below or the code in Classifier.h
// at all, besides adding an #include and updating the "if" statement
// that checks "algo" below.
//
// See assignment handout for command line and project specifications.
//
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>
#include <sys/types.h>
#include <dirent.h>
#include <map>
#include <numeric>
#include <fstream>
#include <sstream>
#include <time.h>
#include <unistd.h>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

// Dataset data structure, set up so that e.g. dataset["bagel"][3] is
// filename of 4th bagel image in the dataset
typedef map<string, vector<string> > Dataset;

#include <Classifier.h>
#include <NearestNeighbor.h>
/* Edited by Nayem*/
#include <SVM.h>

// Figure out a list of files in a given directory.
//
vector<string> files_in_directory(const string &directory, bool prepend_directory = false)
{
    vector<string> file_list;
    DIR *dir = opendir(directory.c_str());
    cout<<directory.c_str()<<endl;
    if(!dir)
        throw std::string("Can't find directory " + directory);
    
    struct dirent *dirent;
    while ((dirent = readdir(dir)))
        if(dirent->d_name[0] != '.')
            file_list.push_back((prepend_directory?(directory+"/"):"")+dirent->d_name);
    
    closedir(dir);
    return file_list;
}

int main(int argc, char **argv)
{
    try {
        if(argc < 3)
            throw string("Insufficent number of arguments");
        
        string mode = argv[1];
        string algo = argv[2];

	/* Edited by Nayem*/
        // command-> ./a3 train svm fooddata2 COLOR
        // command-> ./a3 test svm fooddata2 COLOR
        string FILE_DIRECTORY = argv[3];
        string COLOR_MODE = argv[4]; // COLOR or GRAY
        
        // Scan through the "train" or "test" directory (depending on the
        //  mode) and builds a data structure of the image filenames for each class.
        Dataset filenames;
        
        /* Edited by Nayem*/
        vector<string> class_list = files_in_directory(FILE_DIRECTORY + "/" +mode);
        //vector<string> class_list = files_in_directory(mode);
        
        for(vector<string>::const_iterator c = class_list.begin(); c != class_list.end(); ++c)
            /* Edited by Nayem*/
            filenames[*c] = files_in_directory(FILE_DIRECTORY + "/" +mode + "/" + *c, true);
            //filenames[*c] = files_in_directory(mode + "/" + *c, true);
        
        // set up the classifier based on the requested algo
        Classifier* classifier = NULL;
	classifier = new SVM(class_list,algo);
        
        // now train or test!
        if(mode == "train"){
            /* Edited by Nayem, Part-1.1*/
            classifier->train(filenames, COLOR_MODE);
            int ret;
            
            char cwd[1024];
            chdir("/path/to/change/directory/to");
            getcwd(cwd, sizeof(cwd));
            //printf("Train-> Current working dir: %s\n", cwd);
            
            string directory = "svm_multiclass";
            ret = chdir (directory.c_str());
            //cout<<"chdir():"<<ret<<endl;
            
            chdir("/path/to/change/directory/to");
            getcwd(cwd, sizeof(cwd));
            //printf("Train-> Current working dir: %s\n", cwd);
            
            ret = system("./svm_multiclass_learn -c 5000 train.dat model");
            //cout<<"Train-> run():"<<ret<<endl;

        }
        else if(mode == "test"){
            /* Edited by Nayem, Part-1.3*/
            ((SVM *)classifier)->test(filenames, COLOR_MODE);
            int ret;
            
            char cwd[1024];
            chdir("/path/to/change/directory/to");
            getcwd(cwd, sizeof(cwd));
            //printf("Test-> Current working dir: %s\n", cwd);
            
            string directory = "svm_multiclass";
            ret = chdir (directory.c_str());
            //cout<<"chdir():"<<ret<<endl;
            
            chdir("/path/to/change/directory/to");
            getcwd(cwd, sizeof(cwd));
            //printf("Test-> Current working dir: %s\n", cwd);
            
            ret = system("./svm_multiclass_classify test.dat model predictions");
            //cout<<"Test-> run():"<<ret<<endl;
            
        }
        else
            throw std::string("unknown mode!");
    }
    catch(const string &err)
    {
        cerr << err << endl;
    }
}
#include<set>
#include "extract_haar_feature.h"
class SVM : public Classifier
{
public:
    string featureType;
    SVM(const vector<string> &_class_list,string featureType) : Classifier(_class_list)
    {
        this->featureType = featureType;
    }
    
    // SVM training. All this does is read in all the images, resize
    // them to a common size, and dump them as vectors to a file
    // PART 1.2
    virtual void train(const Dataset &filenames, const string COLOR_MODE="COLOR")
    {
        //create train.dat file
        ofstream train_output_File;
        train_output_File.open("svm_multiclass/train.dat");
        int class_counter = 1;
        
        CImg<double> class_vectors;
        if (COLOR_MODE == "GRAY")
        {
            CImg<double> class_vectors_temp(size*size*1, filenames.size(), 1);
            class_vectors = class_vectors_temp;
        }
        else
        {
            CImg<double> class_vectors_temp(size*size*3, filenames.size(), 1); // COLOR
            class_vectors = class_vectors_temp;
        }
        
        if(featureType == "baseline")
        {
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                for(int i=0; i<c_iter->second.size(); i++)
                {
                    CImg<double> temp = extract_features(c_iter->second[i].c_str(), COLOR_MODE);
                    //class_vectors = class_vectors.draw_image(0, i, 0, 0, temp);
                    writeTrainFile(train_output_File, class_counter, temp);
                }
            }
            //class_vectors.save_png(("nn_model." + c_iter->first + ".png").c_str());
        }
        else if(featureType == "eigen")
        {
            vector< CImg<double> > v;
            cout<<"In the train function for eigen features"<<endl;
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                for(int i=0; i<c_iter->second.size(); i++)
                {
                    CImg<double> temp = extract_features(c_iter->second[i].c_str(), COLOR_MODE);
                    v.push_back(temp);
                }
            }
            //cout<<"Loaded the vector"<<endl;
            CImg<double> matrix(v[0].size(),v.size());
            for(int i=0;i<v.size();i++)
            {
                for(int j=0;j<v[i].size();j++)
                {
                    matrix(j,i) = v[i][j];
                }
            }
            CImg<double> coVariance = matrix.get_transpose()*matrix;
            //cout<<"Created the covariance matrix"<<endl;
            //cout<<"Size of the covariance matrix: "<<coVariance.width()<<" "<<coVariance.height()<<endl;
            CImg<double> U,D,V;
            coVariance.SVD(U,D,V,true);
            //cout<<U.width()<<" "<<U.height()<<endl;
            //cout<<D.width()<<" "<<D.height()<<endl;
            //cout<<V.width()<<" "<<V.height()<<endl;
            
            ofstream pca_File;
            pca_File.open("pca.csv");
            int reducedDimension;
            for (int i=0; i<D.height(); ++i)
            {
                if(D(i,1)>=300.0)
                    reducedDimension = i+1;
            }
            pca_File<<reducedDimension<<"\n";
            for(int i=0;i<U.height();i++)
            {
                for(int j=0;j<reducedDimension;j++)
                {
                    if(j<(reducedDimension-1))
                        pca_File<<U(j,i)<<",";
                    else
                        pca_File<<U(j,i)<<"\n";
                }
            }
            pca_File.close();
            CImg<double> reducedDataMatrix = matrix*U.columns(0,reducedDimension-1);
            //cout<<reducedDataMatrix.width()<<" "<<reducedDataMatrix.height()<<endl;
            //find eigen decomposition
            //from there find the newly reduced vectors
            //then pass it on to the SVM classifiers
            class_counter = 1;
            int j=0;
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                //write every file in train_output_file
                for(int i=0; i<c_iter->second.size(); i++)
                {
                    writeTrainFile(train_output_File, class_counter, reducedDataMatrix.row(j));
                    j++;
                }
            }
        }
        else if(featureType == "haar")
        {
            vector <HaarFilter> haarFilters = createHaarFilters();
            writeHaarFiltersIntoFile(haarFilters);
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                cout << "class " << class_counter << " started"  << endl;
                for(int i=0; i<c_iter->second.size(); i++)
                {
                    CImg<double> input_image(c_iter->second[i].c_str());
                    CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
                    gray = gray.resize(25,25,1,1);
                    CImg<double> int_image(gray.width(), gray.height(), 1, 1);
                    integral_image(gray, int_image);
                    vector<double> feature = haar_feature_all_filters(int_image, haarFilters, 2, 2);
                    CImg<double> temp(feature.size());
                    if(class_counter == 1&& i == 0)
                        cout << "haar feature Size = " << feature.size() << endl;
                    for(int p=0;p<feature.size();p++)
                        temp(p) = feature[p];
                    writeTrainFile(train_output_File, class_counter, temp);
                }
            }
        }
        else if(featureType == "bow")
        {
            vector<vector<float> > v;
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                for(int i=0; i<c_iter->second.size(); i++)
                {
                    CImg<double> input_image(c_iter->second[i].c_str());
                    CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
                    vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
                    for(int j=0;j<descriptors.size();j++)
                    {
                        v.push_back(descriptors[j].descriptor);
                    }
                }
                cout<<"Added sift descriptors: "<<v.size()<<endl;
            }
            cout<<"Size of sift descriptor vector: "<<v.size()<<endl;
            vector<vector<float> > centroids;
            int K = 10;
            //runnig K-Means on V with K = 20
            runKMeans(v,K,centroids);
            cout<<"We have divided the whole sift descriptor vector into "<<centroids.size()<<" clusters"<<endl;
            
            ofstream centroid_File;
            centroid_File.open("centroids.txt");
            for (int i=0; i<K; ++i)
            {
                for(int j=0;j<centroids[i].size();j++)
                    if (j < centroids[i].size()-1)
                        centroid_File<<centroids[i][j]<<" ";
                    else
                        centroid_File<<centroids[i][j]<<"\n";
            }
            centroid_File.close();
            class_counter = 1;
            double* histogram = new double[K];
            //finding histogram for each image and writing it to training file for SVM
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                for(int i=0; i<c_iter->second.size(); i++)
                {
                    CImg<double> input_image(c_iter->second[i].c_str());
                    CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
                    vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
                    memset(histogram,0.0,sizeof(double)*K);
                    for(int j=0;j<descriptors.size();j++)
                    {
                        int index = getNearestCentroid(descriptors[j].descriptor,centroids);
                        histogram[index] = histogram[index]+1;
                    }
                    CImg<double> temp(K);
                    for(int p=0;p<K;p++)
                        temp(p) = histogram[p];
                    writeTrainFile(train_output_File, class_counter, temp);
                }
                cout<<"Done with writing the train output file for class: "<<class_counter<<endl;
            }
            delete histogram;
        }
        else if(featureType == "deep")
        {
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                for(int i=0; i< c_iter->second.size(); i++)
                {
                    char cwd[1024];
                    chdir("/path/to/change/directory/to");
                    getcwd(cwd, sizeof(cwd));
                    //printf("Train-> Current working dir: %s\n", cwd);
                    
                    string directory = "overfeat/bin/linux_64";
                    chdir (directory.c_str());
                    //cout<<"chdir():"<<ret<<endl;
                    
                    chdir("/path/to/change/directory/to");
                    getcwd(cwd, sizeof(cwd));
                    //printf("Train-> Current working dir: %s\n", cwd);
                    
                    // if image is smaller
                    char image_name[300] = "";
                    strcat(image_name, "../../../");
                    strcat(image_name, c_iter->second[i].c_str());
                    CImg<double> temp(image_name);
                    if (temp.width() < 231 && temp.height() < 231 ) {
                        cout<<"Resize W+H: "<<image_name<<endl;
                        temp.resize(DEEP_SUBSAMPLE_SIZE,DEEP_SUBSAMPLE_SIZE,1,3);
                        strcpy (image_name,"../../../resized_image.jpg");
                        // strcpy (image_name,"../../../fooddata2/train/");
                        // strcat (image_name,c_iter->first.c_str());
                        // strcat (image_name,"/resized_image.jpg");
                        temp.save(image_name);
                    }
                    else if (temp.width() < 231 ) {
                        cout<<"Resize W: "<<image_name<<endl;
                        temp.resize(DEEP_SUBSAMPLE_SIZE,DEEP_SUBSAMPLE_SIZE,1,3);
                        strcpy (image_name,"../../../resized_image.jpg");
                        // strcpy (image_name,"../../../fooddata2/train/");
                        // strcat (image_name,c_iter->first.c_str());
                        // strcat (image_name,"/resized_image.jpg");
                        temp.save(image_name);
                    }
                    else if (temp.height() < 231 ){
                        cout<<"Resize H: "<<image_name<<endl;
                        temp.resize(DEEP_SUBSAMPLE_SIZE,DEEP_SUBSAMPLE_SIZE,1,3);
                        strcpy (image_name,"../../../resized_image.jpg");
                        // strcpy (image_name,"../../../fooddata2/train/");
                        // strcat (image_name,c_iter->first.c_str());
                        // strcat (image_name,"/resized_image.jpg");
                        temp.save(image_name);
                    }
                    else{
                        cout<<"No Resize : "<<image_name<<endl;
                        temp.resize(DEEP_SUBSAMPLE_SIZE,DEEP_SUBSAMPLE_SIZE,1,3);
                        strcpy (image_name,"../../../resized_image.jpg");
                        // strcpy (image_name,"../../../fooddata2/train/");
                        // strcat (image_name,c_iter->first.c_str());
                        // strcat (image_name,"/resized_image.jpg");
                        temp.save(image_name);
                    }
                    
                    char command[300] = "";
                    strcat(command, "./overfeat -f ");
                    strcat(command, image_name);
                    strcat(command, " > ../../../feature.txt");
                    cout<<command<<endl;
                    system(command); // command overfeat
                    strcpy (command,"");
                    strcpy (image_name,"");
                    
                    //read feature file
                    ifstream feature_File;
                    string line, word;
                    feature_File.open ("../../../feature.txt");
                    
                    if (feature_File.is_open())
                    {
                        getline (feature_File,line);
                        CImg<double> temp(DEEP_FEATURE_SIZE);
                        for (int k=0; feature_File >> word; k++ )
                            temp[k] = (double)atof(word.c_str());
                        
                        class_vectors = class_vectors.draw_image(0, i, 0, 0, temp);
                        writeTrainFile(train_output_File, class_counter, temp);
                        cout<<"Size:"<<temp.size()<<endl;
                        feature_File.close();
                    }
                    //break;
                    
                }
                //break;
                
            }
            chdir("../../..");
        }
        train_output_File.close();
    }
    double getEuclideanDistance(vector<float> v1,vector<float> v2)
    {
        double distances = 0;
        for(int i=0;i<v1.size();i++)
        {
            distances += pow(v1[i]-v2[i],2);
        }
        return sqrt(distances);
    }
    int getNearestCentroid(vector<float> data,vector<vector<float> > centroids)
    {
        double* distances = new double[centroids.size()];
        for(int i=0;i<centroids.size();i++)
        {
            distances[i] = getEuclideanDistance(data,centroids[i]);
        }
        int minIndex = 0;
        for(int i=1;i<centroids.size();i++)
        {
            if(distances[minIndex]>distances[i])
                minIndex = i;
        }
        return minIndex;
    }
    void runKMeans(vector<vector<float> > data,int K,vector<vector<float> >& centroids)
    {
        set<int> s;
        int index;
        srand(time(NULL));
        for(int i=0;i<K;i++)
        {
            index = rand()%(data.size());
            if(s.find(index)==s.end())
            {
                s.insert(index);
                centroids.push_back(data[index]);
            }
            else
            {
                i--;
            }
        }
        s.clear();
        int iteration = 0;
        int* classLabels = new int[data.size()];
        int temp;
        int numberOfChanges = 0;
        int* count = new int[K];
        while(iteration<100)
        {
            numberOfChanges = 0;
            for(int i=0;i<data.size();i++)
            {
                temp = getNearestCentroid(data[i],centroids);
                if(temp!=classLabels[i])
                    numberOfChanges++;
                classLabels[i] = temp;
            }
            centroids.clear();
            for(int i=0;i<K;i++)
            {
                vector<float> vect(data[0].size(),0.0);
                centroids.push_back(vect);
            }
            memset(count,0,sizeof(int)*K);
            for(int i=0;i<data.size();i++)
            {
                for(int j=0;j<data[i].size();j++)
                {
                    count[classLabels[i]]++;
                    centroids[classLabels[i]][j] = (centroids[classLabels[i]][j]*(count[classLabels[i]]-1) + data[i][j])/count[classLabels[i]];
                }
            }
            iteration++;
            cout<<"Iteration: "<<iteration<<" changes = "<<numberOfChanges<<endl;
            
            for(int i=0;i<K;i++)
            {
                if(count[i]==0)
                {
                    vector<float> newCentroid;
                    cout<<"Found an empty cluster"<<endl;
                    for(int j=0;j<data[0].size();j++)
                    {
                        newCentroid.push_back(rand());
                    }
                    centroids[i] = newCentroid;
                }
            }
            if(numberOfChanges<=100)
                break;
        }
        delete classLabels;
        delete count;
    }
    
    // PART 1.3
    void test(const Dataset &filenames, const string COLOR_MODE="COLOR")
    {
        //create test.dat file
        ofstream train_output_File;
        train_output_File.open ("svm_multiclass/test.dat");
        int class_counter = 1;
        
        CImg<double> class_vectors;
        if (COLOR_MODE == "GRAY")
        {
            CImg<double> class_vectors_temp(size*size*1, filenames.size(), 1);
            class_vectors = class_vectors_temp;
        }
        else
        {
            CImg<double> class_vectors_temp(size*size*3, filenames.size(), 1); // COLOR
            class_vectors = class_vectors_temp;
        }
        
        if(featureType == "baseline"){
            
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                // convert each image to be a row of this "model" image
                for(int i=0; i<c_iter->second.size(); i++){
                    
                    CImg<double> temp = extract_features(c_iter->second[i].c_str(), COLOR_MODE);
                    class_vectors = class_vectors.draw_image(0, i, 0, 0, temp);
                    writeTrainFile(train_output_File, class_counter, temp);
                }
                
            }
        }
        else if(featureType == "eigen")
        {}
        else if(featureType == "haar")
        {
            vector <HaarFilter> haarFilters = readHaarFiltersFromFile();
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                cout << "class " << class_counter << " started"  << endl;
                for(int i=0; i<c_iter->second.size(); i++)
                {
                    CImg<double> input_image(c_iter->second[i].c_str());
                    CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
                    gray = gray.resize(25,25,1,1);
                    CImg<double> int_image(gray.width(), gray.height(), 1, 1);
                    integral_image(gray, int_image);
                    vector<double> feature = haar_feature_all_filters(int_image, haarFilters, 2, 2);
                    CImg<double> temp(feature.size());
                    if(class_counter == 1&& i == 0)
                        cout << "haar feature Size = " << feature.size() << endl;
                    for(int p=0;p<feature.size();p++)
                        temp(p) = feature[p];
                    writeTrainFile(train_output_File, class_counter, temp);
                }
            }
        }
        else if(featureType == "bow")
        {}
        else if(featureType == "deep")
        {
            for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
            {
                for(int i=0; i< c_iter->second.size(); i++)
                {
                    char cwd[1024];
                    chdir("/path/to/change/directory/to");
                    getcwd(cwd, sizeof(cwd));
                    //printf("Train-> Current working dir: %s\n", cwd);
                    
                    string directory = "overfeat/bin/linux_64";
                    chdir (directory.c_str());
                    //cout<<"chdir():"<<ret<<endl;
                    
                    chdir("/path/to/change/directory/to");
                    getcwd(cwd, sizeof(cwd));
                    //printf("Train-> Current working dir: %s\n", cwd);
                    
                    char image_name[300] = "";
                    strcat(image_name, "../../../");
                    strcat(image_name, c_iter->second[i].c_str());
                    CImg<double> temp(image_name);
                    temp.resize(DEEP_SUBSAMPLE_SIZE,DEEP_SUBSAMPLE_SIZE,1,3);
                    strcpy (image_name,"../../../resized_image.jpg");
                    temp.save(image_name);
                    
                    char command[300] = "";
                    strcat(command, "./overfeat -f ");
                    strcat(command, image_name);
                    strcat(command, " > ../../../feature.txt");
                    cout<<command<<endl;
                    system(command); // command overfeat
                    strcpy (command,"");
                    strcpy (image_name,"");
                    
                    //read feature file
                    ifstream feature_File;
                    string line, word;
                    feature_File.open ("../../../feature.txt");
                    
                    if (feature_File.is_open())
                    {
                        getline (feature_File,line);
                        CImg<double> temp(DEEP_FEATURE_SIZE);
                        for (int k=0; feature_File >> word; k++ )
                            temp[k] = (double)atof(word.c_str());
                        
                        class_vectors = class_vectors.draw_image(0, i, 0, 0, temp);
                        writeTrainFile(train_output_File, class_counter, temp);
                        cout<<"Size:"<<temp.size()<<endl;
                        feature_File.close();
                    }
                    //break;
                    
                }
                //break;
                
            }
            chdir("../../..");
            
        }
        
        train_output_File.close();
    }
    
    // Part-1 don't need this. So empty implementation here.
    virtual string classify(const string &filename, const string COLOR_MODE="COLOR"){return NULL;}
    // Part-1 don't need this. So empty implementation here.
    virtual void load_model(){}
    
protected:
    // extract features from an image, which in this case just involves resampling and
    // rearranging into a vector of pixel data.
    CImg<double> extract_features(const string &filename, const string COLOR_MODE)
    {
        if (COLOR_MODE == "GRAY")
            return (CImg<double>(filename.c_str())).resize(size,size,1,1).unroll('x');
        else
            return (CImg<double>(filename.c_str())).resize(size,size,1,3).unroll('x'); //COLOR
    }
    
    // Write file in appropiate formate of API
    void writeTrainFile(ofstream &train_output_File, const int class_counter, const CImg<double> class_vectors){
        train_output_File<<class_counter<<" ";
        
        for (int i=0; i<class_vectors.size(); ++i) {
            
            if (i < class_vectors.size()-1)
                train_output_File<<(i+1)<<":"<<class_vectors[i]<<" ";
            else
                train_output_File<<(i+1)<<":"<<class_vectors[i]<<"\n";
        }
        
    }
    
    static const int size=40;  // subsampled image resolution
    
    static const int DEEP_SUBSAMPLE_SIZE=231;  // Deep Feature subsampled image resolution
    static const int DEEP_FEATURE_SIZE=4096;  // For Subsample 231 is fixed, otherwise multiple of 4096
    // map<string, CImg<double> > models; // trained models
};

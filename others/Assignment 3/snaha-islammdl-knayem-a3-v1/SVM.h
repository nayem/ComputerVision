class SVM : public Classifier
{
public:
    SVM(const vector<string> &_class_list) : Classifier(_class_list) {}
    
    // SVM training. All this does is read in all the images, resize
    // them to a common size, and dump them as vectors to a file
    // PART 1.2
    virtual void train(const Dataset &filenames, const string COLOR_MODE="COLOR")
    {
        //create train.dat file
        ofstream train_output_File;
        train_output_File.open ("svm_multiclass/train.dat");
        
        int class_counter = 1;
        
        for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
        {
            
            CImg<double> class_vectors;
            if (COLOR_MODE == "GRAY"){
                CImg<double> class_vectors_temp(size*size*1, filenames.size(), 1);
                class_vectors = class_vectors_temp;
            }
            else{
                CImg<double> class_vectors_temp(size*size*3, filenames.size(), 1); // COLOR
                class_vectors = class_vectors_temp;
            }
            
            // convert each image to be a row of this "model" image
            for(int i=0; i<c_iter->second.size(); i++){
                
                CImg<double> temp = extract_features(c_iter->second[i].c_str(), COLOR_MODE);
                class_vectors = class_vectors.draw_image(0, i, 0, 0, temp);
                writeTrainFile(train_output_File, class_counter, temp);
            }
            //class_vectors.save_png(("nn_model." + c_iter->first + ".png").c_str());
        }
        
        train_output_File.close();
    }

    // PART 1.3
    void test(const Dataset &filenames, const string COLOR_MODE="COLOR")
    {
        //create test.dat file
        ofstream train_output_File;
        train_output_File.open ("svm_multiclass/test.dat");
        
        int class_counter = 1;
        
        for(Dataset::const_iterator c_iter=filenames.begin(); c_iter != filenames.end(); ++c_iter, ++class_counter)
        {
            
            CImg<double> class_vectors;
            if (COLOR_MODE == "GRAY"){
                CImg<double> class_vectors_temp(size*size*1, filenames.size(), 1);
                class_vectors = class_vectors_temp;
            }
            else{
                CImg<double> class_vectors_temp(size*size*3, filenames.size(), 1); // COLOR
                class_vectors = class_vectors_temp;
            }
            
            // convert each image to be a row of this "model" image
            for(int i=0; i<c_iter->second.size(); i++){
                
                CImg<double> temp = extract_features(c_iter->second[i].c_str(), COLOR_MODE);
                class_vectors = class_vectors.draw_image(0, i, 0, 0, temp);
                writeTrainFile(train_output_File, class_counter, temp);
            }

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
    
    static const int size=150;  // subsampled image resolution
    // map<string, CImg<double> > models; // trained models
};

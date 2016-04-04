/* ==============================================================
// MoCo:  Model-convolution  algorithm   for   3D   mitochondrial
// networks.  In  this  algorithm,  the  graph  representation of
// a mitochondria  is  converted back  into  the microscopy image
// of this  mitochondria by convolving a voxelized version of the
// graph   with  a microscope point-spread-function.
//
// Matheus P. Viana, UC Irvine, vianamp@gmail.com   -  15.10.2014
//                                                     23.05.2015
// --------------------------------------------------------------
// [1] Model Convolution:  A  Computational  Approach to  Digital
// Image Interpretation by MELISSA K. GARDNER et al.
//
// [2] Analyzing  fluorescence microscopy  images with  ImageJ by
// Peter Bankhead Chapter: Simulating Image Formation.
// ============================================================*/

#include "MoCo.h"

int PoissonGen(double mu) {
    int k = 0;
    double p = 1.0, L = exp(-mu);
    do {
        k++;
        p *= (double)(rand()) / RAND_MAX;
    } while (p>L);
    return k;
}


int ScanFolderForThisExtension(std::string RootFolder, const char ext[], std::vector<std::string> *List) {
    DIR *dir;
    int ext_p;
    struct dirent *ent;
    std::string _dir_name;
    if ((dir = opendir(RootFolder.c_str())) != NULL) {
      while ((ent = readdir (dir)) != NULL) {
        _dir_name = std::string(ent->d_name);
        ext_p = (int)_dir_name.find(std::string(ext));
        if (ext_p > 0) {
            #ifdef DEBUG
                printf("File found: %s\n",_dir_name.c_str());
            #endif
            List -> push_back(RootFolder+_dir_name.substr(0,ext_p));
        }
      }
      closedir (dir);
    } else {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int ExportMaxProjection(vtkImageData *Image, const char FileName[]) {

    #ifdef DEBUG
        printf("Saving max projection...\n");
    #endif

    int *Dim = Image -> GetDimensions();
    vtkSmartPointer<vtkImageData> MaxP = vtkSmartPointer<vtkImageData>::New();
    MaxP -> SetDimensions(Dim[0],Dim[1],1);
    vtkIdType N = Dim[0] * Dim[1];

    vtkSmartPointer<vtkFloatArray> MaxPArray = vtkSmartPointer<vtkFloatArray>::New();
    MaxPArray -> SetNumberOfComponents(1);
    MaxPArray -> SetNumberOfTuples(N);

    int x, y, z;
    double v, vproj;
    for (x = Dim[0]; x--;) {
        for (y = Dim[1]; y--;) {
            vproj = 0;
            for (z = Dim[2]; z--;) {
                v = Image -> GetScalarComponentAsFloat(x,y,z,0);
                vproj = (v > vproj) ? v : vproj;
            }
            MaxPArray -> SetTuple1(MaxP->FindPoint(x,y,0),vproj);
        }
    }
    MaxPArray -> Modified();

    MaxP -> GetPointData() -> SetScalars(MaxPArray);

    vtkSmartPointer<vtkTIFFWriter> TIFFWriter = vtkSmartPointer<vtkTIFFWriter>::New();
    TIFFWriter -> SetFileName(FileName);
    TIFFWriter -> SetInputData(MaxP);
    TIFFWriter -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif

    return EXIT_SUCCESS;
}

int ExportSlice(vtkImageData *Image, const char FileName[], int z) {

    #ifdef DEBUG
        printf("Saving slice...\n");
    #endif

    int *Dim = Image -> GetDimensions();
    vtkSmartPointer<vtkImageData> MaxP = vtkSmartPointer<vtkImageData>::New();
    MaxP -> SetDimensions(Dim[0],Dim[1],1);
    vtkIdType N = Dim[0] * Dim[1];

    vtkSmartPointer<vtkFloatArray> MaxPArray = vtkSmartPointer<vtkFloatArray>::New();
    MaxPArray -> SetNumberOfComponents(1);
    MaxPArray -> SetNumberOfTuples(N);

    int x, y;
    double v, vproj;
    for (x = Dim[0]; x--;) {
        for (y = Dim[1]; y--;) {
            v = Image -> GetScalarComponentAsFloat(x,y,z,0);
            MaxPArray -> SetTuple1(MaxP->FindPoint(x,y,0),v);
        }
    }
    MaxPArray -> Modified();

    MaxP -> GetPointData() -> SetScalars(MaxPArray);

    vtkSmartPointer<vtkTIFFWriter> TIFFWriter = vtkSmartPointer<vtkTIFFWriter>::New();
    TIFFWriter -> SetFileName(FileName);
    TIFFWriter -> SetInputData(MaxP);
    TIFFWriter -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif

    return EXIT_SUCCESS;
}

int ExportTIFFSeq(vtkImageData *Image, const char FileName[]) {

    #ifdef DEBUG
        printf("Saving TIFF sequence...\n");
    #endif

    vtkSmartPointer<vtkTIFFWriter> Writer = vtkSmartPointer<vtkTIFFWriter>::New();
    Writer -> SetInputData(Image);
    Writer -> SetFilePattern("%s%04d.tif");
    Writer -> SetFilePrefix(FileName);
    Writer -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif

    return EXIT_SUCCESS;
}

int SaveImageData(vtkImageData *Image, const char FileName[]) {
    #ifdef DEBUG
        printf("Saving ImageData file...\n");
    #endif

    vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    writer -> SetInputData(Image);
    writer -> SetFileType(VTK_BINARY);
    writer -> SetFileName(FileName);
    writer -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif

    return EXIT_SUCCESS;
}

int SavePolyData(vtkPolyData *PolyData, const char FileName[]) {

    #ifdef DEBUG
        printf("Saving PolyData from XYZ list...\n");
    #endif

    #ifdef DEBUG
        printf("\t#Points in PolyData file: %llu.\n",(vtkIdType)PolyData->GetNumberOfPoints());
    #endif

    vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    Writer -> SetFileType(VTK_BINARY);
    Writer -> SetFileName(FileName);
    Writer -> SetInputData(PolyData);
    Writer -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif

    return EXIT_SUCCESS;
}

int Voxelization(std::string _skell_path_prefix, _MoCoControl MoCoControl) {

    #ifdef DEBUG
        printf("Reading skeleton file: %s.vtk\n",_skell_path_prefix.c_str());
    #endif

    vtkSmartPointer<vtkPolyDataReader> PolyReader = vtkSmartPointer<vtkPolyDataReader>::New();
    PolyReader -> SetFileName((_skell_path_prefix+".vtk").c_str());
    PolyReader -> Update();

    vtkPolyData *Skell = PolyReader -> GetOutput();

    vtkPoints *Points = Skell -> GetPoints();

    double *Bds = Points -> GetBounds();

    double r[3];
    vtkIdType idp;
    for (idp = 0; idp < Points -> GetNumberOfPoints(); idp++) {
        Points -> GetPoint(idp,r);
        Points -> SetPoint(idp,r[0],Bds[3]-r[1],r[2]);
    }
    Points -> Modified();

    Bds = Points -> GetBounds();

    #ifdef DEBUG
        printf("\t#Points = %lld\n",Points->GetNumberOfPoints());
        printf("\tBounds\n");
        printf("\t\t X: [%1.2f,%1.2f]um\n",Bds[0],Bds[1]);
        printf("\t\t y: [%1.2f,%1.2f]um\n",Bds[2],Bds[3]);
        printf("\t\t z: [%1.2f,%1.2f]um\n",Bds[4],Bds[5]);
    #endif

    int nx = round((Bds[1]-Bds[0])/MoCoControl.GetXYSpacing()); //70
    int ny = round((Bds[3]-Bds[2])/MoCoControl.GetXYSpacing()); //93
    int nz = round((Bds[5]-Bds[4])/MoCoControl.GetXYSpacing());

    double _off_x = 0.50 * (200-nx);
    double _off_y = 0.50 * (200-ny);
    double _off_z = 30.0;

    int Dim[3];
    Dim[0] = (int)(2*_off_x + nx);
    Dim[1] = (int)(2*_off_y + ny);
    Dim[2] = (int)(2*_off_z + nz);

    #ifdef DEBUG
        printf("\t#Images in X = %d\n",Dim[0]);
        printf("\t#Images in Y = %d\n",Dim[1]);
        printf("\t#Images in Z = %d\n",Dim[2]);
    #endif

    long int N = Dim[0]*Dim[1]*Dim[2];

    vtkSmartPointer<vtkFloatArray> Scalars = vtkSmartPointer<vtkFloatArray>::New();
    Scalars -> SetNumberOfComponents(1);
    Scalars -> SetNumberOfTuples(N);
    Scalars -> FillComponent(0,0.0);

    vtkSmartPointer<vtkImageData> Image = vtkSmartPointer<vtkImageData>::New();
    Image -> SetDimensions(Dim);
    Image -> GetPointData() -> SetScalars(Scalars);

    int x, y, z;
    vtkCell *Line;
    int ir1[3], ir2[3];
    double r1[3], r2[3], vec[3], norm_vec, d, w1, w2;

    vtkSmartPointer<vtkDataArray> Width = Skell -> GetPointData() -> GetArray("Width");

    if (!Width) {
        printf("No diameter found...\n");
    }

    for (long int line = Skell->GetNumberOfCells(); line--;) {

        Line = Skell -> GetCell(line);
        for (idp = Line->GetNumberOfPoints()-1; idp--;) {
            Points -> GetPoint(Line->GetPointId(idp+1),r1);
            Points -> GetPoint(Line->GetPointId(idp  ),r2);
            if (Width) {
                w1 = Width -> GetTuple1(Line->GetPointId(idp+1));
                w2 = Width -> GetTuple1(Line->GetPointId(idp  ));
            } else {
                w1 = w2 = 1.0;
            }
            ir1[0] = (int)(_off_x + (r1[0]-Bds[0])/MoCoControl.GetXYSpacing());
            ir1[1] = (int)(_off_y + (r1[1]-Bds[2])/MoCoControl.GetXYSpacing());
            ir1[2] = (int)(_off_z + (r1[2]-Bds[4])/MoCoControl.GetXYSpacing());

            ir2[0] = (int)(_off_x + (r2[0]-Bds[0])/MoCoControl.GetXYSpacing());
            ir2[1] = (int)(_off_y + (r2[1]-Bds[2])/MoCoControl.GetXYSpacing());
            ir2[2] = (int)(_off_z + (r2[2]-Bds[4])/MoCoControl.GetXYSpacing());

            vec[0] = ir2[0] - ir1[0];
            vec[1] = ir2[1] - ir1[1];
            vec[2] = ir2[2] - ir1[2];
            d = sqrt(pow(vec[0],2)+pow(vec[1],2)+pow(vec[2],2));
            if ((int)d > 1) {
                vec[0] /= d;
                vec[1] /= d;
                vec[2] /= d;
                for (int n = 0; n < (int)d; n++) {
                    x = ir1[0] + (int)(n*vec[0]);
                    y = ir1[1] + (int)(n*vec[1]);
                    z = ir1[2] + (int)(n*vec[2]);
                    Image -> SetScalarComponentFromDouble(x,y,z,0,-0.5*(w1+w2));
                }
            } else {
                Image -> SetScalarComponentFromDouble(ir1[0],ir1[1],ir1[2],0,-0.5*(w1+w2)*1.25); /*Factor 1.25 compensate the underestimation of tubules width*/
            }
        }
    }

    double wrange[2];
    Skell -> GetScalarRange(wrange);

    vtkIdType id;
    double d0, v;
    double dmax = 0.5 * wrange[1] / MoCoControl.GetXYSpacing();
    int dmax_int = round(dmax);

    #ifdef DEBUG
        printf("\tdmax = %1.4f\n",dmax);
        printf("\tdmax_int = %d\n",dmax_int);
    #endif

    std::vector<vtkIdType> Centers;

    for (id = N; id--;) {
        Image -> GetPoint(id,r);
        d0 = -(0.5/MoCoControl.GetXYSpacing()) * Image -> GetPointData() -> GetScalars() -> GetTuple1(id);
        if ( d0 > 0 ) {
            Centers.push_back(id);
            for ( x = -dmax_int; x <= dmax_int; x++ ) {
                for ( y = -dmax_int; y <= dmax_int; y++ ) {
                    for ( z = -dmax_int; z <= dmax_int; z++ ) {
                        d = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
                        if ( d <= dmax && d < d0 ) {
                            Image -> SetScalarComponentFromDouble(r[0]+x,r[1]+y,r[2]+z,0,d+1);
                        }
                    }
                }
            }
        }
    }

    for (long int p = 0; p < Centers.size(); p++) {
        Image -> GetPointData() -> GetScalars() -> SetTuple1(Centers[p],1);
    }
    Centers.clear();

    // Linear regression based on the tubule radius
    // f(d) =  a - (a-b) * (d-1) / dmax;
    double s;
    double a = 1.0;
    double b = 0.5;
    for (id = N; id--;) {
        d = Image -> GetPointData() -> GetScalars() -> GetTuple1(id);
        if ( d > 0.0 ) {
            s = a - (a-b) * (d-1) / dmax;
            Image -> GetPointData() -> GetScalars() -> SetTuple1(id,65535*s);
        }
    }

    /* IMAGE x PSF CONVOLUTION                                   */

    if ( MoCoControl.ExternalPSF() ) {

        ReaderType::Pointer PSFReader = ReaderType::New();
        PSFReader -> SetFileName( MoCoControl.GetPSFFileName() );

        VTK2ITK::Pointer Filter = VTK2ITK::New();
        Filter -> SetInput(Image);
        Filter -> Update();

        FFTConvolutionType::Pointer FFTConv = FFTConvolutionType::New();
        FFTConv -> SetInput( Filter -> GetOutput() );
        FFTConv -> SetKernelImage( PSFReader -> GetOutput() );
        FFTConv -> SetNormalize( true );

        ITK2VTK::Pointer FilterITKVTK = ITK2VTK::New();
        FilterITKVTK -> SetInput(FFTConv->GetOutput());
        FilterITKVTK -> Update();

        Image -> ShallowCopy(FilterITKVTK -> GetOutput());

    } else {

        vtkSmartPointer<vtkImageGaussianSmooth> GaussianKernel = vtkSmartPointer<vtkImageGaussianSmooth>::New();
        GaussianKernel -> SetInputData(Image);
        GaussianKernel -> SetDimensionality(3);
        GaussianKernel -> SetRadiusFactors(MoCoControl.GetPSFRadii(),MoCoControl.GetPSFRadii(),MoCoControl.GetPSFRadii());
        GaussianKernel -> SetStandardDeviations(MoCoControl.GetPSFStdev(),MoCoControl.GetPSFStdev(),2*MoCoControl.GetPSFStdev());
        GaussianKernel -> Update();
        Image -> ShallowCopy(GaussianKernel -> GetOutput());

    }

    /* END OF CONVOLUTION                                         */

    // Define the photon emission at the brightest point
    double range[2];
    Image -> GetScalarRange(range);
    int _exposureTime = MoCoControl.GetExposureTime();
    int _detectorGain = MoCoControl.GetDetectorGain();
    int _detectorOffset = MoCoControl.GetDetectorOffSet();
    int _backgroundPhotons = MoCoControl.GetBackgroundPhotons();
    int _maxPhotonEmission = MoCoControl.GetMaximumPhotonEmission();

    for (id = N; id--;) {
        v = Image -> GetPointData() -> GetScalars() -> GetTuple1(id);
        v = v / range[1] * _maxPhotonEmission + _backgroundPhotons;
        v = _detectorOffset + _detectorGain * _exposureTime * PoissonGen(v);
        v *= 2000.0 / (_detectorOffset + _detectorGain * _exposureTime * _backgroundPhotons);
        Image -> GetPointData() -> GetScalars() -> SetTuple1(id,v);
    }

    vtkSmartPointer<vtkImageResample> Resample = vtkSmartPointer<vtkImageResample>::New();

    Resample -> SetInterpolationModeToLinear();
    Resample -> SetDimensionality(3);
    Resample -> SetInputData(Image);
    Resample -> SetAxisMagnificationFactor(0,1);
    Resample -> SetAxisMagnificationFactor(1,1);
    Resample -> SetAxisMagnificationFactor(2,MoCoControl.GetXYSpacing()/MoCoControl.GetZSpacing());
    Resample -> Update();
    Image -> ShallowCopy(Resample -> GetOutput());
    Image -> SetSpacing(1.0,1.0,1.0);

    vtkSmartPointer<vtkImageCast> Cast = vtkSmartPointer<vtkImageCast>::New();
    Cast -> SetInputData(Image);
    Cast -> SetOutputScalarTypeToUnsignedShort();
    Cast -> Update();

    VTK2ITK_US::Pointer FilterVTKITK = VTK2ITK_US::New();
    FilterVTKITK -> SetInput(Cast->GetOutput());
    FilterVTKITK -> Update();

    WriterUSType::Pointer Writer = WriterUSType::New();
    Writer -> SetInput(FilterVTKITK->GetOutput());
    Writer -> SetFileName((_skell_path_prefix+"_moco.tif").c_str());
    Writer -> Update();

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {     

    srand(getpid());
    std::string Path, PSF;
    _MoCoControl MoCoControl;

    for (int i = 0; i < argc; i++) {
        
        if (!strcmp(argv[i],"-path")) {
            Path = std::string(argv[i+1]);
            if (Path.back() != '/')
                Path = Path + '/';
        }

        if (!strcmp(argv[i],"-psf"))
            MoCoControl.SetPSFFileName(std::string(argv[i+1]));
        
        if (!strcmp(argv[i],"-background_photons"))
            MoCoControl.SetBackgroundPhotons(atoi(argv[i+1]));
        
        if (!strcmp(argv[i],"-snr"))
            MoCoControl.SetSignalToNoiseRatio(atof(argv[i+1]));
        
        if (!strcmp(argv[i],"-exposure_time"))
            MoCoControl.SetExposureTime(atoi(argv[i+1]));
        
        if (!strcmp(argv[i],"-detector_gain"))
            MoCoControl.SetDetectorGain(atoi(argv[i+1]));
        
        if (!strcmp(argv[i],"-detector_offset"))
            MoCoControl.SetDetectorOffSet(atoi(argv[i+1]));
        
        if (!strcmp(argv[i],"-gaussian")) {
            MoCoControl.SetPSFStdev(atof(argv[i+1]));
            MoCoControl.SetPSFRadii(atof(argv[i+2]));
        }
        
        if (!strcmp(argv[i],"-help")) {
            printf("./Moco -path [IMG] -psf [PSF] -background_photons [30] -snr [10.25] -gaussian [2.5 5.0] -exposure_time [200] -detector_gain [1] -detector_offset [100]\n");
            return EXIT_SUCCESS;
        }

    }
    
    if (!MoCoControl.ExternalPSF())
        printf("\tPSF Not Provided. Using Gaussian Kernel.\n");

    std::vector<std::string> Files;
    ScanFolderForThisExtension(Path,".vtk",&Files);

    for (int i = 0; i < Files.size(); i++) {
        Voxelization(Files[i],MoCoControl);
        printf("%s [done]\n",Files[i].c_str());
    }

    return EXIT_SUCCESS;
}

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

#include <list>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/* VTK                                                         */

#include <vtkMath.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkLongArray.h>
#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkTransform.h>
#include <vtkDataObject.h>
#include <vtkPNGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkFloatArray.h>
#include <vtkVectorText.h>
#include <vtkInformation.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkImageResample.h>
#include <vtkStructuredPoints.h>
#include <vtkInformationVector.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkTIFFWriter.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>
#include <vtkPoints.h>

/* ITK                                                         */

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSimpleFilterWatcher.h>
#include <itkVTKImageToImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkConvolutionImageFilter.h>
#include <itkFFTConvolutionImageFilter.h>
#include <itkImageToVTKImageFilter.h>

#define DEBUG

// Add background photons
int _backgroundPhotons = 30;
// Define the photon emission at the brightest point
double _snr = 10.25;
// Define the photon emission at the brightest point
int _maxPhotonEmission = (int) ((_snr-1.0) * _backgroundPhotons);
// Gaussian PSF
double _psfSigma = 2.5;
double _psfRadii = 5.0;
// Simulate photon noise and Multiply by the exposure time
int _exposureTime = 200;
// Simulate the detector gain
// (note this should really add Poisson noise too!)    
int _detectorGain = 1;
// Simulate the detector offset
int _detectorOffset = 100;

typedef float PType;
typedef itk::Image< PType, 3 > IType;
typedef itk::VTKImageToImageFilter<IType> VTK2ITK;
typedef itk::ImageToVTKImageFilter<IType> ITK2VTK;
typedef itk::ImageFileReader< IType > ReaderType;
typedef itk::FFTConvolutionImageFilter< IType > FFTConvolutionType;

int ssdx[26] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
int ssdy[26] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
int ssdz[26] = { 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1};

// Displays parameters that can be changed in this routine
void _help();

// Save VTK ImageData files
void SaveImageData(vtkImageData *Image, const char FileName[]);

// Save VTK PolyData files
void SavePolyData(vtkPolyData *PolyData, const char FileName[]);

// Export the max projection of a given stack
void ExportMaxProjection(vtkImageData *Image, const char FileName[]);

// Export a slice of a given stack
void ExportSlice(vtkImageData *Image, const char FileName[], int z);

// Export a stck a a sequence of TIFF files
// This routine would be ideally replaced by another in which was
// possible to export the whole stack a a single multi-paged TIFF
// file. However, this is not possible at this moment because there
// is a bug in the current version of vtkTIFFWriter class.
void ExportTIFFSeq(vtkImageData *Image, const char FileName[]);

// Main routine in which the vtkPolyData structure representing
// the mitochondrial network is voxelized and convolved with a
// theoretical point-spread-function (PSF).
// Right now we are assuming a Guassian PSF.
int Voxelization(const char _skell_path_prefix[], ReaderType *PSFReader);

// Poisson random number generator
int PoissonGen(double mu);

/* ================================================================
   AUXILIAR
=================================================================*/

int PoissonGen(double mu) {
    int k = 0;
    double p = 1.0, L = exp(-mu);
    do {
        k++;
        p *= (double)(rand()) / RAND_MAX;
    } while (p>L);
    return k;
}


/* ================================================================
   I/O ROUTINES
=================================================================*/

void ExportMaxProjection(vtkImageData *Image, const char FileName[]) {

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

}

void ExportSlice(vtkImageData *Image, const char FileName[], int z) {

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

}

void ExportTIFFSeq(vtkImageData *Image, const char FileName[]) {

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

}

void SaveImageData(vtkImageData *Image, const char FileName[]) {
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
}

void SavePolyData(vtkPolyData *PolyData, const char FileName[]) {

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
}

void ExportGraphProperties(vtkPolyData *Skell, const char FileName[]) {
    vtkCell *Line;
    double r1[3], r2[3], length;
    vtkIdType i, k, n, line, N = 0;
    for (line = Skell -> GetNumberOfCells(); line--;) {
        Line = Skell -> GetCell(line);
        n = Line -> GetNumberOfPoints() - 1;
        for (k = n; k--;) {
            Skell -> GetPoints() -> GetPoint(Line->GetPointId(k+1),r1);
            Skell -> GetPoints() -> GetPoint(Line->GetPointId(k  ),r2);
            length += sqrt(pow(r2[0]-r1[0],2)+pow(r2[1]-r1[1],2)+pow(r2[2]-r1[2],2));
        }
        i = Line -> GetPointId(0); N = (i > N) ? i : N;
        i = Line -> GetPointId(n); N = (i > N) ? i : N;
    }
    N++;
    FILE *f = fopen(FileName,"w");
    fprintf(f,"N\tE\tL\n");
    fprintf(f,"%d\t%d\t%1.3f\n",(int)N,(int)Skell->GetNumberOfCells(),length);
    fclose(f);
}

/* ================================================================
   MAIN ROUTINE
=================================================================*/

int Voxelization(const char _skell_path_prefix[], ReaderType *PSFReader) {

    #ifdef DEBUG
        printf("Reading skeleton file: %s.vtk\n",_skell_path_prefix);
    #endif

    char FileName[256];
    sprintf(FileName,"%s.vtk",_skell_path_prefix);

    vtkSmartPointer<vtkPolyDataReader> PolyReader = vtkSmartPointer<vtkPolyDataReader>::New();
    PolyReader -> SetFileName(FileName);
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

    double dxy = 0.056;
    double dz  = 0.200;
    double _rad = 0.15;

    int nx = round((Bds[1]-Bds[0])/dxy); //70
    int ny = round((Bds[3]-Bds[2])/dxy); //93
    int nz = round((Bds[5]-Bds[4])/dxy);

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

    for (long int line = Skell->GetNumberOfCells(); line--;) {

        Line = Skell -> GetCell(line);
        for (idp = Line->GetNumberOfPoints()-1; idp--;) {
            Points -> GetPoint(Line->GetPointId(idp+1),r1);
            Points -> GetPoint(Line->GetPointId(idp  ),r2);
            w1 = Width -> GetTuple1(Line->GetPointId(idp+1));
            w2 = Width -> GetTuple1(Line->GetPointId(idp  ));

            ir1[0] = (int)(_off_x + (r1[0]-Bds[0])/dxy);
            ir1[1] = (int)(_off_y + (r1[1]-Bds[2])/dxy);
            ir1[2] = (int)(_off_z + (r1[2]-Bds[4])/dxy);

            ir2[0] = (int)(_off_x + (r2[0]-Bds[0])/dxy);
            ir2[1] = (int)(_off_y + (r2[1]-Bds[2])/dxy);
            ir2[2] = (int)(_off_z + (r2[2]-Bds[4])/dxy);

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
                Image -> SetScalarComponentFromDouble(ir1[0],ir1[1],ir1[2],0,-0.5*(w1+w2));
            }
        }
    }

    double wrange[2];
    Skell -> GetScalarRange(wrange);

    vtkIdType id;
    double d0, v;
    double dmax = 0.5 * wrange[1] / dxy;
    int dmax_int = round(dmax);

    #ifdef DEBUG
        printf("\tdmax = %1.4f\n",dmax);
        printf("\tdmax_int = %d\n",dmax_int);
    #endif

    std::vector<vtkIdType> Centers;

    for (id = N; id--;) {
        Image -> GetPoint(id,r);
        d0 = -(0.5/dxy) * Image -> GetPointData() -> GetScalars() -> GetTuple1(id);
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
            Image -> GetPointData() -> GetScalars() -> SetTuple1(id,100*s);
        }
    }

    /* IMAGE x PSF CONVOLUTION                                   */

    if ( PSFReader ) {

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
        GaussianKernel -> SetRadiusFactors(_psfRadii,_psfRadii,_psfRadii);
        GaussianKernel -> SetStandardDeviations(_psfSigma,_psfSigma,2*_psfSigma);
        GaussianKernel -> Update();
        Image -> ShallowCopy(GaussianKernel -> GetOutput());

    }

    /* END OF CONVOLUTION                                         */

    // Define the photon emission at the brightest point
    int _maxPhotonEmission = (int) ((_snr-1.0) * _backgroundPhotons);

    double range[2];
    Image -> GetScalarRange(range);
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
    Resample -> SetAxisMagnificationFactor(2,dxy/dz);
    Resample -> Update();
    Image -> ShallowCopy(Resample -> GetOutput());

    Image -> SetSpacing(1.0,1.0,1.0);

    VTK2ITK::Pointer FilterVTKITK = VTK2ITK::New();
    FilterVTKITK -> SetInput(Image);
    FilterVTKITK -> Update();

    sprintf(FileName,"%s_moco.tif",_skell_path_prefix);

    typedef itk::ImageFileWriter< IType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer -> SetInput( FilterVTKITK->GetOutput() );
    writer -> SetFileName( FileName );
    writer -> Update();

    return 0;
}

int main(int argc, char *argv[]) {     

    srand(getpid());

    int i;
    char _imgpath[128];
    char _psfpath[128];
    long int _cc_value;
    sprintf(_imgpath,"");
    sprintf(_psfpath,"blank");

    // Collecting input parameters
    for (i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            sprintf(_imgpath,"%s//",argv[i+1]);
        }
        if (!strcmp(argv[i],"-psf")) {
            sprintf(_psfpath,"%s",argv[i+1]);
        }
        if (!strcmp(argv[i],"-background_photons")) {
            _backgroundPhotons = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i],"-snr")) {
            _snr = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-gaussian")) {
            _psfSigma = atof(argv[i+2]);
            _psfRadii = atof(argv[i+3]);
        }
        if (!strcmp(argv[i],"-exposure_time")) {
            _exposureTime = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-detector_gain")) {
            _detectorGain = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-detector_offset")) {
            _detectorOffset = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-help")) {
            printf("./Moco -path [IMG] -psf [PSF] -background_photons [30] -snr [10.25] -gaussian [2.5 5.0] -exposure_time [200] -detector_gain [1] -detector_offset [100]\n");
            return 1;
        }        
    }

    // Loading the PSF
    ReaderType::Pointer psf_reader = NULL;
    if (!strcmp(_psfpath,"blank")) {
        #ifdef DEBUG
            printf("\tPSF Not Provided. Using Gaussian Kernel.\n");
        #endif
    } else {
        #ifdef DEBUG
            printf("\tLoading PSF...\n");
        #endif
        psf_reader = ReaderType::New();
        psf_reader -> SetFileName( _psfpath );
    }

    // Generating list of files to run
    char _cmd[256];
    sprintf(_cmd,"ls %s*_skeleton.vtk | sed -e 's/.vtk//' > %smoco.files",_imgpath,_imgpath);
    system(_cmd);

    // Voxelization
    char _skell_path_prefix[256];
    char _skellist_path_filename[256];
    sprintf(_skellist_path_filename,"%smoco.files",_imgpath);
    FILE *f = fopen(_skellist_path_filename,"r");
    while (fgets(_skell_path_prefix,256, f) != NULL) {
        _skell_path_prefix[strcspn(_skell_path_prefix, "\n" )] = '\0';
        Voxelization(_skell_path_prefix,psf_reader);
        printf("%s [done]\n",_skell_path_prefix);
    }
    fclose(f);

    return 1;
}

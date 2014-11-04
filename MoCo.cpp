// ==============================================================
// MoCo:  Model-convolution  algorithm   for   3D   mitochondrial
// networks.  In  this  algorithm,  the  graph  representation of
// a mitochondria  is  converted back  into  the microscopy image
// of this  mitochondria by convolving a voxelized version of the
// graph   with  a theoretical   point-spread-function   of   the
// microscope.
//
// Matheus P. Viana, UC Irvine, vianamp@gmail.com   -  15.10.2014
// --------------------------------------------------------------
// [1] Model Convolution:  A  Computational  Approach to  Digital
// Image Interpretation by MELISSA K. GARDNER et al.
//
// [2] Analyzing  fluorescence microscopy  images with  ImageJ by
// Peter Bankhead Chapter: Simulating Image Formation.
// ==============================================================

#include <list>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

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
int Voxelization(const char _skell_path_prefix[]);

// Poisson random number generator
int PoissonGen(double mu);

/* ================================================================
   HELP
=================================================================*/

void _help() {
    printf("\nMoCo: Model Convolution v1.0\n");
    printf("----------------------------\n");
    printf("\t>Background photons:\t-background_photons\t{30}\n");
    printf("\t>Signal-to-noise ratio:\t-snr\t\t\t{10}\n");
    printf("\t>Point-spread-function:\t-psf\t\t\t{gaussian 2.5 5.0}\n");
    printf("\t>Exposure time:\t\t-exposure_time\t\t{200}\n");
    printf("\t>Detector gain:\t\t-detector_gain\t\t{1}\n");
    printf("\t>Detector offset:\t-detector_offset\t{100}\n");
}

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


/* ================================================================
   MAIN ROUTINE
=================================================================*/

int Voxelization(const char _skell_path_prefix[]) {

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

    int nx = round((Bds[1]-Bds[0])/dxy);
    int ny = round((Bds[3]-Bds[2])/dxy);
    int nz = round((Bds[5]-Bds[4])/dxy);

    double _off_x = 0.0;    //because bounds in _x and _y don't have the same length
    double _off_y = 0.0;
    double _off_z = 5.0/dz;
    double _off_xyz = 5.0;  //global offset

    if (nx>ny) {
        _off_x = 0.0;
        _off_y = 0.5 * (nx-ny);
        ny = nx;
    } else {
        _off_y = 0.0;
        _off_x = 0.5 * (ny-nx);
        nx = ny;        
    }

    int Dim[3];
    Dim[0] = (int)(2*(_off_xyz+_off_x) + nx);
    Dim[1] = (int)(2*(_off_xyz+_off_y) + ny);
    Dim[2] = (int)(2*(_off_xyz+_off_z) + nz);

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
    vtkIdType idp;
    vtkCell *Line;
    int ir1[3], ir2[3];
    double r1[3], r2[3], vec[3], norm_vec, d, w1, w2;
    for (long int line = Skell->GetNumberOfCells(); line--;) {
        Line = Skell -> GetCell(line);
        for (idp = Line->GetNumberOfPoints()-1; idp--;) {
            Points -> GetPoint(Line->GetPointId(idp+1),r1);
            Points -> GetPoint(Line->GetPointId(idp  ),r2);
            w1 = Skell -> GetPointData() -> GetScalars() -> GetTuple1(Line->GetPointId(idp+1));
            w2 = Skell -> GetPointData() -> GetScalars() -> GetTuple1(Line->GetPointId(idp  ));

            ir1[0] = (int)(_off_xyz + _off_x + (r1[0]-Bds[0])/dxy);
            ir1[1] = (int)(_off_xyz + _off_y + (r1[1]-Bds[2])/dxy);
            ir1[2] = (int)(_off_xyz + _off_z + (r1[2]-Bds[4])/dxy);

            ir2[0] = (int)(_off_xyz + _off_x + (r2[0]-Bds[0])/dxy);
            ir2[1] = (int)(_off_xyz + _off_y + (r2[1]-Bds[2])/dxy);
            ir2[2] = (int)(_off_xyz + _off_z + (r2[2]-Bds[4])/dxy);

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
    double d0, r[3], v;
    double dmax = 0.5 * wrange[1] / dxy;
    int dmax_int = round(dmax);

    #ifdef DEBUG
        printf("\tdmax = %1.4f\n",dmax);
        printf("\tdmax_int = %d\n",dmax_int);
    #endif

    for (id = N; id--;) {
        Image -> GetPoint(id,r);
        d0 = -(0.5/dxy) * Image -> GetPointData() -> GetScalars() -> GetTuple1(id);
        if ( d0 > 0 ) {
            for ( x = -dmax_int; x <= dmax_int; x++ ) {
                for ( y = -dmax_int; y <= dmax_int; y++ ) {
                    for ( z = -dmax_int; z <= dmax_int; z++ ) {
                        d = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
                        if ( d <= dmax ) {
                            if ( d < d0 ) {
                                Image -> SetScalarComponentFromDouble(r[0]+x,r[1]+y,r[2]+z,0,d);
                            }
                        }
                    }
                }
            }
        }
    }

    // Linear regression based on the tubule radius
    for (id = N; id--;) {
        d = Image -> GetPointData() -> GetScalars() -> GetTuple1(id);
        if ( d != 0.0 ) {
            Image -> GetPointData() -> GetScalars() -> SetTuple1(id,(dmax-d)/(dmax+1.0));
        }
    }

    vtkSmartPointer<vtkImageGaussianSmooth> GaussianKernel = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    GaussianKernel -> SetInputData(Image);
    GaussianKernel -> SetDimensionality(3);
    GaussianKernel -> SetRadiusFactors(_psfRadii,_psfRadii,_psfRadii);
    GaussianKernel -> SetStandardDeviations(_psfSigma,_psfSigma,2*_psfSigma);
    GaussianKernel -> Update();
    Image -> ShallowCopy(GaussianKernel -> GetOutput());

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

    //SaveImageData(Image,"temp.vtk");

    //ExportMaxProjection(Image,"temp.tif");

    //ExportSlice(Image,"tempz.tif",48);

    char cmd[256];
    sprintf(cmd,"mkdir %s",_skell_path_prefix);
    system(cmd);
    sprintf(FileName,"%s//im",_skell_path_prefix);
    ExportTIFFSeq(Image,FileName);

    return 0;
}

int main(int argc, char *argv[]) {     

    srand(getpid());

    int i;
    char _impath[128];
    long int _cc_value;
    sprintf(_impath,"");

    // Collecting input parameters
    for (i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            sprintf(_impath,"%s//",argv[i+1]);
        }
        if (!strcmp(argv[i],"-help")) {
            _help();
            return 1;
        }
        if (!strcmp(argv[i],"-background_photons")) {
            _backgroundPhotons = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i],"-snr")) {
            _snr = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-psf")) {
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
    }

    // Generating list of files to run
    char _cmd[256];
    sprintf(_cmd,"ls %s*_skeleton.vtk | sed -e 's/.vtk//' > %smoco.files",_impath,_impath);
    system(_cmd);

    // Thinning
    char _skell_path_prefix[256];
    char _skellist_path_filename[256];
    sprintf(_skellist_path_filename,"%smoco.files",_impath);
    FILE *f = fopen(_skellist_path_filename,"r");
    while (fgets(_skell_path_prefix,256, f) != NULL) {
        _skell_path_prefix[strcspn(_skell_path_prefix, "\n" )] = '\0';
        Voxelization(_skell_path_prefix);
    }
    fclose(f);

    return 1;
}

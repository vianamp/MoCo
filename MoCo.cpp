// ==============================================================
// MoCo: ??
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
#include <vtkImageResample.h>
#include <vtkDataArray.h>
#include <vtkTransform.h>
#include <vtkDataObject.h>
#include <vtkPNGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkFloatArray.h>
#include <vtkVectorText.h>
#include <vtkInformation.h>
#include <vtkSphereSource.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkInformationVector.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
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

    int ssdx[26] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int ssdy[26] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int ssdz[26] = { 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1};

// Routine used to save an ImageData
void SaveImageData(vtkImageData *Image, const char FileName[]);

// Routine used to save a PolyData
void SavePolyData(vtkPolyData *PolyData, const char FileName[]);

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
        printf("Saving Max projection...\n");
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
        printf("File Saved!\n");
    #endif

}

void ExportSlice(vtkImageData *Image, const char FileName[], int z) {

    #ifdef DEBUG
        printf("Saving Slice...\n");
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
        printf("File Saved!\n");
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
        printf("File Saved!\n");
    #endif

}

void SaveImageData(vtkImageData *Image, const char FileName[]) {
    #ifdef DEBUG
        printf("Saving ImageData File...\n");
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

void _help() {
    printf("\n MoCo: Model Convolution:\n");
    printf("\t -path, -cc, -graph_off, -thinning_only\n");
}

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
    Dim[0] = (int)(2*_off_xyz + nx);
    Dim[1] = (int)(2*_off_xyz + ny);
    Dim[2] = (int)(2*_off_xyz + nz);

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
    double r1[3], r2[3], vec[3], norm_vec, d;
    for (long int line = Skell->GetNumberOfCells(); line--;) {
        Line = Skell -> GetCell(line);
        for (idp = Line->GetNumberOfPoints()-1; idp--;) {
            Points -> GetPoint(Line->GetPointId(idp+1),r1);
            Points -> GetPoint(Line->GetPointId(idp  ),r2);

            ir1[0] = (int)(_off_xyz + _off_x + (r1[0]-Bds[0])/dxy);
            ir1[1] = (int)(_off_xyz + _off_y + (r1[1]-Bds[2])/dxy);
            ir1[2] = (int)(_off_xyz + (r1[2]-Bds[4])/dxy);

            ir2[0] = (int)(_off_xyz + _off_x + (r2[0]-Bds[0])/dxy);
            ir2[1] = (int)(_off_xyz + _off_y + (r2[1]-Bds[2])/dxy);
            ir2[2] = (int)(_off_xyz + (r2[2]-Bds[4])/dxy);

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
                    Image -> SetScalarComponentFromDouble(x,y,z,0,-1);
                }
            } else {
                Image -> SetScalarComponentFromDouble(ir1[0],ir1[1],ir1[2],0,-1);
            }
        }
    }

    vtkIdType id;
    double d0, r[3], v;
    double dmax = _rad / dxy;
    int dmax_int = round(dmax);

    #ifdef DEBUG
        printf("\tdmax = %1.4f\n",dmax);
        printf("\tdmax_int = %d\n",dmax_int);
    #endif

    for (id = N; id--;) {
        Image -> GetPoint(id,r);
        if ( Image -> GetPointData() -> GetScalars() -> GetTuple1(id) == -1) {
            for ( x = -dmax_int; x <= dmax_int; x++ ) {
                for ( y = -dmax_int; y <= dmax_int; y++ ) {
                    for ( z = -dmax_int; z <= dmax_int; z++ ) {
                        d = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
                        if ( d <= dmax ) {
                            d0 = Image -> FindPoint(r[0]+x,r[1]+y,r[2]+z);
                            if ( d < d0 && d > 0.0 ) {
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

    // Add background photons
    int _backgroundPhotons = 30;
    // Define the photon emission at the brightest point
    double _snr = 2.25;
    // Define the photon emission at the brightest point
    int _maxPhotonEmission = (int) ((_snr-1.0) * _backgroundPhotons);
    // Simulate PSF blurring
    double psfSigma = 2.5;
    // Simulate photon noise and Multiply by the exposure time
    int _exposureTime = 200;
    // Simulate the detector gain
    // (note this should really add Poisson noise too!)    
    int _detectorGain = 1;
    // Simulate the detector offset
    int _detectorOffset = 100;

    //==========================
    vtkSmartPointer<vtkImageGaussianSmooth> GaussianKernel = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    GaussianKernel -> SetInputData(Image);
    GaussianKernel -> SetDimensionality(3);
    GaussianKernel -> SetRadiusFactors(5,5,5);
    GaussianKernel -> SetStandardDeviations(psfSigma,psfSigma,psfSigma);
    GaussianKernel -> Update();
    Image -> ShallowCopy(GaussianKernel -> GetOutput());

    double range[2];
    Image -> GetScalarRange(range);
    for (id = N; id--;) {
        v = Image -> GetPointData() -> GetScalars() -> GetTuple1(id);
        v = v / range[1] * _maxPhotonEmission + _backgroundPhotons;
        v = _detectorOffset + _detectorGain * _exposureTime * PoissonGen(v);
        v *= 2000.0 / (_detectorOffset + _detectorGain * _exposureTime * _backgroundPhotons);
        Image -> GetPointData() -> GetScalars() -> SetTuple1(id,v);
    }

/*
    Image -> GetScalarRange(range);
    #ifdef DEBUG
        printf("\tRange after psf = [%1.4f : %1.4f]\n",range[0],range[1]);
    #endif


    run("Add Specified Noise...", "standard="+readStdDev);
*/


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

    SaveImageData(Image,"temp.vtk");

    //ExportMaxProjection(Image,"temp.tif");

    //ExportSlice(Image,"tempz.tif",48);

    ExportTIFFSeq(Image,"im");

}

int main(int argc, char *argv[]) {     

    int i;
    char _impath[128];
    long int _cc_value;
    sprintf(_impath,"");

    // Collecting input parameters
    for (i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            sprintf(_impath,"%s//",argv[i+1]);
        }
        if (!strcmp(argv[i],"-?")) {
            _help();
            return 1;
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

#include <list>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>

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

typedef float PType;
typedef unsigned short USPType;
typedef itk::Image< PType, 3 > IType;
typedef itk::Image< USPType, 3 > USIType;
typedef itk::VTKImageToImageFilter<IType> VTK2ITK;
typedef itk::VTKImageToImageFilter<USIType> VTK2ITK_US;
typedef itk::ImageToVTKImageFilter<IType> ITK2VTK;
typedef itk::ImageFileReader< IType > ReaderType;
typedef itk::ImageFileWriter< USIType > WriterUSType;
typedef itk::FFTConvolutionImageFilter< IType > FFTConvolutionType;

class _MoCoControl;

int ssdx[26] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
int ssdy[26] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
int ssdz[26] = { 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1};

// Save VTK ImageData files
int SaveImageData(vtkImageData *Image, const char FileName[]);

// Save VTK PolyData files
int SavePolyData(vtkPolyData *PolyData, const char FileName[]);

// Export the max projection of a given stack
int ExportMaxProjection(vtkImageData *Image, const char FileName[]);

// Export a slice of a given stack
int ExportSlice(vtkImageData *Image, const char FileName[], int z);

// Export a stck a a sequence of TIFF files
// This routine would be ideally replaced by another in which was
// possible to export the whole stack a a single multi-paged TIFF
// file. However, this is not possible at this moment because there
// is a bug in the current version of vtkTIFFWriter class.
int ExportTIFFSeq(vtkImageData *Image, const char FileName[]);

// Main routine in which the vtkPolyData structure representing
// the mitochondrial network is voxelized and convolved with a
// theoretical point-spread-function (PSF).
// Right now we are assuming a Guassian PSF.
int Voxelization(const char _skell_path_prefix[], ReaderType *PSFReader);

// Poisson random number generator
int PoissonGen(double mu);

int ScanFolderForThisExtension(std::string RootFolder, const char ext[], std::vector<std::string> *List);

int ExportMaxProjection(vtkImageData *Image, const char FileName[]);

int ExportSlice(vtkImageData *Image, const char FileName[], int z);

int ExportTIFFSeq(vtkImageData *Image, const char FileName[]);

int SaveImageData(vtkImageData *Image, const char FileName[]);

int SavePolyData(vtkPolyData *PolyData, const char FileName[]);

int Voxelization(std::string _skell_path_prefix, _MoCoControl MoCoControl);

class _MoCoControl {

    double _snr;
    double _dxy, _dz;
    int _exposureTime;
    int _detectorGain;
    int _detectorOffset;
    int _backgroundPhotons;
    double _psfSigma, _psfRadii;

    bool _external_psf;
    std::string PSFFileName;

public:
    _MoCoControl() {
        SetSignalToNoiseRatio(10.25);
        SetXYSpacing(0.056);
        SetZSpacing(0.2);
        SetExposureTime(200);
        SetDetectorGain(1);
        SetDetectorOffSet(100);
        SetBackgroundPhotons(30);
        SetPSFStdev(2.5);
        SetPSFRadii(5.0);

        _external_psf = FALSE;
    }
    ~_MoCoControl() {

    }
    
      void SetSignalToNoiseRatio(double snr) { _snr = snr;}
    double GetSignalToNoiseRatio()           {return _snr;}

      void SetXYSpacing(double dxy) { _dxy = dxy;}
    double GetXYSpacing()           {return _dxy;}

      void SetZSpacing(double dz) { _dz = dz;}
    double GetZSpacing()          {return _dz;}

    void SetExposureTime(int exptime) {_exposureTime = exptime;}
     int GetExposureTime()            {   return _exposureTime;}

    void SetDetectorGain(int detctgain) {_detectorGain = detctgain;}
     int GetDetectorGain()              {     return _detectorGain;}

    void SetDetectorOffSet(int detctoff) {_detectorOffset = detctoff;}
     int GetDetectorOffSet()             {   return _detectorOffset;}

    void SetBackgroundPhotons(int backphotons) {_backgroundPhotons = backphotons;}
     int GetBackgroundPhotons()                {       return _backgroundPhotons;}

      void SetPSFStdev(double sigma)            { _psfSigma = sigma;}
    double GetPSFStdev()             {  return _psfSigma;}

      void SetPSFRadii(double radii) { _psfRadii = radii;}
    double GetPSFRadii()             {  return _psfRadii;}

    int GetMaximumPhotonEmission() {return (int) ((_snr-1.0) * _backgroundPhotons);}

    void SetPSFFileName(std::string filename) {
        _external_psf = TRUE;
        PSFFileName = filename;
    }
    std::string GetPSFFileName() { return PSFFileName; }

    bool ExternalPSF() {return _external_psf;}

};

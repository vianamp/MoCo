#include <list>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
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

int PoissonGen(double mu);

int ScanFolderForThisExtension(std::string RootFolder, const char ext[], std::vector<std::string> *List);

int SaveImageData(vtkImageData *Image, const char FileName[]);

int Voxelization(std::string _skell_path_prefix, _MoCoControl MoCoControl);

class _MoCoControl {

    double _snr;
    double _dxy, _dz;
    int _exposureTime;
    int _detectorGain;
    int _detectorOffset;
    int _backgroundPhotons;

    double _radii;
    bool _external_psf;
    std::string PSFFileName;
    double _psfSigma, _psfRadii;

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
        SetRadii(0.15);

        _external_psf = FALSE;
    }
    ~_MoCoControl() {

    }
    
      void SetSignalToNoiseRatio(double snr) { _snr = snr;}
    double GetSignalToNoiseRatio()           {return _snr;}

      void SetRadii(double radii) { _radii = radii;}
    double GetRadii()             {  return _radii;}

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

    void DumpVars() {
        printf("      Spacing = %1.3f,%1.3f\n",_dxy,_dz);
        printf("          SNR = %1.3f\n",_snr);
        printf("Exposure Time = %d\n",_exposureTime);
        printf("     Detector = %d, %d\n",_detectorGain,_detectorOffset);
        printf("      Photons = %d\n",_backgroundPhotons);
        printf("        Radii = %1.3f\n",_radii);
        printf(" External PSF = %s\n",(_external_psf)?"TRUE":"FALSE");
    }

};

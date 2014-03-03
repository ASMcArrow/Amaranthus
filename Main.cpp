#include "itkEllipseSpatialObject.h"
#include "itkImage.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkRegionOfInterestImageFilter.h"

#include "vtkSmartPointer.h"
#include "vtkContourFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleImage.h"
#include "vtkPolyDataNormals.h"
#include "vtkDICOMImageReader.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include "vtkColorTransferFunction.h"
#include "vtkPiecewiseFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkAxesActor.h"
#include "vtkMatrix4x4.h"
#include "vtkInteractorStyleTrackballActor.h"
#include "vtkPropPicker.h"
#include "vtkObjectFactory.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"

#define VTKIS_MOVE   10

class MouseInteractorStyle : public vtkInteractorStyleTrackballActor
{
  public:
    static MouseInteractorStyle* New();
    vtkTypeMacro(MouseInteractorStyle, vtkInteractorStyleTrackballActor)

    void OnMouseMove()
    {
        int x = this->Interactor->GetEventPosition()[0];
        int y = this->Interactor->GetEventPosition()[1];

        switch (this->State)
        {
        case VTKIS_ROTATE:
          this->FindPokedRenderer(x, y);
          this->Rotate();
          this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
          break;

        case VTKIS_PAN:
          this->FindPokedRenderer(x, y);
          this->Pan();
          this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
          break;

        case VTKIS_DOLLY:
          this->FindPokedRenderer(x, y);
          this->Dolly();
          this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
          break;

        case VTKIS_SPIN:
          this->FindPokedRenderer(x, y);
          this->Spin();
          this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
          break;

        case VTKIS_USCALE:
          this->FindPokedRenderer(x, y);
          this->UniformScale();
          this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
          break;

         case VTKIS_MOVE:

            int* clickPos = this->GetInteractor()->GetEventPosition();

                  // Pick from this location.
                  vtkSmartPointer<vtkPropPicker>  picker = vtkSmartPointer<vtkPropPicker>::New();
                  picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

                  double* pos = picker->GetPickPosition();
                  std::cout << "Pick position (world coordinates) is: "
                            << pos[0] << " " << pos[1]
                            << " " << pos[2] << std::endl;

                 // if ( picker->GetActor() == NULL) cout << " nothing picked\n";
                  if((pos[0]!=0)&&(pos[1]!=0)&&(pos[2]!=0))
                  {
                  this->InteractionProp->SetPosition(pos);
                  std::cout << "Actor position (world coordinates) is: "
                            << pos[0] << " " << pos[1]
                            << " " << pos[2] << std::endl;
                  }



                  //std::cout << "Picked actor: " << picker->GetActor() << std::endl;
           // std::cout << "Move!" << std::endl;
            break;
        }
    }

    //----------------------------------------------------------------------------

    void EndMove()
    {
        if (this->State != VTKIS_MOVE)
           {
           return;
           }
         this->StopState();
    }

    void OnRightButtonDown()
    {
        int x = this->Interactor->GetEventPosition()[0];
        int y = this->Interactor->GetEventPosition()[1];

        this->FindPokedRenderer(x, y);
        this->FindPickedActor(x, y);
       // if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
        if (this->InteractionProp == NULL)
        {
            return;
        }

        this->GrabFocus(this->EventCallbackCommand);
        this->StartMove();
    }

    void StartMove()
    {
        if (this->State != VTKIS_NONE)
        {
            return;
        }
        this->StartState(VTKIS_MOVE);
    }

    //----------------------------------------------------------------------------
    void OnRightButtonUp()
    {
        switch (this->State)
          {
          case VTKIS_MOVE:
            this->EndMove();
            break;
          }

        if ( this->Interactor )
        {
        this->ReleaseFocus();
        }
    }
};

vtkStandardNewMacro(MouseInteractorStyle)

int main(int argc, char *argv[])
{
    typedef itk::Image<unsigned char,3> ImageType;

    typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
    ConnectorType::Pointer connector = ConnectorType::New();

    typedef itk::ImageSeriesReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    typedef itk::GDCMImageIO  ImageIOType;
    ImageIOType::Pointer dicomIO = ImageIOType::New();
    reader->SetImageIO(dicomIO);

    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
    nameGenerator->SetDirectory("C:/Users/Nastya/Documents/Ocular grant/ITKHead");

    typedef std::vector<std::string>  SeriesIdContainer;
    const SeriesIdContainer &seriesUID = nameGenerator->GetSeriesUIDs();

    std::string seriesIdentifier = seriesUID.begin()->c_str();

    typedef std::vector< std::string > FileNamesContainer;
    FileNamesContainer fileNames = nameGenerator->GetFileNames(seriesIdentifier);
    reader->SetFileNames(fileNames);
    reader->Update();

    ImageType::Pointer DICOMImage = ImageType::New();
    DICOMImage = reader->GetOutput();

    typedef itk::RegionOfInterestImageFilter<ImageType, ImageType> RegionFilterType;
    RegionFilterType::Pointer regionFilter = RegionFilterType::New();

    ImageType::IndexType start;
    start[0] = 200;
    start[1] = 0;
    start[2] = 100;

    ImageType::SizeType  sizeRegion;
    sizeRegion[0] = 40;
    sizeRegion[1] = 190;
    sizeRegion[2] = 50;

    ImageType::RegionType desiredRegion;
    desiredRegion.SetSize(sizeRegion);
    desiredRegion.SetIndex(start);

    regionFilter->SetRegionOfInterest(desiredRegion);
   // regionFilter->SetRegionOfInterest(DICOMImage->GetLargestPossibleRegion());
    regionFilter->SetInput(reader->GetOutput());

    connector->SetInput(regionFilter->GetOutput());
    connector->Update();

    vtkImageData *vtkImage = connector->GetOutput();

    //Read DICOM image and ray trace it
   /* vtkSmartPointer<vtkDICOMImageReader> DICOMreader = vtkSmartPointer<vtkDICOMImageReader>::New();
    DICOMreader->SetDirectoryName("/home/nastya/Grant/Ocular-build/ITKHead/");
    DICOMreader->Update();*/

    vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> rayCastMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
    vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    vtkSmartPointer<vtkPiecewiseFunction> spwf = vtkSmartPointer<vtkPiecewiseFunction>::New();
    vtkSmartPointer<vtkPiecewiseFunction> gpwf = vtkSmartPointer<vtkPiecewiseFunction>::New();

    rayCastMapper->SetInputData(connector->GetOutput());

    // Set color transfer curve for the volume
    ctf->AddHSVPoint(0, .67, .07, 1);
    ctf->AddHSVPoint(94, .67, .07, 1);
    ctf->AddHSVPoint(139, 0, 0, 0);
    ctf->AddHSVPoint(160, .28, .047, 1);
    ctf->AddHSVPoint(254, .38, .013, 1);

    //Set the opacity curve for the volume
    spwf->AddPoint(0, 0.0);
    spwf->AddPoint(80, 0.00);
    spwf->AddPoint(151, 0.01);
    spwf->AddPoint(240,0.1);
    spwf->AddPoint(255,0.1);

    //Set the gradient curve for the volume
    gpwf->AddPoint(0, .2);
    gpwf->AddPoint(10, .2);
    gpwf->AddPoint(25, 1);

    volume->GetProperty()->SetColor(ctf);
    volume->GetProperty()->SetScalarOpacity(spwf);
    volume->GetProperty()->SetGradientOpacity(gpwf);

    volume->SetMapper(rayCastMapper);

    ImageType::DirectionType d = DICOMImage->GetDirection();
    vtkMatrix4x4 *mat=vtkMatrix4x4::New(); //start with identity matrix
    for (int i=0; i<3; i++)
        for (int k=0; k<3; k++)
            mat->SetElement(i,k, d(i,k));

    //counteract the built-in translation by origin
    ImageType::PointType origin = DICOMImage->GetOrigin();
    volume->SetPosition(-origin[0], -origin[1], -origin[2]);

    //add translation to the user matrix
    for (int i=0; i<3; i++)
        mat->SetElement(i,3, origin[i]);
    volume->SetUserMatrix(mat);

    //Add coordinate system axes, so we have a reference for position and orientation
    vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
    axes->SetTotalLength(250,250,250);
    axes->SetShaftTypeToCylinder();
    axes->SetCylinderRadius(0.01);

    // Create itk::SpatialObject and convert it to vtk::PolyData
    typedef itk::EllipseSpatialObject<3> EllipseType;
    typedef EllipseType::TransformType  TransformType;

    EllipseType::Pointer sclera = EllipseType::New();
    sclera->SetRadius(12.37);
    sclera->GetProperty()->SetName("Sclera");

    TransformType::OffsetType ScleraToWorldOffset;
    ScleraToWorldOffset[0] = 0;
    ScleraToWorldOffset[1] = 0;
    ScleraToWorldOffset[2] = 1;
    sclera->GetObjectToParentTransform()->SetOffset(ScleraToWorldOffset);
    sclera->ComputeObjectToWorldTransform();

    typedef itk::SpatialObjectToImageFilter<itk::EllipseSpatialObject<3>,ImageType> CreateImageFilterType;

    CreateImageFilterType::Pointer createImageFilter = CreateImageFilterType::New();
    createImageFilter->SetOutsideValue(0);
    createImageFilter->SetInsideValue(255);
    createImageFilter->SetInput(sclera);

    ImageType::SizeType size;
    size[0] =  30; //image->GetLargestPossibleRegion().GetSize(0);
    size[1] =  30; //image->GetLargestPossibleRegion().GetSize(1);
    size[2] =  30; //image->GetLargestPossibleRegion().GetSize(1);

    createImageFilter->SetSize(size);
    createImageFilter->SetSpacing(1);

    createImageFilter->Update();
    ImageType::Pointer image = createImageFilter->GetOutput();

    ConnectorType::Pointer connector2 = ConnectorType::New();
    connector2->SetInput(image);
    connector2->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->ComputeNormalsOn();
    contourFilter->SetNumberOfContours(1);
    contourFilter->SetValue(0,255);
    contourFilter->SetInputData(connector2->GetOutput());
    contourFilter->Update();

    vtkSmartPointer<vtkPolyDataNormals> normalsFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalsFilter->SetInputData(contourFilter->GetOutput());
    normalsFilter->SetFeatureAngle(80.0);
    normalsFilter->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(normalsFilter->GetOutput());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->SetPosition(DICOMImage->GetOrigin()[0], DICOMImage->GetOrigin()[1], DICOMImage->GetOrigin()[2]);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    renderer->AddActor(axes);
    renderer->AddVolume(volume);

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
    style->SetDefaultRenderer(renderer);

    renderWindowInteractor->SetInteractorStyle(style);

    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindowInteractor->Initialize();

    renderer->ResetCamera();
    renderWindow->Render();
    renderWindowInteractor->Start();
}

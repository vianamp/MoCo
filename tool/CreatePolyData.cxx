#include "CreatePolyData.h"

double _xy, _dz;
std::string OutPutFile;

void SplitString(std::string s, std::string delimiter, std::vector< std::string > *Array) {
	size_t pos = 0;
	Array -> clear();
	std::string token;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		Array -> push_back(token);
		s.erase(0, pos + delimiter.length());
	}
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

void SavePolyData(vtkPolyData *PolyData, std::string Filename) {

	vtkSmartPointer<vtkPolyDataWriter> W = vtkSmartPointer<vtkPolyDataWriter>::New();
	W -> SetFileName(Filename.c_str());
	W -> SetInputData(PolyData);
	W -> Write();

}

vtkPolyData *GetCircle(std::vector< std::string > Array) {

	double u[3];
	double x = atof((Array[0]).c_str());
	double y = atof((Array[1]).c_str());
	double z = atof((Array[2]).c_str());
	double r = atof((Array[3]).c_str());
	double w = atof((Array[4]).c_str());

	printf("Width = %1.3f\n",w);

	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();

	for (double t = 0; t < 2*_PI; t += 0.01) {
		u[0] = x;
		u[1] = y + r * cos(t);
		u[2] = z + r * sin(t);
		Points -> InsertNextPoint(u);
	}

	vtkIdType N = Points -> GetNumberOfPoints();
	vtkSmartPointer<vtkCellArray> Edges = vtkSmartPointer<vtkCellArray>::New();
	for (vtkIdType i = 0; i < N; i++) {
		Edges -> InsertNextCell(2);
		Edges -> InsertCellPoint(i+0);
		Edges -> InsertCellPoint((i+1)<N?i+1:0);
	}

	vtkSmartPointer<vtkFloatArray> Width = vtkSmartPointer<vtkFloatArray>::New();
	Width -> SetName("Width");
	Width -> SetNumberOfTuples(N);
	Width -> SetNumberOfComponents(1);
	Width -> FillComponent(0,w);

	vtkPolyData *PolyData = vtkPolyData::New();
	PolyData -> SetPoints(Points);
	PolyData -> SetLines(Edges);
	PolyData -> GetPointData() -> SetScalars(Width);
	PolyData -> Modified();

	return PolyData;
}

vtkPolyData *GetPoint(std::vector< std::string > Array) {

	double u[3];
	u[0] = atof((Array[0]).c_str());
	u[1] = atof((Array[1]).c_str());
	u[2] = atof((Array[2]).c_str());

	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();

	Points -> InsertNextPoint(u);

	vtkPolyData *PolyData = vtkPolyData::New();
	PolyData -> SetPoints(Points);
	PolyData -> Modified();

	return PolyData;
}

void GeneratePolyData(std::string Path, std::string FileName) {

	std::string line;

	std::vector< std::string > Input;

	std::ifstream input((FileName+".data").c_str());

	while ( getline( input, line ) ) {

		Input.push_back(line);

	}

	int OutPutFile_Pos, _xy_Pos, _z_Pos;

	for (int i = 0; i < Input.size(); i++) {

		line = Input[i];

		if (*line.begin() == '@') {

			printf("Object: %s\n",line.c_str());

			if (line == "@name") {
				OutPutFile_Pos = i + 1;
			}

			if (line == "@xy") {
				_xy_Pos = i + 1;
			}

			if (line == "@z") {
				_z_Pos = i + 1;
			}

		}

	}

	vtkSmartPointer<vtkAppendPolyData> Append = vtkSmartPointer<vtkAppendPolyData>::New();
	Append -> SetOutputPointsPrecision(vtkAlgorithm::DEFAULT_PRECISION);

	OutPutFile = Input[OutPutFile_Pos];

	float x, y, z, r, w;
	std::vector< std::string > Array;

	for (int i = 0; i < Input.size(); i++) {

		line = Input[i];

		if (*line.begin() == '@') {

			printf("Object: %s\n",line.c_str());

			if (line == "@circle") {

				SplitString(Input[i+1]," ", &Array);

				Append -> AddInputData( GetCircle(Array) );
				Append -> Update();
			}

			if (line == "@point") {

				SplitString(Input[i+1]," ", &Array);

				Append -> AddInputData( GetPoint(Array) );
				Append -> Update();
			}

		}

	}
	
	

	printf("N = %d\n",(int)Append->GetOutput()->GetNumberOfPoints());
	SavePolyData(Append->GetOutput(),Path+OutPutFile);

}

int main(int argc, char *argv[]) {     

	//_MoCoControl MoCoControl;

	std::string Path;
	for (int i = 0; i < argc; i++) {
		if (!strcmp(argv[i],"-path")) {
			Path = std::string(argv[i+1]);
			if (Path.back() != '/')
				Path = Path + '/';
		}
	}

	std::vector<std::string> Files;
	ScanFolderForThisExtension(Path,".data",&Files);

	for (int i = 0; i < Files.size(); i++) {

		GeneratePolyData(Path,Files[i]);

		printf("%s [done]\n",Files[i].c_str());
	}

	printf("Starting MoCo tool...\n");

	return EXIT_SUCCESS;
}

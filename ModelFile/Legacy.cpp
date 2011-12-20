bool ModelFile::ObjRead(const string &filename)
{
	// f v/vt/vn
	//f = face
	//v = vertex
	//vt = vertex texture
	//vn = vertex normal
	cout << "Reading obj file " << filename << "..." << endl;
	ifstream indata;
	indata.open(filename.c_str());
	if(!indata)
	{
		cout << "Could not open file!" << endl;
		return false;
	}
	
	string line;
	while(getline(indata, line))
	{
		if(line.size() > 2)
		{
			if(line[0] == 'v' && line[1] == ' ')//a vertex definition
			{
				stringstream ParsedLine(line);
				string LineType;
			
				ParsedLine >> LineType;//do nothing with this
				
				vector<string> Vertex(3);
				ParsedLine >> Vertex[0] >> Vertex[1] >> Vertex[2];
				
				double x = atof(Vertex[0].c_str());
				double y = atof(Vertex[1].c_str());
				double z = atof(Vertex[2].c_str());
				Points_.push_back(vgl_point_3d<double>(x,y,z));
		
			}
			else if(line[0] == 'v' && line[1] == 'n')//a normal definition
			{
				stringstream ParsedLine(line);
				string LineType;
		
				ParsedLine >> LineType;//do nothing with this
			
				vector<string> Normal(3);
				ParsedLine >> Normal[0] >> Normal[1] >> Normal[2];
		
				double x = atof(Normal[0].c_str());
				double y = atof(Normal[1].c_str());
				double z = atof(Normal[2].c_str());
				Normals_.push_back(vgl_vector_3d<double>(x,y,z));    
			}
			else if(line[0] == 'f') //a face definition
			{
				stringstream ParsedLine(line);
				string LineType;
				ParsedLine >> LineType;//do nothing with this
			
				vector<string> VertexSet(3);
				vector<int> VertexIndex(3);
			
				ParsedLine >> VertexSet[0] >> VertexSet[1] >> VertexSet[2];
				
				for(int i = 0; i < 3; i++)
				{
					stringstream VertexParse(VertexSet.at(i));
					string V;
					getline(VertexParse, V, '/');
					VertexIndex[i] = atoi(V.c_str()) - 1; //obj list of indexes starts at 1, not 0
				}
			
				VertexList_.push_back(VertexIndex);
				//cout << "Face " << TriangleVertexList_.size() - 1 << " : " << VertexIndex[0] << " , " << VertexIndex[1] << " , " << VertexIndex[2] << endl;
			}//end else if line == f
		}//end if line size == 2
	
	}//end while getline
    
	OriginalPoints_ = Points_;
	
	Triangles_.clear();
	Triangles_.resize(VertexList_.size());
	
	UpdateTriangles();
	
	cout << "Finished reading obj file " << filename << "." << endl
			<< *this << endl;

}

bool ModelFile::ObjWrite(const string &Filename) const
{
	ofstream fout(Filename.c_str(), ios::out);

	//f v/vt/vn
	//f = face
	//v = vertex
	//vt = vertex texture
	//vn = vertex normal
	
	fout << "# comment" << endl << endl << "mtllib file.mtl" << endl << "g default" << endl;
	
	//write vertices
	for(int i = 0; i < NumVertices(); i++)
	{
		vgl_point_3d<double> P = Points_[i];
		fout << "v " << P.x() << " " << P.y() << " " << P.z() << endl;
	}
			
	//write normals
	//for(int i = 0; i < Normals_.size(); i++)
	//	fout << "vn " << Vertices_[i] << endl;
	
	//write textures
	//for(int i = 0; i < Vertices_.size(); i++)
	//	fout << "vt " << Vertices_[i] << endl;
	
	//write faces
	if(NumTriangles() > 0)
	{
		for(int i = 0; i < VertexList_.size(); i++)
		{
			vector<int> face = VertexList_[i];
			fout << "f " << face[0] + 1 << "/0/0 " << face[1] + 1 << "/0/0 " << face[2] + 1 << "/0/0" << endl; //the first point is "1", not "0"
		}
	}
	
	fout.close();
	
	cout << "Finished writing obj file " << Filename << "." << endl
			<< *this << endl;
	
	return true;//write finshed ok
}
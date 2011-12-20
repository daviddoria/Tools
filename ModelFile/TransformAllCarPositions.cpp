#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "ModelFile.h"

#include <Geometry/Transformation.h>
#include <Tools/Tools.h>

//Neon.vtp NeonTransforms.xform

int main(int argc, char *argv[])
{
	string ModelFilename = argv[1];
	string TransformFile = argv[2];
	
	ModelFile OrigModel(ModelFilename);
	assert(OrigModel.IsValid());
	
	vector<Transformation> Transformations = ReadAllTransformations(TransformFile);

	for(unsigned int i = 0; i < Transformations.size(); i++)
	{
		ModelFile Model = OrigModel;
		Transformation Trans = Transformations[i];
		Model.Transform(Trans);
		
		cout << "Rotation: " << Trans.getRotation() << endl;
		cout << "Translation: " << Trans.getTranslation() << endl;
		
		stringstream OutputFilename;
		OutputFilename << "Transformed_Scene_" << Tools::FileNameWithoutExtension(TransformFile) << "_Position_" << i << ".vtp";
		Model.Write(OutputFilename.str());
	
	}
	
	return 0;
}


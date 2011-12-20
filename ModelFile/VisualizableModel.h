
class VisualizaableModel : ModelFile
{
		//////////// Displays ///////////
		void DrawPoints(void) const;
		void DrawTriangles(void) const;
		void DrawQuads(void) const;
		void LookAt(const vgl_point_3d<double> &CameraPosition) const;
};
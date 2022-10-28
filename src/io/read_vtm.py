from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


def read_vtm(file_name, reader_type=vtkXMLUnstructuredGridReader):
    reader = reader_type()
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutput()
    return output

from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridWriter


def write_vtm(file_name, data, writer_type=vtkXMLUnstructuredGridWriter):
    # alternative vtkXMLMultiBlockDataWriter
    writer = writer_type()
    writer.SetFileName(file_name)
    writer.SetInputData(data)
    writer.Write()

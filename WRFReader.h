/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtPkNetCDFWRFReader.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME WRFReader - read NetCDF files in parallel with MPI
// .Author Ross Miller 03.14.2011
// .SECTION Description
// vtkNetCDFWRFReader is a source object that reads NetCDF files.
// It should be able to read most any NetCDF file that wants to output a
// rectilinear grid.  The ordering of the variables is changed such that
// the NetCDF x, y, z directions correspond to the vtkRectilinearGrid
// z, y, x directions, respectively.  The striding is done with
// respect to the vtkRectilinearGrid ordering.  Additionally, the
// z coordinates of the vtkRectilinearGrid are negated so that the
// first slice/plane has the highest z-value and the last slice/plane
// has the lowest z-value.

#ifndef __WRFReader_h
#define __WRFReader_h

// #include "vtkIOParallelNetCDFModule.h" // For export macro
#include "vtkStructuredGridReader.h"

class vtkDataArraySelection;
class vtkCallbackCommand;
class vtkMPIController;
class WRFReaderInternal;

class WRFReader : public vtkStructuredGridReader
{
public:
  vtkTypeMacro(WRFReader,vtkStructuredGridReader)
  static WRFReader *New();
  void PrintSelf(ostream& os, vtkIndent indent);

  //Description:
  //The file to open
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Variable array selection.
  virtual int GetNumberOfVariableArrays();
  virtual const char *GetVariableArrayName(int idx);
  virtual int GetVariableArrayStatus(const char *name);
  virtual void SetVariableArrayStatus(const char *name, int status);

  // Set/Get the vtkMultiProcessController which will handle communications
  // for the parallel rendering.
  vtkGetObjectMacro(Controller, vtkMPIController);
  void SetController(vtkMPIController *controller);

protected:
  WRFReader();
  ~WRFReader();

  int RequestData(vtkInformation*,vtkInformationVector**,
                  vtkInformationVector*);
  virtual int RequestInformation(vtkInformation* request,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector);

  static void SelectionModifiedCallback(vtkObject *caller, unsigned long eid,
                                        void *clientdata, void *calldata);

  static void EventCallback(vtkObject* caller, unsigned long eid,
                            void* clientdata, void* calldata);

  vtkCallbackCommand* SelectionObserver;

  char *FileName;
  char *OpenedFileName;
  vtkSetStringMacro(OpenedFileName);

  int NCDFFD; //netcdf file descriptor

  vtkMPIController *Controller;

private:
  WRFReader(const WRFReader&);  // Not implemented.
  void operator=(const WRFReader&);  // Not implemented.

  WRFReaderInternal* Internals;
};
#endif

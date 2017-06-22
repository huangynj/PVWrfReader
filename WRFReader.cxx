/*=========================================================================

Program:   Visualization Toolkit
Module:    WRFReader.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include <iostream>
#include <stdio.h>

#include "WRFReader.h"
#include "vtkCallbackCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkMPI.h"
#include "vtkMPIController.h"

#include "netcdf.h"

#include <gdal.h>
#include <ogr_spatialref.h>

#include <string>
#include <vector>
#include <set>

#ifdef MPI_Comm
    #error MPI_Comm is #define'd somewhere!  That's BAD!  (Try checking netcdf.h.)
#endif

void
dbg()
{
#if 0
	fprintf(stderr, ".");
#else
	fprintf(stderr, "debug %d\n", getpid());
  int dbg = 1;
  while (dbg) 
		sleep(1);
#endif
}


vtkStandardNewMacro(WRFReader);

//============================================================================
#define CALL_NETCDF(call) \
{ \
  int errorcode = call; \
  if (errorcode != NC_NOERR) \
  { \
    vtkErrorMacro(<< "netCDF Error: " << nc_strerror(errorcode)); \
    return 0; \
  } \
}
//============================================================================

struct variable
{	
	char name[NC_MAX_NAME+1];
	int  mapped_index;		// which exposed variable corresponds to this netcdf index
};

class WRFReaderInternal
{
public:
	struct variable *Variables;
  vtkSmartPointer<vtkDataArraySelection> VariableArraySelection;
	int number_of_variables;
	int number_of_exposed_variables;
  int dimensions[3];

  //////////////////////

  WRFReaderInternal()
    {
		this->VariableArraySelection = vtkSmartPointer<vtkDataArraySelection>::New();
		this->Variables = NULL;
    }

	~WRFReaderInternal()
	  {
		if (this->Variables == NULL) delete[] this->Variables;
		}
};

//----------------------------------------------------------------------------
//set default values
WRFReader::WRFReader()
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  this->FileName = NULL;
  this->OpenedFileName = NULL;
  this->SelectionObserver = vtkCallbackCommand::New();
  this->SelectionObserver->SetCallback(&WRFReader::SelectionModifiedCallback);
  this->SelectionObserver->SetClientData(this);
  this->Internals = new WRFReaderInternal;
  this->Internals->VariableArraySelection->AddObserver(vtkCommand::ModifiedEvent, this->SelectionObserver);
  this->Controller = NULL;
  this->SetController(vtkMPIController::SafeDownCast(vtkMultiProcessController::GetGlobalController()));
  this->NCDFFD = -1;
}

//----------------------------------------------------------------------------
//delete filename and netcdf file descriptor
WRFReader::~WRFReader()
{
  this->SetController(NULL);
  this->SetFileName(0);
  if(this->OpenedFileName)
    {
    nc_close(this->NCDFFD);
    }
  this->SetOpenedFileName(0);
  if(this->SelectionObserver)
    {
    this->SelectionObserver->Delete();
    this->SelectionObserver = NULL;
    }
  if(this->Internals)
    {
    delete this->Internals;
    this->Internals = NULL;
    }
}

//----------------------------------------------------------------------------
void WRFReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(NULL)") << endl;
  os << indent << "OpenedFileName: "
     << (this->OpenedFileName ? this->OpenedFileName : "(NULL)") << endl;
  if(this->Controller)
    {
    os << indent << "Controller: " << this->Controller << endl;
    }
  else
    {
    os << indent << "Controller: (NULL)" << endl;
    }
  os << indent << "NCDFFD: " << this->NCDFFD << endl;

  this->Internals->VariableArraySelection->PrintSelf(os, indent.GetNextIndent());
}

//----------------------------------------------------------------------------
// RequestInformation supplies global meta information
// This should return the reality of what the reader is going to supply.
// This retrieve the extents for the rectilinear grid
// NC_MAX_VAR_DIMS comes from the nc library
int WRFReader::RequestInformation(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector)
{
	// fprintf(stderr, "RequestInfo\n");
	// dbg();

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  if(this->FileName == NULL)
    {
    vtkErrorMacro("FileName not set.");
    return 0;
    }

	if(this->OpenedFileName)
		{
		nc_close(this->NCDFFD);
		}
	int retval = nc_open(this->FileName, NC_NOWRITE, &this->NCDFFD);//read file
	if (retval != NC_NOERR)//checks if read file error
		{
		vtkErrorMacro(<< "Can't read file " << nc_strerror(retval));
		this->SetOpenedFileName(0);
		return 0;
		}
	this->SetOpenedFileName(this->FileName);
  
	
	int ndims;
  nc_inq_ndims(this->NCDFFD, &ndims);
  for (int i = 0; i < ndims; i++)
	  {
		size_t len;
		char *name = new char[NC_MAX_NAME+1];
		nc_inq_dim(this->NCDFFD, i, name, &len);
		if (! strcmp(name, "west_east"))        this->Internals->dimensions[0] = len;
		else if (! strcmp(name, "south_north")) this->Internals->dimensions[1] = len;
		else if (! strcmp(name, "bottom_top"))  this->Internals->dimensions[2] = len;
		}

  int extent[6];
	extent[0] = extent[2] = extent[4] = 0;
	extent[1] = this->Internals->dimensions[0]-1;
	extent[3] = this->Internals->dimensions[1]-1;
	extent[5] = this->Internals->dimensions[2]-1;

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),extent,6);
  // outInfo->Set(CAN_PRODUCE_SUB_EXTENT(), 1);
	
	nc_inq_nvars(this->NCDFFD, &this->Internals->number_of_variables);
	this->Internals->number_of_variables += 1; // FOR WIND
	this->Internals->Variables = new variable[this->Internals->number_of_variables];
	this->Internals->number_of_exposed_variables = 0;

	// dbg();
	
	// For every variable in the file
	for(int i=0;i< this->Internals->number_of_variables-1; i++)
		{
    int nd;

		CALL_NETCDF(nc_inq_varname(this->NCDFFD, i, this->Internals->Variables[i].name));
		CALL_NETCDF(nc_inq_varndims(this->NCDFFD, i, &nd));

		// Expose the 4-D variables (remember time)  EXCEPT PH and PHB
		if (nd == 4 && strcmp(this->Internals->Variables[i].name, "PH") && strcmp(this->Internals->Variables[i].name, "PHB"))
			{
			this->Internals->Variables[i].mapped_index = this->Internals->number_of_exposed_variables++;
			this->Internals->VariableArraySelection->AddArray(this->Internals->Variables[i].name);
			}
		else
			this->Internals->Variables[i].mapped_index = -1;
		}

	// and expose wind
	strcpy(this->Internals->Variables[this->Internals->number_of_variables-1].name, "wind");
	this->Internals->Variables[this->Internals->number_of_variables-1].mapped_index = -2;	// no matching netcdf variable
	this->Internals->number_of_exposed_variables++;

  return 1;
}

int WRFReader::RequestData(vtkInformation* request,
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector  )
{
  this->UpdateProgress(0);

  int outputPort = request->Get(vtkDemandDrivenPipeline::FROM_OUTPUT_PORT());
  if (outputPort == -1) outputPort = 0;

  // get the data object
  vtkInformation *outInfo = outputVector->GetInformationObject(outputPort);
  vtkDataObject* output = outInfo->Get(vtkDataObject::DATA_OBJECT());

  int subext[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),subext);
  // fprintf(stderr, "extent %d %d %d %d %d %d\n", subext[0], subext[1], subext[2], subext[3], subext[4], subext[5]);

  vtkStructuredGrid *sgrid = vtkStructuredGrid::SafeDownCast(output);
  sgrid->SetExtent(subext);

	size_t start[]= { static_cast<size_t>(subext[4]),			// Z - 0 232
										static_cast<size_t>(subext[0]),			// Y
										static_cast<size_t>(subext[2]) };		// X

  size_t count[]= { static_cast<size_t>(subext[5]-subext[4]+1),
                    static_cast<size_t>(subext[1]-subext[0]+1),
                    static_cast<size_t>(subext[3]-subext[2]+1) };

	// fprintf(stderr, "start: %d %d %d count: %d %d %d\n", start[0], start[1], start[2], count[0], count[1], count[2]);

  int lat_id = -1;
  int lon_id = -1;
  int ph_id  = -1;
  int phb_id = -1;
	for (int i = 0; i < this->Internals->number_of_variables; i++)
		if      (! strcmp(this->Internals->Variables[i].name, "XLAT"))  lat_id = i;
		else if (! strcmp(this->Internals->Variables[i].name, "XLONG")) lon_id = i;
		else if (! strcmp(this->Internals->Variables[i].name, "PH"))    ph_id  = i;
		else if (! strcmp(this->Internals->Variables[i].name, "PHB"))   phb_id = i;

	double *X_buf = new double[count[1] * count[2]];
	double *Y_buf = new double[count[1] * count[2]];

	size_t read_corner_2d[3] = {0,  start[2], start[1]};
	size_t read_count_2d[3]  = {1, count[2], count[1]};

	nc_get_vara_double(this->NCDFFD, lon_id, read_corner_2d, read_count_2d, X_buf);
	nc_get_vara_double(this->NCDFFD, lat_id,  read_corner_2d, read_count_2d, Y_buf);

	float lat[count[1] * count[2]];
	float lon[count[1] * count[2]];
	for (int i = 0; i < count[1] * count[2]; i++) 
	{
		lat[i] = (float)Y_buf[i];
		lon[i] = (float)X_buf[i];
	}

	int utm = (lon[0] - (-180)) / 6;

	char buf[256];
	sprintf(buf, "UTM %d (WGS84)", utm);
	// std::cerr << "using " << buf << "\n";

	OGRSpatialReference oUTM, *oLL;
	oUTM.SetProjCS(buf);
	oUTM.SetWellKnownGeogCS( "WGS84" );
	oUTM.SetUTM(utm, TRUE);
	
	oLL = oUTM.CloneGeogCS();
	OGRCoordinateTransformation *proj = OGRCreateCoordinateTransformation(oLL,  &oUTM);

	// dbg();

	proj->Transform((int)(count[1] * count[2]), X_buf, Y_buf);
	delete proj;

	vtkPoints *points = vtkPoints::New();
	points->SetNumberOfPoints(count[0] * count[1] * count[2]);

	size_t read_corner[4] = {0, start[0], start[2], start[1]};
	size_t read_count[4]  = {1, count[0]+1, count[2], count[1]};	// Going to get the Z-staggered pressure variables

	float *PH_buf  = new float[read_count[1] * read_count[2] * read_count[3]];
	float *PHB_buf = new float[read_count[1] * read_count[2] * read_count[3]];

	nc_get_vara_float(this->NCDFFD, ph_id,  read_corner, read_count, PH_buf);
	nc_get_vara_float(this->NCDFFD, phb_id, read_corner, read_count, PHB_buf);

  int stride = count[2] * count[1];

	float *latlon = new float[2 * count[0] * count[1] * count[2]];

  int next_point = 0;
  int xyz_offset = 0;
	for (int kk = 0; kk < count[0]; kk++)
  {
    int xy_offset = 0;
		for (int jj = 0; jj < count[1]; jj++)
			for (int ii = 0; ii < count[2]; ii++, xyz_offset++, xy_offset ++)
				{
				double x = X_buf[xy_offset];
				double y = Y_buf[xy_offset];
				
				latlon[(next_point<<1) + 0] = lon[xy_offset];
				latlon[(next_point<<1) + 1] = lat[xy_offset];

				double ph  = (PH_buf[xyz_offset] + PH_buf[xyz_offset+stride]) / 2.0;
				double phb = (PHB_buf[xyz_offset] + PHB_buf[xyz_offset+stride]) / 2.0;
				double z = (ph + phb) / 9.81;

				points->SetPoint(next_point++, x, y, z);
				}
	}

	sgrid->SetPoints(points);
	points->Delete();

  vtkFloatArray *latlon_array = vtkFloatArray::New();
	latlon_array->SetNumberOfComponents(2);
	latlon_array->SetName("LatLon");
	latlon_array->SetArray(latlon, 2 * count[0] * count[1] * count[2], 0, 1);	// Hand buffer to VTK array; it'll handle deletion.   We
	sgrid->GetPointData()->AddArray(latlon_array);
	latlon_array->Delete();

	delete[] X_buf;
	delete[] Y_buf;
	delete[] PH_buf;
	delete[] PHB_buf;

  // int do_wind = this->Internals->VariableArraySelection->GetArraySetting(this->Internals->number_of_exposed_variables-1) != 0;
  int do_wind = 1;

	float *u_buf = NULL; int delete_u_buf = 1;
	float *v_buf = NULL; int delete_v_buf = 1;
	float *w_buf = NULL; int delete_w_buf = 1;

  for(size_t i = 0; i < this->Internals->number_of_variables;i++)
    {
		int  mi   = this->Internals->Variables[i].mapped_index;
		char *name = this->Internals->Variables[i].name;

		// int  selected = (mi >= 0) && (this->Internals->VariableArraySelection->GetArraySetting(mi) != 0);
		// int selected = !strcmp(name, "REFL_10CM") || !strcmp(name, "QRAIN") ||
									 // !strcmp(name, "U") || !strcmp(name, "V") || !strcmp(name, "W") ||
									 // !strcmp(name, "PH") || !strcmp(name, "PHB");

		int selected = !strcmp(name, "REFL_10CM") || !strcmp(name, "QRAIN");

		// fprintf(stderr, "name: %s selected: %d\n", name, selected);

		if (! strcmp(name, "U") && (selected || do_wind))
			{
			size_t read_corner[4] = {0, start[0], start[2], start[1]};
			size_t read_count[4]  = {1, count[0], count[2], count[1]+1};	// U is X-staggered

			int   src_size = count[0] * count[2] * (count[1]+1);
			float *tmp = new float[src_size];

			int   dst_size = count[0] * count[2] * count[1];
			u_buf = new float[src_size];

			nc_get_vara_float(this->NCDFFD, i, read_corner, read_count, tmp);

			int stride = 1;

			float *sptr = tmp;
			float *dptr = u_buf;
			for (int ii = 0; ii < count[0]; ii++)	
				for (int jj = 0; jj < count[2]; jj++)
				{
					for (int kk = 0; kk < count[1]; kk++, sptr++, dptr++) 
						*dptr = (sptr[0] + sptr[stride]) / 2.0;
					sptr += stride;
				}

			delete[] tmp;

			if (selected)
				{
				vtkFloatArray *array = vtkFloatArray::New();
				array->SetNumberOfComponents(1);
				array->SetName("U");

				array->SetArray(u_buf, dst_size, 0, 1);	// Hand buffer to VTK array; it'll handle deletion.   We
			  delete_u_buf = 0;												// can still use it for wind later

				sgrid->GetPointData()->AddArray(array);
				array->Delete();
				}
			}
		else if (! strcmp(name, "V") && (selected || do_wind))
			{
			size_t read_corner[4] = {0, start[0], start[2], start[1]};
			size_t read_count[4]  = {1, count[0], count[2]+1, count[1]};	// V is Y staggered

			int   src_size = count[0] * (count[2]+1) * count[1];
			float *tmp = new float[src_size];

			int   dst_size = count[0] * count[2] * count[1];
			v_buf = new float[src_size];

			nc_get_vara_float(this->NCDFFD, i, read_corner, read_count, tmp);

			int stride = count[1];

			float *sptr = tmp;
			float *dptr = v_buf;
			for (int ii = 0; ii < count[0]; ii++)	
				{
				for (int jj = 0; jj < count[2]; jj++)
					for (int kk = 0; kk < count[1]; kk++, sptr++, dptr++) 
						*dptr = (sptr[0] + sptr[stride]) / 2.0;
				sptr += stride;
				}

			delete[] tmp;

			if (selected)
				{
				vtkFloatArray *array = vtkFloatArray::New();
				array->SetNumberOfComponents(1);
				array->SetName("V");

				array->SetArray(v_buf, dst_size, 0, 1);	// Hand buffer to VTK array; it'll handle deletion.   We
			  delete_v_buf = 0;												// can still use it for wind later

				sgrid->GetPointData()->AddArray(array);
				array->Delete();
				}
			}
		else if (!strcmp(name, "W") && (selected || do_wind))
			{
			size_t read_corner[4] = {0, start[0], start[2], start[1]};
			size_t read_count[4]  = {1, count[0]+1, count[2], count[1]};	// W is Z staggered

			int   src_size = (count[0]+1) * count[2] * count[1];
			float *tmp = new float[src_size];

			int   dst_size = count[0] * count[2] * count[1];
			w_buf = new float[src_size];

			nc_get_vara_float(this->NCDFFD, i, read_corner, read_count, tmp);

			int stride = count[1]*count[2];

			float *sptr = tmp;
			float *dptr = w_buf;
			for (int ii = 0; ii < count[0]; ii++)	
				for (int jj = 0; jj < count[2]; jj++)
					for (int kk = 0; kk < count[1]; kk++, sptr++, dptr++) 
						*dptr = (sptr[0] + sptr[stride]) / 2.0;

			delete[] tmp;

			if (selected)
				{
				vtkFloatArray *array = vtkFloatArray::New();
				array->SetNumberOfComponents(1);
				array->SetName(name);

				array->SetArray(w_buf, dst_size, 0, 1);	// Hand buffer to VTK array; it'll handle deletion.   We
			  delete_w_buf = 0;												// can still use it for wind later

				sgrid->GetPointData()->AddArray(array);
				array->Delete();
				}
			}
		else if (strcmp(name, "wind") && selected)
			{
			size_t read_corner[4] = {0, start[0], start[2], start[1]};
			size_t read_count[4]  = {1, count[0], count[2], count[1]};	// NO stagger

			int   size = count[0] * count[2] * count[1];
			float *buf = new float[size];

			nc_get_vara_float(this->NCDFFD, i, read_corner, read_count, buf);

			vtkFloatArray *array = vtkFloatArray::New();
			array->SetNumberOfComponents(1);
			array->SetName(name);
			array->SetArray(buf, size, 0, 1);	// Hand buffer to VTK array; it'll handle deletion
			sgrid->GetPointData()->AddArray(array);
			array->Delete();
			}
		}

	if (do_wind)
		{
		int knt = count[0] * count[1] * count[2];
		float *dst = new float[3*knt];
		float *dptr = dst;
		float *uptr = u_buf;
		float *vptr = v_buf;
		float *wptr = w_buf;

		for (int i = 0; i < knt; i++)
			{
			*dptr++ = *uptr++;
			*dptr++ = *vptr++;
			*dptr++ = *wptr++;
			}

		vtkFloatArray *array = vtkFloatArray::New();
		array->SetNumberOfComponents(3);
		array->SetName("wind");
		array->SetArray(dst, 3*knt, 0, 1);	// Hand buffer to VTK array; it'll handle deletion
		sgrid->GetPointData()->SetVectors(array);
		array->Delete();
		}
			
	if (delete_u_buf && u_buf) delete[] u_buf;
	if (delete_v_buf && v_buf) delete[] v_buf;
	if (delete_w_buf && w_buf) delete[] w_buf;

  return 1;
}

//----------------------------------------------------------------------------
//following 5 functions are used for paraview user interface
void WRFReader::SelectionModifiedCallback(vtkObject*, unsigned long,
                                                   void* clientdata, void*)
{
  static_cast<WRFReader*>(clientdata)->Modified();
}

//-----------------------------------------------------------------------------
int WRFReader::GetNumberOfVariableArrays()
{
  return this->Internals->VariableArraySelection->GetNumberOfArrays();
}

//-----------------------------------------------------------------------------
const char* WRFReader::GetVariableArrayName(int index)
{
  if(index < 0 || index >= this->GetNumberOfVariableArrays())
    {
    return NULL;
    }
  return this->Internals->VariableArraySelection->GetArrayName(index);
}

//-----------------------------------------------------------------------------
int WRFReader::GetVariableArrayStatus(const char* name)
{
  return this->Internals->VariableArraySelection->ArrayIsEnabled(name);
}

//-----------------------------------------------------------------------------
void WRFReader::SetVariableArrayStatus(const char* name, int status)
{
  vtkDebugMacro("Set cell array \"" << name << "\" status to: " << status);
  if(this->Internals->VariableArraySelection->ArrayExists(name) == 0)
    {
    vtkErrorMacro(<< name << " is not available in the file.");
    return;
    }
  int enabled = this->Internals->VariableArraySelection->ArrayIsEnabled(name);
  if(status != 0 && enabled == 0)
    {
    this->Internals->VariableArraySelection->EnableArray(name);
    this->Modified();
    }
  else if(status == 0 && enabled != 0)
    {
    this->Internals->VariableArraySelection->DisableArray(name);
    this->Modified();
    }
}
//
//----------------------------------------------------------------------------
//
void WRFReader::SetController(vtkMPIController *controller)
{
  if(this->Controller != controller)
    {
    this->Controller = controller;
    }
}

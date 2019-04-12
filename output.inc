#include "globals.h"

#include <H5Cpp.h>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

using namespace H5;

namespace output {

//function to output one BZ-index (f_, ky_, kx_) as a 3D matrix
int output_BZ(const std::complex<double>* data, const int size, const std::string file_name, const std::string data_name)
{

    std::cout << "writing to file : " << file_name << '\n';

    const H5std_string filename(file_name);
    const H5std_string dataname(data_name);

    //# of frequencies to be outputted
    const int Nw_internal = size / (Nx * Ny);
    //write as 3D array
    constexpr int rank = 3;

    try {
        //in case an exception is thorwn don't print data but abort
        Exception::dontPrint();

        //if file exists it will be overwritten - meaning of H5F_ACC_TRUNC
        H5File* file = new H5File(file_name, H5F_ACC_TRUNC);

        //set dimensions
        hsize_t dims[rank];
        dims[0] = Nw_internal;
        dims[1] = Ny;
        dims[2] = Nx;

        //create complex data type
        CompType complex_data_type(sizeof(std::complex<double>(0.0, 0.0)));
        complex_data_type.insertMember("r", 0, PredType::NATIVE_DOUBLE);
        complex_data_type.insertMember("i", sizeof(double), PredType::NATIVE_DOUBLE);

        //create a dataspace - defines layout of data
        DataSpace dataspace(rank, dims);

        //create actual dataset - buffer
        DataSet* data_set = new DataSet(file->createDataSet(data_name, complex_data_type, dataspace));
        data_set->write(data, complex_data_type);

        delete data_set;
        delete file;

        //handle exceptions
    } catch (FileIException error) {

        error.printErrorStack();
        return -1;

    } catch (DataSetIException error) {

        error.printErrorStack();
        return -1;

    } catch (DataSpaceIException error) {

        error.printErrorStack();
        return -1;

    } catch (DataTypeIException error) {

        error.printErrorStack();
        return -1;
    }
    return 0;
}

} //end: namespace output

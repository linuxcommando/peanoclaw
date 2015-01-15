#include <cstring>
#include <cmath>
#include"dem.h"
#include<vector>

// ver 1.0.2

const float DEM::INVALID = -1e38f;

DEM::DEM(void) : m_data(NULL), m_boundary_data(NULL), m_bound_bottom(NULL), m_bound_top(NULL), m_bound_left(NULL), m_bound_right(NULL) {
	m_dimension[0]	= m_dimension[1]	= 0;
	m_LowerLeft[0]	= m_LowerLeft[1]	= m_LowerLeft[2] = INVALID;
	m_UpperRight[0]	= m_UpperRight[1]	= m_UpperRight[2]= INVALID;
	m_boundary_size = 0;
}

DEM::~DEM(void) {
	clear();
}

DEM::DEM(const DEM& other) : m_data(NULL), m_boundary_data(NULL), m_bound_bottom(NULL), m_bound_top(NULL), m_bound_left(NULL), m_bound_right(NULL) {
	m_dimension[0] = m_dimension[1] = 0;
	*this = other;
}

bool DEM::save(const std::string& name) const {
	FILE* fptr = NULL;
#ifdef _MSC_VER
	if (fopen_s(&fptr,name.c_str(),"wb")) return false;
#else 
	fptr = fopen(name.c_str(),"wb");
	if (fptr==NULL) return false;
#endif
	bool bRes = save(fptr);
	fclose(fptr);
	return bRes;
}

bool DEM::exportOBJ(const std::string& name) const {
	FILE* fptr = NULL;
#ifdef _MSC_VER
	if (fopen_s(&fptr,name.c_str(),"wt")) return false;
#else 
	fptr = fopen(name.c_str(),"wt");
	if (fptr==NULL) return false;
#endif
	bool bRes = exportOBJ(fptr);
	fclose(fptr);
	return bRes;
}

bool DEM::exportCSV(const std::string& name, const char separator) const {
	FILE* fptr = NULL;
#ifdef _MSC_VER
	if (fopen_s(&fptr,name.c_str(),"wt")) return false;
#else 
	fptr = fopen(name.c_str(),"wt");
	if (fptr==NULL) return false;
#endif
	bool bRes = exportCSV(fptr,separator);
	fclose(fptr);
	return bRes;
}

bool DEM::exportESRI(const std::string& name) const {
	FILE* fptr = NULL;
#ifdef _MSC_VER
	if (fopen_s(&fptr,name.c_str(),"wt")) return false;
#else 
	fptr = fopen(name.c_str(),"wt");
	if (fptr==NULL) return false;
#endif
	bool bRes = exportESRI(fptr);
	fclose(fptr);
	return bRes;
}

bool DEM::save(FILE* stream) const {
	if (stream==NULL) return false;
	int dims[3];
	dims[0] = int(m_dimension[0]);
	dims[1] = int(m_dimension[1]);
	dims[2] = int(m_boundary_size);
	if (fwrite(dims,sizeof(int)*3,1,stream)!=1)					return false;
	if (fwrite(m_LowerLeft,sizeof(double)*3,1,stream)!=1)		return false;
	if (fwrite(m_UpperRight,sizeof(double)*3,1,stream)!=1)		return false;
	
	const size_t chunk = (8<<20)/sizeof(float);	// 8MB chunks
	size_t size = nPixels();
	float* ptr = m_data;
	while (size>0) {
		size_t toWrite = std::min<size_t>(chunk,size);
		if (fwrite(ptr,sizeof(float)*toWrite,1,stream)!=1)		return false;
		size-=toWrite;
		ptr+=toWrite;
	}

	if (dims[2]>0) {
		float *ptr = m_boundary_data;
		size = nBoundaryPixels();
		while (size>0) {
			size_t toWrite = std::min<size_t>(chunk,size);
			if (fwrite(ptr,sizeof(float)*toWrite,1,stream)!=1)	return false;
			size-=toWrite;
			ptr+=toWrite;
		}
	}
	return true;
}

bool DEM::exportOBJ(FILE* stream) const {
#if 0
	if (!stream) return false;
	if (!isValid()) return false;
	double x,y;
	for (int j=0; j<dimension(1); j++) {
		for (int i=0; i<dimension(0); i++) {
			pixelspace_to_worldspace(x,y,i,j);
			double z = operator()(i,j);
			fprintf(stream,"v %g %g %g\n",x,y,z);
		}
	}
	for (int j=0; j<dimension(1)-1; j++) {
		for (int i=0; i<dimension(0)-1; i++) {
			int a = i+j*dimension(0);
			int b = a+1;
			int c = a+dimension(0);
			int d = c+1;
			if (std::abs(operator[](a)-operator[](d))<std::abs(operator[](b)-operator[](c))) { // use diagonal ad
				fprintf(stream,"f %i %i %i\n",a+1,b+1,d+1);
				fprintf(stream,"f %i %i %i\n",a+1,d+1,c+1);
			} else { // use diagonal bc
				fprintf(stream,"f %i %i %i\n",a+1,b+1,c+1);
				fprintf(stream,"f %i %i %i\n",b+1,d+1,c+1);
			}
		}
	}
	return true;
#else
    return false;
#endif
}

bool DEM::exportCSV(FILE* stream, const char separator) const {
	if (!stream) return false;
	if (!isValid()) return false;
	fprintf(stream,"DEM east-major\n");
	fprintf(stream,"width,%i,height,%i\n",dimension(0),dimension(1));
	fprintf(stream,"minEast,%g,maxEast,%g\n",lower_left(0),upper_right(0));
	fprintf(stream,"minNorth,%g,maxNorth,%g\n\n",lower_left(1),upper_right(1));
	double dx = (upper_right(0)-lower_left(0))/double(dimension(0)-1);
	double dy = (upper_right(1)-lower_left(1))/double(dimension(1)-1);
	fprintf(stream,"scaleEast,%g,scaleNorth,%g\n",dx,dy);
	fprintf(stream,"\n\n\n\n");
	for (int j=0; j<dimension(1); j++) {
		for (int i=0; i<dimension(0); i++) {
			fprintf(stream,"%g%c",operator()(i,j),i+1==dimension(0) ? separator :'\n');
		}
	}
	return false;
}

bool DEM::exportESRI(FILE* stream) const {
	if (!stream) return false;
	if (!isValid()) return false;
	fprintf(stream,"ncols     %i\n",dimension(0));
	fprintf(stream,"nrows     %i\n",dimension(1));
	fprintf(stream,"xllcorner %g\n",lower_left(0));
	fprintf(stream,"yllcorner %g\n",lower_left(1));
	double dx = (upper_right(0)-lower_left(0))/double(dimension(0)-1);
	double dy = (upper_right(1)-lower_left(1))/double(dimension(1)-1);
	fprintf(stream,"cellsize  %g\n",0.5*(dx+dy));
	fprintf(stream,"NODATA_value -9999\n");
	for (int j=dimension(1)-1; j>=0; j--) {
		for (int i=0; i<dimension(0); i++) {
			fprintf(stream,"%g%c",operator()(i,j),i+1==dimension(0) ? '\n' : ' ');
		}
	}
	return true;
}



bool DEM::importESRI(const std::string& name) {
	/*FILE* fptr = NULL;
#ifdef _MSC_VER
	if (fopen_s(&fptr,name.c_str(),"rt")) return false;
#else
	fptr = fopen(name.c_str(),"rt");
	if (fptr==NULL) return false;
#endif
	bool bRes = importESRI(fptr);
	fclose(fptr);
	return bRes;*/

    return false;
}

bool DEM::importESRI(FILE* stream) {
	if (stream==NULL) return false;
	clear();
	int width=-1, height=-1;
	lower_left(0)=0.0;
	lower_left(1)=0.0;
	double cellsize=1.0;
	double NODATA = DEM::INVALID;
	size_t addr = 0;
	std::vector<char> chLine;
	chLine.resize(4096);
	while(!feof(stream)) {
#ifdef _MSC_VER
		int result = 0;
		while(result==0) {
			long long ptr = _ftelli64(stream);
			result = fscanf_s(stream,"%[^\n]\n",&chLine[0],chLine.size());
			if (result==0) {
				_fseeki64(stream,ptr,SEEK_SET);
				chLine.resize(2*chLine.size());
			}
		}
#else
		unsigned char bummer[0]; // this has still to be resolved!!
		//fscanf(stream,"%[^\n]\n",chLine);
#endif
		std::string sLine(&chLine[0]);
		while(!sLine.empty()) {
			std::string token = nextToken(sLine);
			if (token==std::string("ncols")) {
				width = atoi(nextToken(sLine).c_str());
				sLine.clear();
			}
			else if (token==std::string("nrows")) {
				height= atoi(nextToken(sLine).c_str());
				sLine.clear();
			}
			else if (token==std::string("xllcorner")) {
				lower_left(0) = atof(nextToken(sLine).c_str());
				sLine.clear();
			}
			else if (token==std::string("yllcorner")) {
				lower_left(1) = atof(nextToken(sLine).c_str());
				sLine.clear();
			}
			else if (token==std::string("cellsize")) {
				cellsize = atof(nextToken(sLine).c_str());
				sLine.clear();
			}
			else if (token==std::string("NODATA_value") || token==std::string("nodata_value")) {
				NODATA = atof(nextToken(sLine).c_str());
				sLine.clear();
			}
			else {
				if (!isValid()) {
					if (width>0 && height>0) { if (!resize(width,height,0)) return false; }
					else return false;
				} 
				if (addr>=size_t(width)*size_t(height)) {
					printf(">>> WARNING: gridfile contains superfluous data -- skipping!\n");
					printf("    Remaining line: \"%s\"\n",sLine.c_str());
					sLine.clear();
				} else {
					int i = int(addr%size_t(width));
					int j = height-1-int(addr/size_t(height));
					double val = atof(token.c_str());
					operator()(i,j)=float(val==NODATA ? DEM::INVALID : val);
					addr++;
				}
			}
		} // end while line empty
	} // end if !feof
	upper_right(0) = lower_left(0) + double(dimension(0)-1)*cellsize;
	upper_right(1) = lower_left(1) + double(dimension(1)-1)*cellsize;
	minmax();
	return true;

}


bool DEM::load(const std::string& name) {
	FILE* fptr = NULL;
#ifdef _MSC_VER
	if (fopen_s(&fptr,name.c_str(),"rb")) return false;
#else
	fptr = fopen(name.c_str(),"rb");
	if (fptr==NULL) {
        return false;
    }
#endif
	bool bRes = load(fptr);
	fclose(fptr);
	return bRes;
}

bool DEM::load(FILE* stream) {
	if (stream==NULL) return false;
	clear();
	int dims[3];
	if (fread(dims,sizeof(int)*3,1,stream)!=1)					return false;
	m_dimension[0] = dims[0];
	m_dimension[1] = dims[1];
	m_boundary_size= dims[2];
	if (fread(m_LowerLeft,sizeof(double)*3,1,stream)!=1)		return false;
	if (fread(m_UpperRight,sizeof(double)*3,1,stream)!=1)		return false;
	m_data = new float[nPixels()];
	if (!m_data) return false;
	const size_t chunk = (8<<20)/sizeof(float); // 8MB chunks
	size_t size = nPixels();
	float* ptr = m_data;
	while (size>0) {
		size_t toRead = std::min<size_t>(chunk,size);
		if (fread(ptr,sizeof(float)*toRead,1,stream)!=1)		return false;
		size-=toRead;
		ptr+=toRead;
	}
	if (m_boundary_size>0) {
		size_t size = nBoundaryPixels();
		m_boundary_data = new float[size];
		if (!m_boundary_data) return false;
		float* ptr = m_boundary_data;
		while (size>0) {
			size_t toRead = std::min<size_t>(chunk,size);
			if (fread(ptr,sizeof(float)*toRead,1,stream)!=1)	return false;
			size-=toRead;
			ptr+=toRead;
		}
		m_bound_bottom	=  m_boundary_data;
		m_bound_left	= &m_bound_bottom[m_boundary_size*(m_dimension[0]+2*m_boundary_size)];
		m_bound_right	= &m_bound_left[m_boundary_size*m_dimension[1]];
		m_bound_top		= &m_bound_right[m_boundary_size*m_dimension[1]];
	}
	return true;
}

std::string DEM::nextToken(std::string& input) const {
	const std::string blanks(" \t\n");
	size_t p = input.find_first_not_of(blanks);
	if (p==std::string::npos) {
		input=std::string("");
		return input;
	}
	input = input.substr(p,input.size());
	size_t end = input.find_first_of(blanks);
	if (end==std::string::npos) end=input.size();
	std::string s = input.substr(0,end);
	input = input.substr(s.size(),input.size()-s.size());
	return s;
}

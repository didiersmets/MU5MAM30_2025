#include <utility>
#include <stdio.h>

#include "sparse_matrix.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


CSRMatrix::CSRMatrix(){
	rows = 0;
	cols = 0;
	symmetric = false;
	nnz = 0;
	row_start = nullptr;
	col = nullptr;
	data = TArray<double>();
}

CSRMatrix::CSRMatrix(CSRPattern &pattern){
	rows = pattern.rows;
	cols = pattern.cols;
	symmetric = pattern.symmetric;
	nnz = pattern.nnz;
	// memcpy(row_start,pattern.row_start.data,pattern.row_start.size);
	// memcpy(col,pattern.col.data,pattern.col.size);
	row_start = pattern.row_start.data;
	col = pattern.col.data;
	data = std::move(TArray<double>(nnz));
}

CSRMatrix::CSRMatrix(CSRPattern &pattern,double default_val){
	rows = pattern.rows;
	cols = pattern.cols;
	symmetric = pattern.symmetric;
	nnz = pattern.nnz;
	// memcpy(row_start,pattern.row_start.data,pattern.row_start.size);
	// memcpy(col,pattern.col.data,pattern.col.size);
	row_start = pattern.row_start.data;
	col = pattern.col.data;
	data = TArray<double>(nnz,default_val);
}

void CSRMatrix::mvp(const double *__restrict x, double *__restrict y) const{
	for (size_t row = 0; row<rows;row++){
		y[row] = 0;
		for(uint32_t k = row_start[row]; k<row_start[row+1];k++){
			uint32_t collumn = col[k];
			y[row]+= data[k] * x[collumn];
		}
	}
}
void CSRMatrix::add_mvp(const double *__restrict x, double *__restrict y) const{
	for (size_t row = 0; row<rows;row++){
		for(uint32_t k = row_start[row]; k<row_start[row+1];k++){
			uint32_t collumn = col[k];
			y[row]+= data[k] * x[collumn];
		}
	}
}
double CSRMatrix::sum() const{
	double acc = 0;
	for (size_t i=0;i<nnz;i++){
		acc += data[i];
	}
	return acc;
}

uint32_t find_dichotomic(uint32_t val,const uint32_t *buf,uint32_t start,uint32_t stop,bool& found){
	do {
		uint32_t mid_pt = (start + stop)/2;
		uint32_t test_val = buf[mid_pt];
		if (val == test_val){
			found = true;
			return mid_pt;
		}
		if (val < test_val){
			stop = mid_pt;
		}else {
			start = mid_pt+1;
		}
	}while (start < stop);
	found = false;
	return 0;
}

double &CSRMatrix::operator()(uint32_t i, uint32_t j){
	bool found = false;
	uint32_t k = find_dichotomic(j,col,row_start[i],row_start[i+1],found);
	assert(found); // tried to access 0 element in sparse matrix
	return data[k];
}

void CSRMatrix::print(){
	for (uint32_t r=0;r<rows;r++){
		uint32_t non_zero_col = row_start[r];
		for (uint32_t y=0; y<cols; y++){
			if (y==col[non_zero_col]){
				printf("%lf ",data[non_zero_col]);
				non_zero_col++;
			} else {
				printf("%lf ",0.0);
			}
		}
		printf("\n");
	}
}
void CSRMatrix::print_pattern(){
	for (uint32_t r=0;r<rows;r++){
		uint32_t non_zero_col = row_start[r];
		for (uint32_t y=0; y<cols; y++){
			if (non_zero_col<row_start[r+1] && y==col[non_zero_col]){
				printf("* ");
				non_zero_col++;
			} else {
				printf("  ");
			}
		}
		printf("\n");
	}
}

/* Visualisation of non zeros */
void spy(const CSRPattern &P, uint32_t width, const char *fname){
	TArray<uint8_t> img(width*width,0);
	for(size_t r=0;r<P.rows;r++){
		for(size_t i=P.row_start[r];i<P.row_start[r+1];i++){
			size_t c = P.col[i];
			size_t x = r * width/P.rows;
			size_t y = c * width/P.cols;
			img[x*width+y] = 255;
		}
	}
	stbi_write_png(fname, width, width, 1, (void *)img.data, width);
}

void dump(CSRMatrix &M, const char *fname){
	FILE *f = fopen(fname,"w");
	for (uint32_t r=0;r<M.rows;r++){
		uint32_t non_zero_col = M.row_start[r];
		for (uint32_t y=0; y<M.cols; y++){
			if (y==M.col[non_zero_col]){
				fprintf(f,"%lf ",M.data[non_zero_col]);
				non_zero_col++;
			} else {
				fprintf(f,"%lf ",0.0);
			}
		}
		fprintf(f,"\n");
	}
	fclose(f);
}

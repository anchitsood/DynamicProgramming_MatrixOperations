/* ////////////////////////////////////////////////////////////

File Name: ps11.cpp
Copyright (c) 2016 Anchit Sood (sood.anchit@gmail.com).  All rights reserved.


Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
   this list of conditions and the following disclaimer in the documentation 
   and/or other materials provided with the distribution.


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

//////////////////////////////////////////////////////////// */




/*
Author: Anchit Sood
24-787: Engineering Computation
Andrew ID: anchits

This program performs some matrix operations. Specifically, it calculates transpose and inverse of a 4x4 matrix.
*/

#include <stdio.h>


//////////////////////////////////////////////////////////
//MatrixTemplate  template class
//////////////////////////////////////////////////////////

template <const int NR, const int NC> class MatrixTemplate
{
protected:
	double matrix[NC*NR];

public:
	void Set(int row, int col, double value);	//Sets the matrix element at (row,col) to 'value'; indexed at 1
	double Value(int row, int col);	//Returns the matrix element at (row,col); indexed at 1
	void Print(void);	//Prints the current state of the matrix
};

template <const int NR, const int NC> void MatrixTemplate<NR,NC>::Set(int row, int col, double value)
{
	if (row >= 1 && row <= NR && col >= 1 && col <= NC)
	{
		matrix[(row - 1)*NC + (col - 1)] = value;
	}
}

template <const int NR, const int NC> double MatrixTemplate<NR,NC>::Value(int row, int col)
{
	if (row >= 1 && row <= NR && col >= 1 && col <= NC)
	{
		return matrix[(row - 1)*NC + (col - 1)];
	}
	else
	{
		return 0.0;
	}
}

template <const int NR, const int NC> void MatrixTemplate<NR,NC>::Print(void)
{
	for (int rowcount = 1; rowcount <= NR; ++rowcount)
	{
		for (int colcount = 1; colcount <= NC; ++colcount)
		{
			printf("%0.4f ", Value(rowcount, colcount));
		}
		printf("\n");
	}
}

//////////////////////////////////////////////////////////
//end of MatrixTemplate template class
//////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
//Matrix4x4 class
//////////////////////////////////////////////////////////

class Matrix4x4 : public MatrixTemplate <4,4>
{
public:
	void Transpose(void);	//Replaces the class instance 4x4 matrix with its transpose
	MatrixTemplate<3,3> MinorMatrix(int row, int col);	//Returns the 3x3 minor matrix linked to the element at (row,col)
	double Determinant3x3(MatrixTemplate <3,3> Matrix3x3);	//Returns the determinant of a 3x3 matrix
	double Determinant4x4(void);	//Returns the determinant of the class instance 4x4 matrix
	MatrixTemplate<4,4> CofactorMatrix(void);	//Returns the cofactor matrix of the class instance 4x4 matrix
	void Invert(void);	//Replaces the class instance 4x4 matrix with its inverse
};

void Matrix4x4::Transpose(void)
{
	for (int rowcount = 1; rowcount <= 4; ++rowcount)
	{
		for (int colcount = rowcount + 1; colcount <= 4; ++colcount)
		{
			double temp = Value(rowcount, colcount);
			//double temp2 = Value(colcount, rowcount);
			Set(rowcount, colcount, Value(colcount, rowcount));
			Set(colcount, rowcount, temp);
		}
	}
}

MatrixTemplate<3,3> Matrix4x4::MinorMatrix(int row, int col)
{
	MatrixTemplate<3,3> SignedMinorMat;

	if (row <= 4 && row >= 1 && col <= 4 && col >= 1)
	{
		
		int rowoffset = 0, coloffset = 0;
		
		for (int rowcount = 1; rowcount <= 4; ++rowcount)
		{
			for (int colcount = 1; colcount <= 4; ++colcount)
			{
				if (rowcount == row)
				{
					++rowcount;
					rowoffset = 1;
				}
				if (colcount == col)
				{
					++colcount;
					coloffset = 1;
				}
				SignedMinorMat.Set(rowcount - rowoffset, colcount - coloffset, Value(rowcount, colcount));
				//coloffset = 0;
			}
			coloffset = 0;
		}
	}
	else
	{
		for (int rowcount = 1; rowcount <= 3; ++rowcount)
		{
			for (int colcount = 1; colcount <= 3; ++colcount)
			{
				if (rowcount == colcount)
				{
					SignedMinorMat.Set(rowcount, colcount, 1);
				}
				else
				{
					SignedMinorMat.Set(rowcount, colcount, 0);
				}
			}
		}
	}
	return SignedMinorMat;
}

double Matrix4x4::Determinant3x3(MatrixTemplate<3,3> Matrix3x3)
{
	return 
		(
			((Matrix3x3.Value(1,1)) * (((Matrix3x3.Value(2,2))*(Matrix3x3.Value(3,3))) - ((Matrix3x3.Value(2,3))*(Matrix3x3.Value(3,2)))))
		-	((Matrix3x3.Value(1,2)) * (((Matrix3x3.Value(2,1))*(Matrix3x3.Value(3,3))) - ((Matrix3x3.Value(2,3))*(Matrix3x3.Value(3,1)))))
		+	((Matrix3x3.Value(1,3)) * (((Matrix3x3.Value(2,1))*(Matrix3x3.Value(3,2))) - ((Matrix3x3.Value(2,2))*(Matrix3x3.Value(3,1)))))
		);
}

double Matrix4x4::Determinant4x4(void)
{
	return (
			(Value(1, 1)*Value(2, 2)*Value(3, 3)*Value(4, 4)) + (Value(1, 1)*Value(2, 3)*Value(3, 4)*Value(4, 2))
		+	(Value(1, 1)*Value(2, 4)*Value(3, 2)*Value(4, 3)) + (Value(1, 2)*Value(2, 1)*Value(3, 4)*Value(4, 3))
		+	(Value(1, 2)*Value(2, 3)*Value(3, 1)*Value(4, 4)) + (Value(1, 2)*Value(2, 4)*Value(3, 3)*Value(4, 1))
		+	(Value(1, 3)*Value(2, 1)*Value(3, 2)*Value(4, 4)) + (Value(1, 3)*Value(2, 2)*Value(3, 4)*Value(4, 1))
		+	(Value(1, 3)*Value(2, 4)*Value(3, 1)*Value(4, 2)) + (Value(1, 4)*Value(2, 1)*Value(3, 3)*Value(4, 2))
		+	(Value(1, 4)*Value(2, 2)*Value(3, 1)*Value(4, 3)) + (Value(1, 4)*Value(2, 3)*Value(3, 2)*Value(4, 1))
		-	(Value(1, 1)*Value(2, 2)*Value(3, 4)*Value(4, 3)) - (Value(1, 1)*Value(2, 3)*Value(3, 2)*Value(4, 4))
		-	(Value(1, 1)*Value(2, 4)*Value(3, 3)*Value(4, 2)) - (Value(1, 2)*Value(2, 1)*Value(3, 3)*Value(4, 4))
		-	(Value(1, 2)*Value(2, 3)*Value(3, 4)*Value(4, 1)) - (Value(1, 2)*Value(2, 4)*Value(3, 1)*Value(4, 3))
		-	(Value(1, 3)*Value(2, 1)*Value(3, 4)*Value(4, 2)) - (Value(1, 3)*Value(2, 2)*Value(3, 1)*Value(4, 4))
		-	(Value(1, 3)*Value(2, 4)*Value(3, 2)*Value(4, 1)) - (Value(1, 4)*Value(2, 1)*Value(3, 2)*Value(4, 3))
		-	(Value(1, 4)*Value(2, 2)*Value(3, 3)*Value(4, 1)) - (Value(1, 4)*Value(2, 3)*Value(3, 1)*Value(4, 2))
		);
}

MatrixTemplate<4,4> Matrix4x4::CofactorMatrix(void)
{
	Matrix4x4 preCofactorMatrix;
	int signbit = 0;
	
	for (int rowcount = 1; rowcount <= 4; ++rowcount)
	{
		for (int colcount = 1; colcount <= 4; ++colcount)
		{
			if (((rowcount + colcount) % 2) == 0)
			{
				signbit = 1;
			}
			else
			{
				signbit = -1;
			}
			preCofactorMatrix.Set(rowcount, colcount, signbit*Determinant3x3(MinorMatrix(rowcount, colcount)));
		}
	}
	preCofactorMatrix.Transpose();
	return preCofactorMatrix;
}

void Matrix4x4::Invert(void)
{
	MatrixTemplate<4,4> Inverted = CofactorMatrix();
	bool flag = true;
	double det = Determinant4x4();

	if (det == 0)
	{
		flag = false;
	}

	for (int rowcount = 1; rowcount <= 4; ++rowcount)
	{
		for (int colcount = 1; colcount <= 4; ++colcount)
		{
			if (flag == true)
			{
				Set(rowcount, colcount, Inverted.Value(rowcount, colcount)/det);
			}
			else
			{
				Set(rowcount, colcount, 0.0);
			}
		}
	}
}

//////////////////////////////////////////////////////////
//end of Matrix4x4 class
//////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
//main() function
//////////////////////////////////////////////////////////

int main(void)
{
	Matrix4x4 mat;
	const double v[] =
	{
		6.0, 1.0, 4.0, 9.0,
		9.0, 8.0, 6.0, 1.0,
		7.0, 2.0, 9.0, 4.0,
		1.0, 7.0, 5.0, 9.0
	};
	for (int i = 0; i<16; ++i)
	{
		const int r = 1 + i / 4;
		const int c = 1 + i % 4;
		mat.Set(r, c, v[i]);
	}
	mat.Print();
	mat.Transpose();
	printf("\n");
	mat.Print();

	printf("\n");
	mat.Invert();
	mat.Print();
	return 0;

	//printf("\n");
	//printf("%0.2f\n", mat.Determinant4x4());
	//printf("%0.2f\n", mat.CofactorMatrix().Value(3,1));
	
	//printf("\n");
	//mat.CofactorMatrix().Print();
	//printf("\n");
	//mat.Invert();
	//mat.Print();

}

//////////////////////////////////////////////////////////
//end of main() function
//////////////////////////////////////////////////////////
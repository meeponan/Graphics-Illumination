#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SSD_util.h"
#define MAXLINELENGTH 1000
#define MAXFILELEN    256

typedef struct {
	double rgba[4];
	double z;
} HIDDEN;

typedef struct {
	double sPoint[4];
	double ePoint[4];
} Line;

char saved_fname[256];
SCENE thescene;
CAMERA vcamera;

void display(void);
int Render_SSD(SCENE *ascene, CAMERA *acamera);
void vecCross(double firstVec[], double secondVec[], double result[]);
void vecUnitization(double vec[],double result[]);
double vecDotProduct(double firstVec[], double secondVec[]);
void matrixInitial(double matrix[][4]);
void matrixMultiply(double firstMatrix[][4], double secondMatrix[][4], int inverse);
void matrixApply(double matrix[][4], double coordinates[]);
double decision(double sPoint[], double ePoint[], double x, double y);
void triangleNormal(double a[], double b[], double c[], double normal[]);
void Illumination(double normal[],double diffuse[],double specular[],double illuPoint[]);
void triRendering(double v0[],double v1[],double v2[],float c0[],float c1[],float c2[],double d[],double s[],int shading,double mInverse[][4]);
void toScreen(double v1[], double v2[], double v3[]);
void toHomo(double v1[], double v2[], double v3[]);

void drawFloor(double xmin, double xmax, double ymin, double ymax, double nX, double nY, double floorEdge, double matrixFinal[][4],Line floor[]);
void getFinalTransformMatrix(double anglePers, double nearPers, double farPers, double rightOrtho, double topOrtho, double farOrtho, double nearOrtho, double pjType, double screenWidth, double screenHeight, double mCam[][4], double matrixFinal[][4], double mInverse[][4]);
void rotateMatrix(double axis[], double mRotate[][4], double mTransform[][4], double radian);
void setBuffer(COLOR_VERTEX colorVertices[MAXLINELENGTH],double matrixFinal[][4],int nM,double d[],double s[],int shading,double mInverse[][4]);


HIDDEN *buffer;
double illuColor[3];



void display(void)
{
  Render_SSD(&thescene, &vcamera);
}

void init (void)
{
  /* select clearing color  to the specified background  */
  glClearColor(thescene.bcolor.rgba[0], thescene.bcolor.rgba[1], 
	       thescene.bcolor.rgba[2], thescene.bcolor.rgba[3]);
  glClear (GL_COLOR_BUFFER_BIT);
  /* initialize viewing values  */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, thescene.screen_w-1.0, 
	  0.0, thescene.screen_h-1.0, -1.0, 1.0);
  
}

double min(double a, double b)
{
  	double min;
  	if(a<=b)
    		min=a;
  	else
    		min=b;
  	return min;  
}

double max(double a, double b)
{
  	double max;
  	if(a>=b)
    		max=a;
  	else
    		max=b;
  	return max;  
}

void vecCross(double firstVec[], double secondVec[], double result[])
{
	result[0] = firstVec[1] * secondVec[2] - firstVec[2] * secondVec[1];
	result[1] = firstVec[2] * secondVec[0] - firstVec[0] * secondVec[2];
	result[2] = firstVec[0] * secondVec[1] - firstVec[1] * secondVec[0];
}

void vecUnitization(double vec[],double result[])
{
	double vecSqrt = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	result[0] = vec[0]/vecSqrt;
	result[1] = vec[1]/vecSqrt;
	result[2] = vec[2]/vecSqrt;
}

double vecDotProduct(double firstVec[], double secondVec[])
{
	double product;
	product = firstVec[0]*secondVec[0]+firstVec[1]*secondVec[1]+firstVec[2]*secondVec[2];
	if(product < 0)
		product = 0;
	return product; 
}

void matrixMultiply(double firstMatrix[][4], double secondMatrix[][4], int inverse)
{
	int i,j;
	double result[4][4];
	for(i = 0; i < 4; i++){
		for(j = 0; j < 4; j++){
			result[i][j] = firstMatrix[i][0] * secondMatrix[0][j] + firstMatrix[i][1] * secondMatrix[1][j] + firstMatrix[i][2] * secondMatrix[2][j] + firstMatrix[i][3] * secondMatrix[3][j];
		}
	}
	if(inverse == 0){
		for(i = 0; i < 4; i++){
			for(j = 0; j < 4; j++){
				secondMatrix[i][j] = result[i][j];
			}
		}
	}
	else{
		for(i = 0; i < 4; i++){
			for(j = 0; j < 4; j++){
				firstMatrix[i][j] = result[i][j];
			}
		}
	}

}

void matrixApply(double matrix[][4], double coordinates[])
{
	int i;
	double tmp[4];
	for(i = 0; i < 4; i++){
		tmp[i] = matrix[i][0] * coordinates[0] + matrix[i][1] * coordinates[1] + matrix[i][2] * coordinates[2] + matrix[i][3] * coordinates[3];
		
	}
	for(i = 0; i < 4; i++){
		coordinates[i] = tmp[i];
	}

}

void matrixInitial(double matrix[][4])
{
	int i,j;
	for(i = 0; i < 4; i++){
		for(j = 0; j < 4; j++){
			if(i == j){matrix[i][j] = 1;}
			else{matrix[i][j] = 0;}
		}
	}
}

double decision(double sPoint[], double ePoint[], double x, double y)
{
			double result;
			result = (sPoint[1]-ePoint[1])*x + (ePoint[0] - sPoint[0])*y + sPoint[0]*ePoint[1] - ePoint[0]*sPoint[1];
			return result;

}

void triangleNormal(double a[], double b[], double c[], double normal[])
{
	int i;
	double b_a[3],c_a[3];
	for(i = 0; i < 3; i++){
		b_a[i] = b[i] - a[i];
		c_a[i] = c[i] - a[i];
	}
	vecCross(b_a,c_a,normal);
	vecUnitization(normal,normal);
}

void toScreen(double v1[], double v2[], double v3[])
{
	int i=0;
	for(i = 0; i < 2; i++){
		v1[i] = v1[i]/v1[3];
		v2[i] = v2[i]/v2[3];
		v3[i] = v3[i]/v3[3];
	}
}

void toHomo(double v1[], double v2[], double v3[])
{
	int i=0;
	for(i = 0; i < 2; i++){
		v1[i] = v1[i]*v1[3];
		v2[i] = v2[i]*v2[3];
		v3[i] = v3[i]*v3[3];
	}		
}

void Illumination(double normal[],double diffuse[],double specular[],double illuPoint[])
{
	int i,l,d;
	illuColor[0] = illuColor[1] = illuColor[2] = 0;
	/*ambient*/
	for(i = 0; i < 3; i++){
		illuColor[i] += diffuse[i] * thescene.ambient[i];
	}

	for(l = 0; l < thescene.nlights; l++){
		for(d = 0; d < thescene.lights[l].ndirections; d++){
			double direction[3] = {thescene.lights[l].directions[d][0],
														 thescene.lights[l].directions[d][1],
														 thescene.lights[l].directions[d][2]};
			vecUnitization(direction,direction);

			double eye[3] = {vcamera.eye.xyzw[0]-illuPoint[0],vcamera.eye.xyzw[1]-illuPoint[1],vcamera.eye.xyzw[2]-illuPoint[2]};
			vecUnitization(eye,eye);
			double half[3] = {eye[0]+direction[0],eye[1]+direction[1],eye[2]+direction[2]};
			vecUnitization(half,half);
			
			/*diffuse*/
			for(i = 0; i < 3; i++){
				illuColor[i] += diffuse[i] * thescene.lights[l].light[i] * vecDotProduct(normal,direction);

			}
			/*specular*/
			for(i = 0; i < 3; i++){
				illuColor[i] += specular[i] * thescene.lights[l].light[i] * pow(vecDotProduct(normal,half),specular[3]);

			}
		}
	}

	for(i = 0; i < 3; i++){
		if(illuColor[i]>1)
			illuColor[i] = 1;
	}
}

void triRendering(double v0[],double v1[],double v2[],float c0[],float c1[],float c2[],double d[],double s[],int shading,double mInverse[][4])
{
	/*triangle constant*/
	int xmax,xmin,ymax,ymin;
	double l0_12,l1_02,l2_01;
	double x_incr_alpha,x_incr_beta,x_incr_gamma;
	double y_incr_alpha,y_incr_beta,y_incr_gamma;
	double alpha0,beta0,gamma0,flag_12,flag_02,flag_01;

	xmin=min(min(v0[0],v1[0]),v2[0]);
  	xmax=max(max(v0[0],v1[0]),v2[0]);
  	ymin=min(min(v0[1],v1[1]),v2[1]);
  	ymax=max(max(v0[1],v1[1]),v2[1]);
	
	/*const*/
	l0_12 = decision(v1,v2,v0[0],v0[1]);
	l1_02 = decision(v0,v2,v1[0],v1[1]);
	l2_01 = decision(v0,v1,v2[0],v2[1]);

	x_incr_alpha = (v1[1]-v2[1])/l0_12;
	x_incr_beta = (v0[1]-v2[1])/l1_02;
	x_incr_gamma = (v0[1]-v1[1])/l2_01;

	y_incr_alpha = (v2[0]-v1[0])/l0_12;
	y_incr_beta = (v2[0]-v0[0])/l1_02;
	y_incr_gamma = (v1[0]-v0[0])/l2_01;

	alpha0 = decision(v1,v2,xmin,ymin)/l0_12;
	beta0 = decision(v0,v2,xmin,ymin)/l1_02;
	gamma0 = decision(v0,v1,xmin,ymin)/l2_01;

	/* (-1,-1) or (-2,-1) */
	flag_12 = decision(v1,v2,-1,-1);
	flag_02 = decision(v0,v2,-1,-1);
	flag_01 = decision(v0,v1,-1,-1);

	if(flag_12 == 0)
			flag_12 = decision(v1,v2,-2,-1);
	if(flag_02 == 0)
			flag_02 = decision(v0,v2,-2,-1);
	if(flag_01 == 0)
			flag_01 = decision(v0,v1,-2,-1);

	int i,x,y;
	double alpha,beta,gamma;	 
	for (y = ymin; y <= ymax; y++){
		alpha = alpha0;
		beta = beta0;
		gamma = gamma0;
		for (x = xmin; x <= xmax; x++){
			if(alpha >= 0 && beta >= 0 && gamma >= 0){
				if( (alpha > 0 || l0_12 * flag_12 >0) && (beta > 0 || l1_02 * flag_02 >0) && (gamma > 0 || l2_01 * flag_01 > 0)){				
					double z = alpha * v0[2] + beta * v1[2] + gamma * v2[2];
					if(z < buffer[y*thescene.screen_w+x].z){
						if(shading == 0){
							for(i = 0; i < 3; i++){
								buffer[y*thescene.screen_w+x].rgba[i] = illuColor[i];

							}			
				              buffer[y*thescene.screen_w+x].z = z;
						}
						else if(shading == 1){
							for(i = 0; i < 3; i++){
								buffer[y*thescene.screen_w+x].rgba[i] = alpha * c0[i] + beta * c1[i] + gamma * c2[i];
							}
						 buffer[y*thescene.screen_w+x].z = z;
						}
						else{	
							double normal[3];
							for(i = 0; i < 3; i++){
								normal[i] = alpha * c0[i] + beta * c1[i] + gamma * c2[i];
							}
							vecUnitization(normal,normal);

							double w = alpha * v0[3] + beta * v1[3] + gamma * v2[3];
							double p[4] = {x*w,y*w,z,w};
							
							matrixApply(mInverse,p);
							Illumination(normal,d,s,p);
							for(i = 0; i < 3; i++){
								buffer[y*thescene.screen_w+x].rgba[i] = illuColor[i];

							}		
							buffer[y*thescene.screen_w+x].z = z;
						}
					}
				}	
			}
			alpha += x_incr_alpha;
			beta += x_incr_beta;
			gamma += x_incr_gamma;
			
		}
		alpha0 += y_incr_alpha;
		beta0 += y_incr_beta;
		gamma0 += y_incr_gamma;
	}
}

//draw floor
void drawFloor(double xmin, double xmax, double ymin, double ymax, double nX, double nY, double floorEdge, double matrixFinal[][4],Line floor[])
{
	int i,j;
	for(i = 0; i < nX; i++){
		floor[i].sPoint[0] = xmin;
		floor[i].sPoint[1] = ymin + i * floorEdge;
		floor[i].sPoint[2] = 0;
		floor[i].sPoint[3] = 1;

		floor[i].ePoint[0] = xmax;
		floor[i].ePoint[1] = floor[i].sPoint[1];
		floor[i].ePoint[2] = 0;
		floor[i].ePoint[3] = 1;	
	}
	
	for(i = i,j = 0; i < nX + nY; i++,j++){
		floor[i].sPoint[0] = xmin + j * floorEdge;
		floor[i].sPoint[1] = ymin;
		floor[i].sPoint[2] = 0;
		floor[i].sPoint[3] = 1;

		floor[i].ePoint[0] = floor[i].sPoint[0];
		floor[i].ePoint[1] = ymax;
		floor[i].ePoint[2] = 0;
		floor[i].ePoint[3] = 1;
	}
	for(i = 0; i < nX + nY; i++){
		matrixApply(matrixFinal,floor[i].sPoint);
		matrixApply(matrixFinal,floor[i].ePoint);
    				  		glVertex2d(floor[i].sPoint[0]/floor[i].sPoint[3],floor[i].sPoint[1]/floor[i].sPoint[3]);
    		glVertex2d(floor[i].ePoint[0]/floor[i].ePoint[3],floor[i].ePoint[1]/floor[i].ePoint[3]);
	}   
}

void getFinalTransformMatrix(double anglePers, double nearPers, double farPers, double rightOrtho, double topOrtho, double farOrtho, double nearOrtho, double pjType, double screenWidth, double screenHeight, double mCam[][4], double matrixFinal[][4], double mInverse[][4])
{
	double mPersp[4][4],mOrtho[4][4],mVp[4][4];
	matrixInitial(mPersp);
	matrixInitial(mOrtho);
	matrixInitial(mVp);
	if(pjType == 1){
		double Pi = 3.141592653;
		double radian = (anglePers/(double)180) * Pi;
		double top = tan(radian/(double)2) * (-nearPers);
		double right = (double)screenWidth/(double)screenHeight * top;
		mPersp[0][0] = nearPers/right;
		mPersp[1][1] = nearPers/top;
		mPersp[2][2] = (farPers + nearPers)/(nearPers - farPers);
		mPersp[2][3] = (2 * nearPers * farPers)/(farPers- nearPers);
		mPersp[3][2] = 1;
		mPersp[3][3] = 0;

		double persInverse[4][4];
		double orthoInverse[4][4];
		matrixInitial(persInverse);
		matrixInitial(orthoInverse);

		persInverse[0][0] = (double)1/nearPers;
		persInverse[1][1] = (double)1/nearPers;
		persInverse[2][2] = 0;
		persInverse[2][3] = 1;
		persInverse[3][2] = (double)-1/(nearPers*farPers);
		persInverse[3][3] = (nearPers+farPers)/(nearPers*farPers);
		
		orthoInverse[0][0] = right;
		orthoInverse[1][1] = top;
		orthoInverse[2][2] = (nearPers-farPers)/(double)2;
		orthoInverse[2][3] = (nearPers+farPers)/(double)2;

	 	matrixMultiply(mInverse,persInverse,1);
	 	matrixMultiply(mInverse,orthoInverse,1);

	}
	else{
		mOrtho[0][0] = (double)1/rightOrtho; //r=-l
		mOrtho[1][1] = (double)1/topOrtho;
		mOrtho[2][2] = (double)2/(farOrtho - nearOrtho);
		mOrtho[2][3] = -(farOrtho + nearOrtho)/(farOrtho - nearOrtho);
	}

	/*View point matrix*/
	mVp[0][0] = screenWidth/(double)2;
	mVp[0][3] = (screenWidth - 1)/(double)2;
	mVp[1][1] = screenHeight/(double)2;
	mVp[1][3] = (screenHeight - 1)/(double)2;
	/*Final transform matrix*/
	matrixMultiply(mCam,matrixFinal,0);
	matrixMultiply(mPersp,matrixFinal,0);
	matrixMultiply(mOrtho,matrixFinal,0);
	matrixMultiply(mVp,matrixFinal,0);
}

void rotateMatrix(double axis[], double mRotate[][4], double mTransform[][4], double radian)
{
	int i,j;
	double u[3],v[3],t[3];
	double tmp[4][4];
	vecUnitization(axis,axis);

	t[0] = axis[0];
	t[1] = axis[1]+1;
	t[2] = axis[2];
	vecCross(t,axis,u);
	vecUnitization(u,u);
	vecCross(axis,u,v);
	vecUnitization(v,v);

	matrixInitial(tmp);
	matrixInitial(mRotate);				
	tmp[0][0] = cos(radian);
	tmp[0][1] = -sin(radian);
	tmp[1][0] = sin(radian);
	tmp[1][1] = cos(radian);
	for(i = 0; i < 3; i++){
		mRotate[0][i] = u[i];
		mRotate[1][i] = v[i];
		mRotate[2][i] = axis[i];

	}
	matrixMultiply(tmp,mRotate,0);
	
	matrixInitial(tmp);
	for(i = 0; i < 3; i++){
		tmp[i][0] = u[i];
		tmp[i][1] = v[i];
		tmp[i][2] = axis[i];

	}									
	matrixMultiply(tmp,mRotate,0);
	matrixMultiply(mRotate,mTransform,0);
}


int Render_SSD(SCENE *ascene, CAMERA *acamera)
{
  /* We clear all pixels  */
  glClearColor(ascene->bcolor.rgba[0], ascene->bcolor.rgba[1],
	       ascene->bcolor.rgba[2], ascene->bcolor.rgba[3]);
  glClear (GL_COLOR_BUFFER_BIT);
  	int i,j;
	double matrixFinal[4][4],mTransform[4][4];
	double intermediaM[4][4],mCam[4][4],mInverse[4][4],tmpMatrix[4][4];
	double mTranslate[4][4],mRotate[4][4],mScale[4][4];
	double homo_coordinates[4] = {0,0,0,1};
	//double mInverse[4][4];
	
	matrixInitial(matrixFinal);
	matrixInitial(mTransform);
	matrixInitial(mTranslate);
	matrixInitial(mRotate);
	matrixInitial(mScale);
	
	int screenW = ascene->screen_w;
	int screenH = ascene->screen_h;

	buffer = (HIDDEN *)malloc(sizeof(HIDDEN) * screenW * screenH);


	for(i = 0;i < screenW * screenH;i++)
	{
		buffer[i].rgba[0] = ascene->bcolor.rgba[0];
		buffer[i].rgba[1] = ascene->bcolor.rgba[1];
		buffer[i].rgba[2] = ascene->bcolor.rgba[2];
		buffer[i].z = 9999;
	}

	/* Camera View */
	int ii;
	matrixInitial(mInverse);
	matrixInitial(tmpMatrix);
	matrixInitial(intermediaM);
	matrixInitial(mCam);
	double u[3],v[3],w[3];
	double gaze[3] = {acamera->gaze.xyzw[0],acamera->gaze.xyzw[1],acamera->gaze.xyzw[2]};
	double upVector[3] = {acamera->upVector.xyzw[0],acamera->upVector.xyzw[1],acamera->upVector.xyzw[2]};
	
	vecUnitization(gaze,w);

	for(ii = 0; ii < 3; ii++){
		w[ii] = -w[ii];
	}
	vecCross(upVector,w,u);
	vecUnitization(u,u);
	vecCross(w,u,v);
	ii = 0;
	for(ii = 0; ii < 3; ii++){
		intermediaM[0][ii] = u[ii];
		intermediaM[1][ii] = v[ii];
		intermediaM[2][ii] = w[ii];
		mCam[ii][3] = -acamera->eye.xyzw[ii];
	}	
	matrixMultiply(intermediaM,mCam,0);
	
	/*inverse  matrixMultiply() 0:normal  1:inverse */
	matrixInitial(intermediaM);
	ii = 0;
	for(ii = 0; ii < 3; ii++){
		intermediaM[ii][0] = u[ii];
		intermediaM[ii][1] = v[ii];
		intermediaM[ii][2] = w[ii];
		tmpMatrix[ii][3] = acamera->eye.xyzw[ii];
	}	
	matrixMultiply(mInverse,tmpMatrix,1);
	matrixMultiply(mInverse,intermediaM,1);


	/*Persective(1) and Orthographic(0) Projection matrix*/
	double anglePers, nearPers, farPers, rightOrtho, topOrtho, farOrtho, nearOrtho,pjType;
	double screenWidth,screenHeight;
	screenWidth = ascene->screen_w;
	screenHeight = ascene->screen_h;
	anglePers = ascene->persp.angle;
	nearPers = ascene->persp.near;
	farPers = ascene->persp.far; 
	rightOrtho = ascene->ortho.right;
	topOrtho = ascene->ortho.top;
	farOrtho = ascene->ortho.far; 
	nearOrtho = ascene->ortho.near;
	pjType = ascene->pjType;

	getFinalTransformMatrix(anglePers, nearPers, farPers, rightOrtho, topOrtho, farOrtho, nearOrtho, pjType, screenWidth, screenHeight, mCam, matrixFinal,mInverse);

	double vpInverse[4][4];
	matrixInitial(vpInverse);
	vpInverse[0][0] = (double)2/ascene->screen_w;
	vpInverse[0][3] = (double)(1-ascene->screen_w)/ascene->screen_w;
	vpInverse[1][1] = (double)2/ascene->screen_h;
	vpInverse[1][3] = (double)(1-ascene->screen_h)/ascene->screen_h;
	matrixMultiply(mInverse,vpInverse,1);

	/* Draw floor */
	double xmin,xmax,ymin,ymax,floorEdge;
	int nX,nY;

	xmin = ascene->floor.xmin;
	xmax = ascene->floor.xmax;
	ymin = ascene->floor.ymin;
	ymax = ascene->floor.ymax;
	floorEdge = ascene->floor.size;
	nY = ((xmax-xmin)/floorEdge) + 1;
	nX = ((ymax-ymin)/floorEdge) + 1;
	Line floor[nX + nY];
	
	glLineWidth(2);
	glBegin(GL_LINES);
	glColor3f(ascene->floor.color.rgba[0],ascene->floor.color.rgba[1],ascene->floor.color.rgba[2]);
	drawFloor(xmin,xmax,ymin,ymax,nX,nY,floorEdge,matrixFinal,floor);
	glEnd();

	/* draw axis*/
	glLineWidth(ascene->axis.width);
	glBegin(GL_LINES);
	if(ascene->isAxis == 1){
		double origin[4] = {0,0,0,1};
		double axisX[4] = {ascene->axis.length,0,0,1};
		double axisY[4] = {0,ascene->axis.length,0,1};
		double axisZ[4] = {0,0,ascene->axis.length,1};
		
		matrixApply(matrixFinal,origin);
		matrixApply(matrixFinal,axisX);
		matrixApply(matrixFinal,axisY);
		matrixApply(matrixFinal,axisZ);

		glColor3f(1,0,0);
		glVertex2d(origin[0]/origin[3],origin[1]/origin[3]);
		glVertex2d(axisX[0]/axisX[3],axisX[1]/axisX[3]);
		glColor3f(0,1,0);
		glVertex2d(origin[0]/origin[3],origin[1]/origin[3]);
		glVertex2d(axisY[0]/axisY[3],axisY[1]/axisY[3]);
		glColor3f(0,0,1);
		glVertex2d(origin[0]/origin[3],origin[1]/origin[3]);
		glVertex2d(axisZ[0]/axisZ[3],axisZ[1]/axisZ[3]);
	}
	glEnd();

	/* implement objects*/
	int nT = 0;
	int nR = 0;
	int nS = 0;
	int nM = 0;
	for(i = 0; i < ascene->nidentities; i++){
		matrixInitial(mTransform);
		for(j = 0; j < ascene->identities[i].inStr_num; j++){
			if(ascene->identities[i].instr[j] == TRANSLATE_KEY){
					matrixInitial(mTranslate);
					mTranslate[0][3] = ascene->translate[nT].xyz[0];
					mTranslate[1][3] = ascene->translate[nT].xyz[1];
					mTranslate[2][3] = ascene->translate[nT].xyz[2];
					matrixMultiply(mTranslate,mTransform,0);
					nT++;
			}
			else if(ascene->identities[i].instr[j] == ROTATE_KEY){
					double axis[3];
					axis[0] = ascene->rotate[nR].xyz[0];
					axis[1] = ascene->rotate[nR].xyz[1];
					axis[2] = ascene->rotate[nR].xyz[2];
					double Pi = 3.141592653;
					double radian = (ascene->rotate[nR].angle/(double)180) * Pi;
					rotateMatrix(axis,mRotate,mTransform,radian);
					nR++;
					
			}
			else if(ascene->identities[i].instr[j] == SCALE_KEY){
					matrixInitial(mScale);
					mScale[0][0] = ascene->scale[nS].xyz[0];
					mScale[1][1] = ascene->scale[nS].xyz[1];
					mScale[2][2] = ascene->scale[nS].xyz[2];
					matrixMultiply(mScale,mTransform,0);
					nS++;
			}
			else if(ascene->identities[i].instr[j] == MESH_KEY){			
					double tM[4][4];
					matrixInitial(tM);
					matrixMultiply(mTransform,tM,0);
					matrixMultiply(matrixFinal,tM,0);
					int k,l,d;

					/*apply transform matrix to all vertices in world coordinates*/
					COLOR_VERTEX colorVertices[ascene->mesh[nM].nvertices];
					for(k = 0; k < ascene->mesh[nM].nvertices; k++){
						colorVertices[k] = ascene->mesh[nM].vertices[k];
						matrixApply(mTransform,colorVertices[k].xyzw);
					}

					//flat shading
					if(ascene->mesh[nM].shading == 0){
						for(k = 0; k < ascene->mesh[nM].npolygons; k++){				
								COLOR_VERTEX vertices[3] = {colorVertices[ascene->mesh[nM].polygons[k].num[0]],
																						colorVertices[ascene->mesh[nM].polygons[k].num[1]],
																						colorVertices[ascene->mesh[nM].polygons[k].num[2]]};
								
								double normal[3]; 
								triangleNormal(vertices[0].xyzw,vertices[1].xyzw,vertices[2].xyzw,normal);
								double center[3];
								center[0] = (vertices[0].xyzw[0] + vertices[1].xyzw[0] + vertices[2].xyzw[0])/(double)3;
								center[1] = (vertices[0].xyzw[1] + vertices[1].xyzw[1] + vertices[2].xyzw[1])/(double)3;
								center[2] = (vertices[0].xyzw[2] + vertices[1].xyzw[2] + vertices[2].xyzw[2])/(double)3;
										Illumination(normal,ascene->mesh[nM].diffuse,ascene->mesh[nM].specular,center);	
	int i;
								for(i = 0; i < 3; i++){			matrixApply(matrixFinal,vertices[i].xyzw);
								}				
								toScreen(vertices[0].xyzw,vertices[1].xyzw,vertices[2].xyzw);								triRendering(vertices[0].xyzw,vertices[1].xyzw,vertices[2].xyzw,0,0,0,0,0,0,mInverse);

						}
					}

					//phong shading
					else if (ascene->mesh[nM].shading == 2){
						double vertex_normal[ascene->mesh[nM].nvertices][3];
						for(k = 0; k < ascene->mesh[nM].nvertices; k++){
						double composeNormal[3] = {0,0,0};
							for(l = 0; l < ascene->mesh[nM].npolygons; l++){
								for(d = 0; d < ascene->mesh[nM].polygons[l].nvertices; d++){
									if(ascene->mesh[nM].polygons[l].num[d] == k){
									COLOR_VERTEX vertices[3] = {colorVertices[ascene->mesh[nM].polygons[l].num[0]],
																								colorVertices[ascene->mesh[nM].polygons[l].num[1]],
																								colorVertices[ascene->mesh[nM].polygons[l].num[2]]};
									double normal[3];
										triangleNormal(vertices[0].xyzw,vertices[1].xyzw,vertices[2].xyzw,normal);
													composeNormal[0] += normal[0];
										composeNormal[1] += normal[1];
										composeNormal[2] += normal[2];
									}
								}
							}										vecUnitization(composeNormal,composeNormal);
						int i=0;			
						for(i = 0; i < 3; i++){
							colorVertices[k].rgba[i] = composeNormal[i];
							}
						}
	setBuffer(colorVertices,matrixFinal,nM,ascene->mesh[nM].diffuse,ascene->mesh[nM].specular,2,mInverse);				
					}
					//smooth shading
					else{
					   for(k = 0; k < ascene->mesh[nM].nvertices; k++){
						double composeNormal[3] = {0,0,0};
						for(l = 0; l < ascene->mesh[nM].npolygons; l++){
						   for(d = 0; d < ascene->mesh[nM].polygons[l].nvertices; d++){
																		       if(ascene->mesh[nM].polygons[l].num[d] == k){
							     COLOR_VERTEX vertices[3] = {colorVertices[ascene->mesh[nM].polygons[l].num[0]],
																								colorVertices[ascene->mesh[nM].polygons[l].num[1]],
																								colorVertices[ascene->mesh[nM].polygons[l].num[2]]};
							     double normal[3];
					triangleNormal(vertices[0].xyzw,vertices[1].xyzw,vertices[2].xyzw,normal);	
														composeNormal[0] += normal[0];								composeNormal[1] += normal[1];								composeNormal[2] += normal[2];
							  }
						      }
						  }
	vecUnitization(composeNormal,composeNormal);
	Illumination(composeNormal,ascene->mesh[nM].diffuse,ascene->mesh[nM].specular,colorVertices[k].xyzw);

	colorVertices[k].rgba[0] = illuColor[0];
	colorVertices[k].rgba[1] = illuColor[1];
	colorVertices[k].rgba[2] = illuColor[2];
						}	
	setBuffer(colorVertices,matrixFinal,nM,0,0,1,mInverse);		
					}
					glLineWidth(ascene->mesh[nM].width);
					glBegin(GL_POINTS);
					for(k = 0; k < ascene->screen_h; k++){
						for(l = 0; l < ascene->screen_w; l++){
						   	if(buffer[k*(ascene->screen_w)+l].z < 9999){										glColor3f(buffer[k*(ascene->screen_w) + l].rgba[0],buffer[k*(ascene->screen_w) + l].rgba[1],buffer[k*(ascene->screen_w) + l].rgba[2]);										glVertex2i(l,k);
							}
						} 					
					}
					glEnd();
					nM++;

			}		
		}
	}


  for(i = 0; i < ascene->nlines; i++){
	ascene->lines[i].vertices[0].xyzw[3] = 1;
	ascene->lines[i].vertices[1].xyzw[3] = 1;
	COLOR_VERTEX vertices[2];
	vertices[0] = ascene->lines[i].vertices[0];
	vertices[1] = ascene->lines[i].vertices[1];
	matrixApply(matrixFinal,vertices[0].xyzw);
	matrixApply(matrixFinal,vertices[1].xyzw);

	glLineWidth(ascene->lines[i].width);
	glBegin(GL_LINES);
	glColor3f(vertices[0].rgba[0],vertices[0].rgba[1],vertices[0].rgba[2]);
	glVertex2d(vertices[0].xyzw[0]/vertices[0].xyzw[3],vertices[0].xyzw[1]/vertices[0].xyzw[3]);
	glColor3f(vertices[1].rgba[0],vertices[1].rgba[1],vertices[1].rgba[2]);
	glVertex2d(vertices[1].xyzw[0]/vertices[1].xyzw[3],vertices[1].xyzw[1]/vertices[1].xyzw[3]);
	glEnd();
	}

  free(buffer);
  glFlush ();
  glutSwapBuffers();
  return 0;
}

void setBuffer(COLOR_VERTEX colorVertices[MAXLINELENGTH],double matrixFinal[][4],int nM,double d[],double s[],int shading,double mInverse[][4])
{
	int k;
	for(k = 0; k < thescene.mesh[nM].npolygons; k++){
		COLOR_VERTEX vertices[3] = {colorVertices[thescene.mesh[nM].polygons[k].num[0]],
																							colorVertices[thescene.mesh[nM].polygons[k].num[1]],
																							colorVertices[thescene.mesh[nM].polygons[k].num[2]]};
		int i=0;
		for(i = 0; i < 3; i++){										matrixApply(matrixFinal,vertices[i].xyzw);
		}																			toScreen(vertices[0].xyzw,vertices[1].xyzw,vertices[2].xyzw);
								triRendering(vertices[0].xyzw,vertices[1].xyzw,vertices[2].xyzw,vertices[0].rgba,vertices[1].rgba,vertices[2].rgba,d,s,shading,mInverse);
	
	}
}

int main(int argc, char** argv)
{
  int ii, jj, kk, argc_1;
  char **my_argv;
  char ssd_fname[MAXFILELEN];
  if (argc < 2) {
    printf("%s:%d Usage: %s SSD_file\n", 
	   __FILE__, __LINE__, argv[0]);
    return 0;
  }
  strcpy(ssd_fname, argv[1]);
  strcpy(saved_fname,"graphics_tmp.ppm");
  argc_1 = argc - 1;
  my_argv = (char **)malloc(sizeof(char *) * argc);
  my_argv[0] = argv[0];
  for (ii=2; ii <= argc; ii++) {
    my_argv[ii-1] = argv[ii];
  }
  glutInit(&argc_1, my_argv);
  free(my_argv);
  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
  /* Set the default size and background */
  Read_SSD_Scene(ssd_fname, &thescene,&vcamera, saved_fname);
  glutInitWindowSize (thescene.screen_w, thescene.screen_h);
  glutInitWindowPosition (50, 50);
  glutCreateWindow (argv[0]);
  init ();
  glutDisplayFunc(display);  
 // glutMouseFunc(mouse);
 // glutKeyboardFunc(ssd_keyboard);
  glutMainLoop();
  return 0;   /* ANSI C requires main to return int; it will never be 
		 reached as glutMainLoop() does not return. */
}


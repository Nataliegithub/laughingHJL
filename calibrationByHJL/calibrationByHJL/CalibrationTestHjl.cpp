/*单个相机标定函数：
输入参数：
const char* imageList      IN保存图片名的txt文件
CvMat* object_points       OUT世界坐标系点的矩阵
CvMat* image_points        OUT图像坐标系矩阵
CvMat* intrinsic_matrix    OUT返回的内参数矩阵
CvMat* distortion_coeffs   OUT返回的畸变向量
int n_boards               IN图片的数量
int board_w                IN每张图片x方向角点数量
int board_h                IN每张图片y方向角点数量
CvSize* imgSize            OUT每张图片的像素尺寸
*/
static void SingleCalib(const char* imageList, CvMat* object_points, CvMat* image_points, CvMat* intrinsic_matrix, CvMat* distortion_coeffs,
	int n_boards, int board_w, int board_h, CvSize* imgSize)
{
	//定义文件类  
	FILE* f;
	fopen_s(&f, imageList, "rt");

	int board_n = board_w*board_h;//每张图片中角点总数量  
	CvSize board_sz = cvSize(board_w, board_h);//角点尺寸矩阵  

	CvPoint2D32f* corners = new CvPoint2D32f[board_w*board_h];//定义用于存放每张图片角点数量的一维点数组  
	CvMat* point_counts = cvCreateMat(n_boards, 1, CV_32SC1);//向量，每个元素代表每张图片角点的数量  
	int successes = 0;//找到全部角点的图片数量  
	int step = 0;//用于记录每张图片角点的起始位置  

				 //文件读取不成功:  
	if (!f)
	{
		fprintf(stderr, "can not open file %s\n", imageList);//要读写, 得知道从哪里读, 往哪里写吧?stderr -- 标准错误输出设备  
		return;
	}
	//利用i循环读取文件中的字符，然后用于读取图片  
	for (int i = 0;; i++)
	{
		//读取图片  
		char buf[1024];//存放读取的字符数组  
		int count = 0, result = 0;//count找的的角点数量，result找角点结果标志，全部角点找到非零，否者为0；  

		if (!fgets(buf, sizeof(buf) - 3, f))//提取文件的字符存放到buf  
			break;
		size_t len = strlen(buf);//len为字符数组的长度  
		while (len > 0 && isspace(buf[len - 1]))//int isspace(int c)检查参数c是否为空格字符，也就是判断是否为空格(' ')、定位字符(' \t ')、CR(' \r ')、换行(' \n ')  
			buf[--len] = '\0';                  //、垂直定位字符(' \v ')或翻页(' \f ')的情况。，既在非空白字符的后面以为添加‘\0’  
		if (buf[0] == '#')//开头为'#',结束此次循环  
			continue;
		IplImage* img = cvLoadImage(buf, 0);//读取图片  
		if (!img)
			break;
		//获取图片尺寸  
		*imgSize = cvGetSize(img);
		//提取角点  
		result = cvFindChessboardCorners(img, cvSize(board_w, board_h), corners, &count, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_NORMALIZE_IMAGE);

		if (result)
		{
			//Calibration will suffer without subpixel interpolation  
			//函数 cvFindCornerSubPix 通过迭代来发现具有亚象素精度的角点位置  
			cvFindCornerSubPix(img, corners, count, cvSize(11, 11), cvSize(-1, -1), cvTermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));
			/*TermCriteria迭代算法的终止准则:
			typedef struct CvTermCriteria0d
			{
			int    type;  CV_TERMCRIT_ITER 和CV_TERMCRIT_EPS二值之一，或者二者的组合
			int    max_iter; 最大迭代次数
			double epsilon; 结果的精确性
			}
			一般表示迭代终止的条件，如果为CV_TERMCRIT_ITER，则用最大迭代次数作为终止条件，如果为CV_TERMCRIT_EPS 则用精度作为迭代条件，
			如果为CV_TERMCRIT_ITER+CV_TERMCRIT_EPS则用最大迭代次数或者精度作为迭代条件，看哪个条件先满足。*/

			//开始保存角点的其实位置  
			step = successes*board_n;
			//将角点从数组corners压入矩阵image_points；以及给对应的object_points赋值  
			for (int i = step, j = 0; j < board_n; ++i, ++j)
			{
				CV_MAT_ELEM(*image_points, float, i, 0) = corners[j].x;
				CV_MAT_ELEM(*image_points, float, i, 1) = corners[j].y;
				CV_MAT_ELEM(*object_points, float, i, 0) = (j / board_w);
				CV_MAT_ELEM(*object_points, float, i, 1) = (j % board_w);
				CV_MAT_ELEM(*object_points, float, i, 2) = 0.0f;
			}
			//给对应图片的point_counts赋值  
			CV_MAT_ELEM(*point_counts, int, successes, 0) = board_n;
			successes++;
		}
		//释放该角点图像  
		cvReleaseImage(&img);
	}
	//关闭文件  
	fclose(f);

	//初始化相机内参矩阵  
	CV_MAT_ELEM(*intrinsic_matrix, float, 0, 0) = 1.0f;
	CV_MAT_ELEM(*intrinsic_matrix, float, 1, 1) = 1.0f;
	//标定相机的内参矩阵和畸变系数向量  
	cvCalibrateCamera2(object_points, image_points, point_counts, *imgSize, intrinsic_matrix, distortion_coeffs, NULL, NULL, 0);

}

static void StereoCalib(CvMat* _left_object_points, CvMat* _left_image_points, CvMat* _right_image_points, CvMat* left_intrinsic, CvMat* right_intrinsic,
	CvMat* left_distortion, CvMat* right_distortion, CvMat* _R, CvMat* _T, CvMat* _E, CvMat* _F, int _n_boards, int _board_w, int _board_h, CvSize _imgSize)
{
	int board_n = _board_w*_board_h;//每张图片中角点总数量  

									//初始化Number_perImg  
	CvMat* Number_perImg = cvCreateMat(1, _n_boards, CV_32SC1);
	int* pInt;
	pInt = (int*)(Number_perImg->data.ptr);
	for (int i = 0; i < _n_boards; ++i)
	{
		*pInt = board_n;
		pInt++;
	}
	//Show(用于调试)  
	/*for (int i = 0; i < _n_boards; i++)
	{
	cout << CV_MAT_ELEM(*Number_perImg, int, 0, i) << endl;
	}*/


	//立体标定.计算 _R, _T, _E, _F,，并同时调整Number_perImg, left_intrinsic, D_left,right_intrinsic, D_right  
	cvStereoCalibrate(_left_object_points, _left_image_points, _right_image_points, Number_perImg, left_intrinsic, left_distortion,
		right_intrinsic, right_distortion, _imgSize, _R, _T, _E, _F,
		cvTermCriteria(CV_TERMCRIT_ITER + CV_TERMCRIT_EPS, 100, 1e-5), CV_CALIB_FIX_ASPECT_RATIO + CV_CALIB_ZERO_TANGENT_DIST + CV_CALIB_SAME_FOCAL_LENGTH);
}

/*计算标定误差函数
CvMat* left_image_points     IN左相机图像坐标系点的矩阵
CvMat* right_image_points    IN右相机图像坐标系点的矩阵
CvMat* left_intrinsic        INandOUT左相机的内参矩阵，经立体标定调整后输出
CvMat* right_intrinsic       INandOUT右相机的内参矩阵，经立体标定调整后输出
CvMat* left_distortion       INandOUT左相机的畸变向量，经立体标定调整后输出
CvMat* right_distortion      INandOUT右相机的畸变向量，经立体标定调整后输出
CvMat* _F                    OUT基础矩阵
int n_boards                 IN图片的数量
int board_w                  IN每张图片x方向角点数量
int board_h                  IN每张图片y方向角点数量
*/
double Calib_Quality_Check(CvMat* _left_image_points, CvMat* _right_image_points, CvMat* left_intrinsic, CvMat* right_intrinsic,
	CvMat* left_distortion, CvMat* right_distortion, CvMat* _F, int _n_boards, int _board_w, int _board_h)
{
	int board_n = _board_w*_board_h;//每张图片中角点总数量  
									//整理LeftPoints。因为_left_image_points数据类型是_n_boards*board_n*3*CV_32FC1，left数据类型是1*_n_boards*board_n*CV_32FC2，对参数进行转换  
	CvMat* left = cvCreateMat(1, _n_boards*board_n, CV_32FC2);
	float *p = (float*)cvPtr2D(left, 0, 0);
	float x, y;
	for (int i = 0; i < _n_boards*board_n; ++i)
	{
		x = CV_MAT_ELEM(*_left_image_points, float, i, 0);
		y = CV_MAT_ELEM(*_left_image_points, float, i, 1);
		*p = x;
		p++;
		*p = y;
		p++;
	}
	//整理LeftPoints。因为_left_image_points数据类型是_n_boards*board_n*3*CV_32FC1，left数据类型是1*_n_boards*board_n*CV_32FC2，对参数进行转换  
	CvMat* right = cvCreateMat(1, _n_boards*board_n, CV_32FC2);
	p = (float*)cvPtr2D(right, 0, 0);
	for (int i = 0; i < _n_boards*board_n; ++i)
	{
		x = CV_MAT_ELEM(*_right_image_points, float, i, 0);
		y = CV_MAT_ELEM(*_right_image_points, float, i, 1);
		*p = x;
		p++;
		*p = y;
		p++;
	}
	//调整畸变向量。因为left_distortion数据类型是5*1的向量，D_left数据类型是1*5向量，对参数进行转换  
	CvMat* D_left = cvCreateMat(1, 5, CV_32FC1);
	for (int i = 0; i < 5; ++i)
	{
		CV_MAT_ELEM(*D_left, float, 0, i) = CV_MAT_ELEM(*left_distortion, float, i, 0);
	}
	//调整畸变向量。因为right_distortion数据类型是5*1的向量，D_right数据类型是1*5向量，对参数进行转换  
	CvMat* D_right = cvCreateMat(1, 5, CV_32FC1);
	for (int i = 0; i < 5; ++i)
	{
		CV_MAT_ELEM(*D_right, float, 0, i) = CV_MAT_ELEM(*right_distortion, float, i, 0);
	}
	//定义极线  
	CvMat* Lines1 = cvCreateMat(1, _n_boards*board_n, CV_32FC3);
	CvMat* Lines2 = cvCreateMat(1, _n_boards*board_n, CV_32FC3);
	//校正点的畸变  
	cvUndistortPoints(_left_image_points, _left_image_points, left_intrinsic, D_left, 0, left_intrinsic);
	cvUndistortPoints(_right_image_points, _right_image_points, right_intrinsic, D_right, 0, right_intrinsic);
	//计算极线  
	cvComputeCorrespondEpilines(_left_image_points, 1, _F, Lines1);
	cvComputeCorrespondEpilines(_right_image_points, 2, _F, Lines2);
	//计算左相机图像点和极线的距离  
	float a, b, c;//极线系数ax+by+c=0  
	double err, avgErr = 0;//误差变量  
	float* p_temp = (float*)cvPtr2D(_left_image_points, 0, 0);//用于临时计算的点  
	float* l_temp = (float*)cvPtr2D(Lines2, 0, 0, 0);//用于临时计算的极线  
	for (int i = 0; i < _n_boards*board_n; i++)
	{
		//提取点坐标  
		x = *p_temp;
		y = *(p_temp + 1);
		p_temp += 2;
		//提取极线系数  
		a = *l_temp;
		b = *(l_temp + 1);
		c = *(l_temp + 2);
		l_temp += 3;
		//计算点到直线的距离  
		err = fabs(a*x + b*y + c) / sqrt(a*a + b*b);
		//累加误差  
		avgErr += err;
	}
	//计算右相机图像点和极线的距离  
	p_temp = (float*)cvPtr2D(_right_image_points, 0, 0);//用于临时计算的点  
	l_temp = (float*)cvPtr2D(Lines1, 0, 0);//用于临时计算的极线  
	for (int i = 0; i < _n_boards*board_n; i++)
	{
		//提取点坐标  
		x = *p_temp;
		y = *(p_temp + 1);
		p_temp += 2;
		//提取极线系数  
		a = *l_temp;
		b = *(l_temp + 1);
		c = *(l_temp + 2);
		l_temp += 3;
		//计算点到直线的距离  
		err = fabs(a*x + b*y + c) / sqrt(a*a + b*b);
		//累加误差  
		avgErr += err;
	}
	//求误差的平均值  
	avgErr /= (_n_boards*board_n);
	return avgErr;
}


/*立体校正的映射矩阵mapx和mapy函数
CvMat* left_image_points     IN左相机图像坐标系点的矩阵
CvMat* right_image_points    IN右相机图像坐标系点的矩阵
CvMat* left_intrinsic        INandOUT左相机的内参矩阵，经立体标定调整后输出
CvMat* right_intrinsic       INandOUT右相机的内参矩阵，经立体标定调整后输出
CvMat* left_distortion       INandOUT左相机的畸变向量，经立体标定调整后输出
CvMat* right_distortion      INandOUT右相机的畸变向量，经立体标定调整后输出
CvMat* R                     OUT相机相对旋转矩阵
CvMat* T                     OUT相机相对平移矩阵
CvMat* F                     OUT基础矩阵
CvSize* imgSize              IN每张图片的像素尺寸
int n_boards                 IN图片的数量
int board_w                  IN每张图片x方向角点数量
int board_h                  IN每张图片y方向角点数量
bool useUncalibrated         IN是否采用非标定算法Hartly计算立体校正的映射矩阵mapx和mapy
CvMat* mxl                   OUT左相机立体校正的映射矩阵mapx
CvMat* myl                   OUT左相机立体校正的映射矩阵mapy
CvMat* mxr                   OUT右相机立体校正的映射矩阵mapx
CvMat* myr                   OUT右相机立体校正的映射矩阵mapy
*/
/*mapx和mapy是函数返回的查找映射表，。对目标图像上的每一个像素来说，映射表说明应该从什么位置插值源像素；
并且映射表可以直接被Remap()嵌入使用。函数cvInitUndistortRectifyMap（）被左右摄像机分别调用，这样可以获得它们的各自不同的重映射参数mapx和mapy，
每次我们有新的左右立体图像需要校正时可以使用左右映射表*/
static void ComputeRectification(CvMat* imagePoints_left, CvMat* imagePoints_right, CvMat* _left_intrinsic, CvMat* _right_intrinsic, CvMat* _left_distortion, CvMat* _right_distortion,
	CvMat* R, CvMat* T, CvMat* F, CvSize imageSize, int _n_boards, int _board_w, int _board_h, CvMat* mxl, CvMat* myl, CvMat* mxr, CvMat* myr)
{
	/*将参数变量重新整理成CV_64F类型的，我也不知道为什么，cvStereoRectify如果传入CV_32FC1的参数会中断，改为CV_64FC1则成功*/
	//整理左相机内参矩阵。因为_left_intrinsic数据类型为CV_32FC1，M_left为CV_64FC1  
	CvMat* M_left = cvCreateMat(3, 3, CV_64FC1);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			CV_MAT_ELEM(*M_left, double, i, j) = CV_MAT_ELEM(*_left_intrinsic, float, i, j);
	//整理右相机内参矩阵。因为_left_intrinsic数据类型为CV_32FC1，M_left为CV_64FC1  
	CvMat* M_right = cvCreateMat(3, 3, CV_64FC1);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			CV_MAT_ELEM(*M_right, double, i, j) = CV_MAT_ELEM(*_right_intrinsic, float, i, j);
	//调整畸变向量。因为left_distortion数据类型是5*1的向量，D_left数据类型是1*5向量，同时转换为CV_64FC1。  
	CvMat* D_left = cvCreateMat(1, 5, CV_64FC1);
	for (int i = 0; i < 5; ++i)
	{
		CV_MAT_ELEM(*D_left, double, 0, i) = CV_MAT_ELEM(*_left_distortion, float, i, 0);
	}
	//调整畸变向量。因为right_distortion数据类型是5*1的向量，D_right数据类型是1*5向量，同时转换为CV_64FC1。  
	CvMat* D_right = cvCreateMat(1, 5, CV_64FC1);
	for (int i = 0; i < 5; ++i)
	{
		CV_MAT_ELEM(*D_right, double, 0, i) = CV_MAT_ELEM(*_right_distortion, float, i, 0);
	}
	//转换旋转矩阵。CV_32FC1转换为CV_64FC1  
	CvMat* RR = cvCreateMat(3, 3, CV_64FC1);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			CV_MAT_ELEM(*RR, double, i, j) = CV_MAT_ELEM(*R, float, i, j);
	//定义平移向量。CV_32FC1转换为CV_64FC1  
	CvMat* TT = cvCreateMat(3, 1, CV_64FC1);
	for (int i = 0; i < 3; ++i)
	{
		CV_MAT_ELEM(*TT, double, i, 0) = CV_MAT_ELEM(*T, float, i, 0);
	}
	//定义和初始化左右相机的行对准旋转矩阵  
	CvMat* Rl = cvCreateMat(3, 3, CV_64FC1);
	cvZero(Rl);
	CvMat* Rr = cvCreateMat(3, 3, CV_64FC1);
	cvZero(Rr);
	//计算相机共面和行对准的校正项  
	if (useUncalibrated == false)//useUncalibrated == 0，使用Bouguet方法计算校正项  
	{
		//定义和初始化做右相机投影矩阵，该矩阵可以计算校正后的相机内参数矩阵  
		CvMat* Pl = cvCreateMat(3, 4, CV_64FC1);
		cvZero(Pl);
		CvMat* Pr = cvCreateMat(3, 4, CV_64FC1);
		cvZero(Pr);
		//计算相机校正项，计算得到 Rl, Rr, Pl, Pr  
		cvStereoRectify(M_left, M_right, D_left, D_right, imageSize,
			RR, TT, Rl, Rr, Pl, Pr, 0, 0, -1.0, cvSize(0, 0), (CvRect*)0, (CvRect*)0);
		//Show（用于调试）  
		/*cout << endl << "Pr:" << endl;
		for (int i = 0; i < 3; i++)
		{
		const double *p = (const double*)(Pr->data.ptr + i*Pr->step);
		cout << *p << "    " << *(p + 1) << "       " << *(p + 2)<<"        "<<*(p+3) << endl;
		}*/
		/*通过比较转换矩阵P位置[0][3]和位置[1][3]的大小 判断相机是水平布置还是垂直布置。
		因为这两个位置的元素分别表示右相机坐标系原点移动到左相机坐标系原点的x和y方向移动距离。*/
		isVerticalStereo = fabs(CV_MAT_ELEM(*Pr, double, 1, 3)) > fabs(CV_MAT_ELEM(*Pr, double, 0, 3));
		//Precompute maps for cvRemap()  
		//下面通过函数cvInitUndistortRectifyMap()计算立体校正的映射矩阵mx和my  
		cvInitUndistortRectifyMap(_left_intrinsic, D_left, Rl, Pl, mxl, myl);
		cvInitUndistortRectifyMap(_right_intrinsic, D_right, Rr, Pr, mxr, myr);
	}
	else//useUncalibrated！ = 0，使用Hartley方法计算校正项  
	{
		int board_n = _board_w*_board_h;//每张图片中角点总数量  
		CvMat* Hl = cvCreateMat(3, 3, CV_64FC1);//Hartley方法计算得到的校正项  
		CvMat* Hr = cvCreateMat(3, 3, CV_64FC1);//Hartley方法计算得到的校正项  
		CvMat *invM = cvCreateMat(3, 3, CV_64FC1);//用于临时存放相机内参数的逆  
												  //整理imagePoints_left和imagePoints_right  
		float x, y;//用于传值的临时变量  
		double* p;//用于传值的临时变量  
				  //LeftPoints  
		CvMat* left = cvCreateMat(1, _n_boards*board_n, CV_64FC2);
		p = (double*)cvPtr2D(left, 0, 0);
		for (int i = 0; i < _n_boards*board_n; ++i)
		{
			x = CV_MAT_ELEM(*imagePoints_left, float, i, 0);
			y = CV_MAT_ELEM(*imagePoints_left, float, i, 1);
			*p = x;
			p++;
			*p = y;
			p++;
		}
		//Show  
		/*p = (float*)cvPtr2D(left, 0, 0);
		for (int i = 0; i < n_boards*board_n; i++)
		{
		cout << *p << "    " << *(p + 1) << endl;
		p += 2;
		}*/
		//RightPoints  
		CvMat* right = cvCreateMat(1, _n_boards*board_n, CV_64FC2);
		p = (double*)cvPtr2D(right, 0, 0);
		for (int i = 0; i < _n_boards*board_n; ++i)
		{
			x = CV_MAT_ELEM(*imagePoints_right, float, i, 0);
			y = CV_MAT_ELEM(*imagePoints_right, float, i, 1);
			*p = x;
			p++;
			*p = y;
			p++;
		}
		//Show  
		/*p = (float*)cvPtr2D(right, 0, 0);
		for (int i = 0; i < n_boards*board_n; i++)
		{
		cout << *p << "    " << *(p + 1) << endl;
		p += 2;
		}*/
		//通过两个相机的图像点计算基础矩阵F  
		cvFindFundamentalMat(left, right, F);
		//使用非标定算法Hartley立体校正，获得做右相机的校正单应性矩阵H1,H2  
		cvStereoRectifyUncalibrated(left, right, F, imageSize, Hl, Hr, 3);
		/*使用cvStereoRectifyUncalibrated()的结果H1,H2来校正立体摄像机，必须先预处理一下单应性矩阵Hl和我Hr，
		如果没有先前标定中获得Mrect，令Mrect=M,分别计算左右校正旋转矩阵Rl=inv(Mrect_l)Hl(Ml)和Rr=inv(Mrect_r)Hr(Mr)，
		最后还需要分别计算左右相机的畸变系数向量*/
		//Show（用于调试）  
		/*cout << endl << "M1:";
		for (int i = 0; i < 3; i++)
		{
		cout << endl;
		for (int j = 0; j < 3; j++)
		{
		cout << CV_MAT_ELEM(*_left_intrinsic, float, i, j) << "     ";
		}
		}*/
		//预处理Rl=inv(Mrect_l)Hl(Ml)  
		cvInvert(M_left, invM);
		cvMatMul(Hl, M_left, Rl);
		cvMatMul(invM, Rl, Rl);
		//预处理Rr=inv(Mrect_r)Hr(Mr)  
		cvInvert(M_right, invM);//_iM为_M1的逆矩阵  
		cvMatMul(Hr, M_right, Rr);//_R1=_H1*_M1  
		cvMatMul(invM, Rr, Rr);//_R2=_iM*_R1  
							   //Precompute map for cvRemap()，用预处理后的_R1和_R2，以及前面计算的相机参数代入cvInitUndistortRectifyMap计算mapx和mapy  
		cvInitUndistortRectifyMap(M_left, D_left, Rl, M_left, mxl, myl);
		cvInitUndistortRectifyMap(M_right, D_right, Rr, M_right, mxr, myr);

	}
}

/*校正图片和计算视察映射的函数
const char* imageList_left       IN保存左相机图片名的txt文件
const char* imageList_right      IN保存右相机图片名的txt文件
CvSize* imgSize                  IN每张图片的像素尺寸
CvMat* mxl                       IN左相机立体校正的映射矩阵mapx
CvMat* myl                       IN左相机立体校正的映射矩阵mapy
CvMat* mxr                       IN右相机立体校正的映射矩阵mapx
CvMat* myr                       IN右相机立体校正的映射矩阵mapy
*/
static void Rectify_and_ComputeDisparity(const char* imageList_left, const char* imageList_right, CvSize imSize, CvMat* mxl, CvMat* myl, CvMat* mxr, CvMat* myr)
{
	FILE *f_left, *f_right;//定义两个文件操作对象  
	fopen_s(&f_left, imageList_left, "rt");//打开左相机图像的图像名文件  
	fopen_s(&f_right, imageList_right, "rt");//打开右相机图像的图像名文件  
	IplImage *img_left, *img_right;//定义两个图像对象  
								   //定义两个矩阵用于存放校正后的图像  
	CvMat* img1r = cvCreateMat(imSize.height, imSize.width, CV_8U);
	CvMat* img2r = cvCreateMat(imSize.height, imSize.width, CV_8U);

	CvMat* disp = cvCreateMat(imSize.height, imSize.width, CV_16S);//视察映射  
	CvMat* vdisp = cvCreateMat(imSize.height, imSize.width, CV_8U);//归一化后的视差映射  
																   /*创建和初始化一个结构体CvStereoBMState用于cvFindStereoCorrespondenceBM()的匹配过程的参数。因为如果每次调用cvFindStereoCorrespondenceBM()
																   都离散分配内部离散缓存会影响程序的运行速度。*/
	CvStereoBMState *BMState = cvCreateStereoBMState();
	assert(BMState != 0);//如果创建失败则终止程序  
	BMState->preFilterSize = 41;//预过滤窗口的尺寸，使图像亮度归一化并加强纹理  
	BMState->preFilterCap = 31;//预过滤的阈值  
	BMState->SADWindowSize = 41;//"绝对误差累计"小窗口(SAD)的尺寸  
	BMState->minDisparity = -64;//匹配搜索起始点位置  
	BMState->numberOfDisparities = 128;//匹配搜索 的视差范围  
	BMState->textureThreshold = 10;//后过滤阈值，响应小于该值的匹配不予考虑  
	BMState->uniquenessRatio = 15;//后过滤比例阈值，uniquenessRatio>(match_val-min_match)/min_match的地方过滤掉匹配值  

	cvNamedWindow("disparity");//创建一个窗口用于显示视差映射  
	cvNamedWindow("rectified", 1);//创建窗口，用于显示校正后的图像  
	CvMat* pair;//定义一个矩阵用于存放校正后的两张图片  
	CvMat part;//定义一个矩阵用于存放校正后的一张图片  
	if (!isVerticalStereo)//相机水平布置则创建尺寸为（imageSize.height, imageSize.width * 2）的图片用于显示左右相机立体校正后的图片  
		pair = cvCreateMat(imSize.height, imSize.width * 2, CV_8UC3);
	else//相机竖直布置则创建尺寸为（imageSize.height * 2, imageSize.width）的图片用于显示上下相机立体校正后的图片  
		pair = cvCreateMat(imSize.height * 2, imSize.width, CV_8UC3);

	//循环读取左右相机的图像  
	for (int i = 0;; i++)
	{
		//读取图片  
		char buf_left[1024], buf_right[1024];//存放读取的字符数组  
		if ((!fgets(buf_left, sizeof(buf_left) - 3, f_left)) || (!fgets(buf_right, sizeof(buf_right) - 3, f_right)))//提取文件的字符存放到buf  
			break;
		size_t len = strlen(buf_left);//len为字符数组的长度  
		while (len > 0 && isspace(buf_left[len - 1]))//int isspace(int c)检查参数c是否为空格字符，也就是判断是否为空格(' ')、定位字符(' \t ')、CR(' \r ')、换行(' \n ')  
			buf_left[--len] = '\0';                 //、垂直定位字符(' \v ')或翻页(' \f ')的情况。，既在非空白字符的后面以为添加‘\0’  
		len = strlen(buf_right);//len为字符数组的长度  
		while (len > 0 && isspace(buf_right[len - 1]))//int isspace(int c)检查参数c是否为空格字符，也就是判断是否为空格(' ')、定位字符(' \t ')、CR(' \r ')、换行(' \n ')  
			buf_right[--len] = '\0';                    //、垂直定位字符(' \v ')或翻页(' \f ')的情况。，既在非空白字符的后面以为添加‘\0’  
		if ((buf_left[0] == '#') || (buf_right[0] == '#'))//开头为'#',结束此次循环  
			continue;
		img_left = cvLoadImage(buf_left, 0);//读取图片  
		img_right = cvLoadImage(buf_right, 0);//读取图片  
		if ((!img_left) || (!img_right))
			break;
		//使用用前面计算的mx和my立体校正图片，两个相机的图片处于同一个相机平面且行对准  
		cvRemap(img_left, img1r, mxl, myl);
		cvRemap(img_right, img2r, mxr, myr);
		/*接下来用立体校准后的图片计算双目视觉模型的视察映射。分为相机水平布置和垂直布置两种情况，
		对于摄像机垂直布置的情形，如果没有自己添加图像转置的代码，函数cvFindStereoCorrespondenceBM()
		只能计算非标定校正的视察，对于摄像机水平对齐的情形，函数cvFindStereoCorrespondenceBM()可以计算
		标定和非标定的校正立体像对的视差。*/
		if (!isVerticalStereo || useUncalibrated != 0)
		{
			cvFindStereoCorrespondenceBM(img1r, img2r, disp, BMState);//计算视差映射  
			cvNormalize(disp, vdisp, 0, 256, CV_MINMAX);//归一化视差映射  
														//显示视差映射  
			cvShowImage("disparity", vdisp);
		}
		//下面显示立体校正后的两个相机对应的图片  
		if (!isVerticalStereo)//水平布置相机的情形  
		{
			//pair的左侧复制左相机立体校正后的照片  
			cvGetCols(pair, &part, 0, imSize.width);/*CvMat* cvGetCols( const CvArr* arr, CvMat* submat, int start_col, int end_col );
													返回数组的一定跨度内的列
													arr-输入数组
													submat-指向结果子数组头指针.
													start_col-跨度的开始列（包括该列）索引下标，该下标以0为基准。
													end_col-跨度的结束列（不包括该列）索引下标，该下标以0为基准。*/
			cvCvtColor(img1r, &part, CV_GRAY2BGR);//GRAY->RGB  
												  //pair右侧复制右相机立体校正后的照片  
			cvGetCols(pair, &part, imSize.width, imSize.width * 2);
			cvCvtColor(img2r, &part, CV_GRAY2BGR);//GRAY->RGB  
												  //在一定距离间隔处绘制水平直线  
			for (int j = 0; j < imSize.height; j += 16)
				cvLine(pair, cvPoint(0, j), cvPoint(imSize.width * 2, j), CV_RGB(0, 255, 0));
		}
		else//竖直布置相机的情况  
		{
			//pair的上侧复制上相机立体校正后的照片  
			cvGetRows(pair, &part, 0, imSize.height);
			cvCvtColor(img1r, &part, CV_GRAY2BGR);//GRAY->RGB  
												  //pair的下侧复制下相机立体校正后的照片  
			cvGetRows(pair, &part, imSize.height, imSize.height * 2);
			cvCvtColor(img2r, &part, CV_GRAY2BGR);//GRAY->RGB  
												  //在一定距离间隔处绘制竖直的直线  
			for (int j = 0; j < imSize.width; j += 16)
				cvLine(pair, cvPoint(j, 0), cvPoint(j, imSize.height * 2), CV_RGB(0, 255, 0));
		}
		//显示pair  
		cvShowImage("rectified", pair);
		//等待用户输入。。。  
		if (cvWaitKey() == 27)
			break;
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n_boards = 13;//标定板图片数量  
	int board_w = 8;//x方向角点数目  
	int board_h = 6;//y方向角点数目  
	int board_n = board_w*board_h;//n_boards张图片中角点总数量  

								  //标定左相机  
	const char* left_imageList = "left_ImgName.txt";
	CvMat* left_object_points = cvCreateMat(n_boards*board_n, 3, CV_32FC1);//左相机对应的世界坐标中的角点  
	CvMat* left_image_points = cvCreateMat(n_boards*board_n, 2, CV_32FC1);//左相机图片坐标中的角点  
	CvMat* left_intrinsic_matrix = cvCreateMat(3, 3, CV_32FC1);//左相机内参数矩阵  
	CvMat* left_distortion_coeffs = cvCreateMat(5, 1, CV_32FC1);//左相机畸变系数向量  
	CvSize* left_imgSize = new CvSize;//图片尺寸  
									  //标定左相机内参数和畸变向量  
	SingleCalib(left_imageList, left_object_points, left_image_points, left_intrinsic_matrix, left_distortion_coeffs,
		n_boards, board_w, board_h, left_imgSize);
	//Show(用于调试)  
	/*int cout_step = left_intrinsic_matrix->step / sizeof(float);
	float* data = left_intrinsic_matrix->data.fl;
	for (int i = 0; i < 3; i++)
	{
	cout << "\n";
	for (int j = 0; j < 3; j++)
	{
	cout << (data + i*cout_step)[j] << "      ";
	}
	}*/

	//标定右相机  
	const char* right_imageList = "right_ImgName.txt";
	CvMat* right_object_points = cvCreateMat(n_boards*board_n, 3, CV_32FC1);//右相机对应的世界坐标中的角点  
	CvMat* right_image_points = cvCreateMat(n_boards*board_n, 2, CV_32FC1);//右相机图片坐标中的角点  
	CvMat* right_intrinsic_matrix = cvCreateMat(3, 3, CV_32FC1);//右相机内参数矩阵  
	CvMat* right_distortion_coeffs = cvCreateMat(5, 1, CV_32FC1);//右相机畸变系数向量  
	CvSize* right_imgSize = new CvSize;//图片尺寸  
									   //标定右相机内参数和畸变向量  
	SingleCalib(right_imageList, right_object_points, right_image_points, right_intrinsic_matrix, right_distortion_coeffs,
		n_boards, board_w, board_h, right_imgSize);
	//输出标定的内参数矩阵  
	/*cout_step = right_intrinsic_matrix->step / sizeof(float);
	data = right_intrinsic_matrix->data.fl;
	for (int i = 0; i < 3; i++)
	{
	cout << "\n";
	for (int j = 0; j < 3; j++)
	{
	cout << (data + i*cout_step)[j] << "      ";
	}
	}*/
	//图片尺寸（取左右相机的都一样）  
	CvSize imgSize;
	imgSize = *right_imgSize;
	CvMat* R = cvCreateMat(3, 3, CV_32FC1);//旋转矩阵  
	CvMat* T = cvCreateMat(3, 1, CV_32FC1);//平移矩阵  
	CvMat* E = cvCreateMat(3, 3, CV_32FC1);//本征矩阵  
	CvMat* F = cvCreateMat(3, 3, CV_32FC1);//基础矩阵  
										   /*立体标定，得到R, T, E, F，同时可以校正right_image_points, left_intrinsic_matrix, right_intrinsic_matrix,
										   left_distortion_coeffs*/
	StereoCalib(left_object_points, left_image_points, right_image_points, left_intrinsic_matrix, right_intrinsic_matrix,
		left_distortion_coeffs, right_distortion_coeffs, R, T, E, F, n_boards, board_w, board_h, imgSize);

	//输出立体标定结果  
	cout << endl << "M1:";
	for (int i = 0; i < 3; i++)
	{
		cout << endl;
		for (int j = 0; j < 3; j++)
		{
			cout << CV_MAT_ELEM(*left_intrinsic_matrix, float, i, j) << "     ";
		}
	}
	cout << endl << "D1:";
	for (int i = 0; i < 5; i++)
	{
		cout << endl;
		cout << CV_MAT_ELEM(*left_distortion_coeffs, float, i, 0) << "     ";
	}
	cout << endl << "M2:";
	for (int i = 0; i < 3; i++)
	{
		cout << endl;
		for (int j = 0; j < 3; j++)
		{
			cout << CV_MAT_ELEM(*right_intrinsic_matrix, float, i, j) << "     ";
		}
	}
	cout << endl << "D2:";
	for (int i = 0; i < 5; i++)
	{
		cout << endl;
		cout << CV_MAT_ELEM(*right_distortion_coeffs, float, i, 0) << "     ";
	}
	cout << endl << "R:";
	for (int i = 0; i < 3; i++)
	{
		cout << endl;
		for (int j = 0; j < 3; j++)
		{
			cout << CV_MAT_ELEM(*R, float, i, j) << "     ";
		}
	}
	cout << endl << "T:";
	for (int i = 0; i < 3; i++)
	{
		cout << endl;
		cout << CV_MAT_ELEM(*T, float, i, 0) << "     ";
	}
	cout << endl << "E:";
	for (int i = 0; i < 3; i++)
	{
		cout << endl;
		for (int j = 0; j < 3; j++)
		{
			cout << CV_MAT_ELEM(*E, float, i, j) << "     ";
		}
	}
	cout << endl << "F:";
	for (int i = 0; i < 3; i++)
	{
		cout << endl;
		for (int j = 0; j < 3; j++)
		{
			cout << CV_MAT_ELEM(*F, float, i, j) << "     ";
		}
	}

	//计算标定结果的误差  
	double Error = 0;
	Error = Calib_Quality_Check(left_image_points, right_image_points, left_intrinsic_matrix, right_intrinsic_matrix,
		left_distortion_coeffs, right_distortion_coeffs, F, n_boards, board_w, board_h);
	//输出标定误差  
	cout << endl << "Error:" << Error << endl;

	//保存标定结果  
	CvFileStorage* fs = cvOpenFileStorage("Result.xml", 0, CV_STORAGE_WRITE);
	cvWrite(fs, "left_intrinsic_matrix", left_intrinsic_matrix, cvAttrList(0, 0));
	cvWrite(fs, "left_distortion", right_distortion_coeffs, cvAttrList(0, 0));
	cvWrite(fs, "right_intrinsic_matrix", right_intrinsic_matrix, cvAttrList(0, 0));
	cvWrite(fs, "right_distortion", right_distortion_coeffs, cvAttrList(0, 0));
	cvWrite(fs, "R", R, cvAttrList(0, 0));
	cvWrite(fs, "T", T, cvAttrList(0, 0));
	cvWrite(fs, "E", E, cvAttrList(0, 0));
	cvWrite(fs, "F", F, cvAttrList(0, 0));
	cvWriteReal(fs, "Error", Error);
	cvReleaseFileStorage(&fs);

	//计算立体校正的映射矩阵mapx和mapy  
	CvMat* mxl = cvCreateMat(imgSize.height, imgSize.width, CV_32F);
	CvMat* myl = cvCreateMat(imgSize.height, imgSize.width, CV_32F);
	CvMat* mxr = cvCreateMat(imgSize.height, imgSize.width, CV_32F);
	CvMat* myr = cvCreateMat(imgSize.height, imgSize.width, CV_32F);
	ComputeRectification(left_image_points, right_image_points, left_intrinsic_matrix, right_intrinsic_matrix, left_distortion_coeffs, right_distortion_coeffs,
		R, T, F, imgSize, n_boards, board_w, board_h, mxl, myl, mxr, myr);

	//校正图片和计算视察映射并显示  
	Rectify_and_ComputeDisparity(left_imageList, right_imageList, imgSize, mxl, myl, mxr, myr);


	return 0;
}
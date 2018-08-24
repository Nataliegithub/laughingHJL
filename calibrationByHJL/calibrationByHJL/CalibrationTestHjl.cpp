/*��������궨������
���������
const char* imageList      IN����ͼƬ����txt�ļ�
CvMat* object_points       OUT��������ϵ��ľ���
CvMat* image_points        OUTͼ������ϵ����
CvMat* intrinsic_matrix    OUT���ص��ڲ�������
CvMat* distortion_coeffs   OUT���صĻ�������
int n_boards               INͼƬ������
int board_w                INÿ��ͼƬx����ǵ�����
int board_h                INÿ��ͼƬy����ǵ�����
CvSize* imgSize            OUTÿ��ͼƬ�����سߴ�
*/
static void SingleCalib(const char* imageList, CvMat* object_points, CvMat* image_points, CvMat* intrinsic_matrix, CvMat* distortion_coeffs,
	int n_boards, int board_w, int board_h, CvSize* imgSize)
{
	//�����ļ���  
	FILE* f;
	fopen_s(&f, imageList, "rt");

	int board_n = board_w*board_h;//ÿ��ͼƬ�нǵ�������  
	CvSize board_sz = cvSize(board_w, board_h);//�ǵ�ߴ����  

	CvPoint2D32f* corners = new CvPoint2D32f[board_w*board_h];//�������ڴ��ÿ��ͼƬ�ǵ�������һά������  
	CvMat* point_counts = cvCreateMat(n_boards, 1, CV_32SC1);//������ÿ��Ԫ�ش���ÿ��ͼƬ�ǵ������  
	int successes = 0;//�ҵ�ȫ���ǵ��ͼƬ����  
	int step = 0;//���ڼ�¼ÿ��ͼƬ�ǵ����ʼλ��  

				 //�ļ���ȡ���ɹ�:  
	if (!f)
	{
		fprintf(stderr, "can not open file %s\n", imageList);//Ҫ��д, ��֪���������, ������д��?stderr -- ��׼��������豸  
		return;
	}
	//����iѭ����ȡ�ļ��е��ַ���Ȼ�����ڶ�ȡͼƬ  
	for (int i = 0;; i++)
	{
		//��ȡͼƬ  
		char buf[1024];//��Ŷ�ȡ���ַ�����  
		int count = 0, result = 0;//count�ҵĵĽǵ�������result�ҽǵ�����־��ȫ���ǵ��ҵ����㣬����Ϊ0��  

		if (!fgets(buf, sizeof(buf) - 3, f))//��ȡ�ļ����ַ���ŵ�buf  
			break;
		size_t len = strlen(buf);//lenΪ�ַ�����ĳ���  
		while (len > 0 && isspace(buf[len - 1]))//int isspace(int c)������c�Ƿ�Ϊ�ո��ַ���Ҳ�����ж��Ƿ�Ϊ�ո�(' ')����λ�ַ�(' \t ')��CR(' \r ')������(' \n ')  
			buf[--len] = '\0';                  //����ֱ��λ�ַ�(' \v ')��ҳ(' \f ')������������ڷǿհ��ַ��ĺ�����Ϊ��ӡ�\0��  
		if (buf[0] == '#')//��ͷΪ'#',�����˴�ѭ��  
			continue;
		IplImage* img = cvLoadImage(buf, 0);//��ȡͼƬ  
		if (!img)
			break;
		//��ȡͼƬ�ߴ�  
		*imgSize = cvGetSize(img);
		//��ȡ�ǵ�  
		result = cvFindChessboardCorners(img, cvSize(board_w, board_h), corners, &count, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_NORMALIZE_IMAGE);

		if (result)
		{
			//Calibration will suffer without subpixel interpolation  
			//���� cvFindCornerSubPix ͨ�����������־��������ؾ��ȵĽǵ�λ��  
			cvFindCornerSubPix(img, corners, count, cvSize(11, 11), cvSize(-1, -1), cvTermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));
			/*TermCriteria�����㷨����ֹ׼��:
			typedef struct CvTermCriteria0d
			{
			int    type;  CV_TERMCRIT_ITER ��CV_TERMCRIT_EPS��ֵ֮һ�����߶��ߵ����
			int    max_iter; ����������
			double epsilon; ����ľ�ȷ��
			}
			һ���ʾ������ֹ�����������ΪCV_TERMCRIT_ITER������������������Ϊ��ֹ���������ΪCV_TERMCRIT_EPS ���þ�����Ϊ����������
			���ΪCV_TERMCRIT_ITER+CV_TERMCRIT_EPS�����������������߾�����Ϊ�������������ĸ����������㡣*/

			//��ʼ����ǵ����ʵλ��  
			step = successes*board_n;
			//���ǵ������cornersѹ�����image_points���Լ�����Ӧ��object_points��ֵ  
			for (int i = step, j = 0; j < board_n; ++i, ++j)
			{
				CV_MAT_ELEM(*image_points, float, i, 0) = corners[j].x;
				CV_MAT_ELEM(*image_points, float, i, 1) = corners[j].y;
				CV_MAT_ELEM(*object_points, float, i, 0) = (j / board_w);
				CV_MAT_ELEM(*object_points, float, i, 1) = (j % board_w);
				CV_MAT_ELEM(*object_points, float, i, 2) = 0.0f;
			}
			//����ӦͼƬ��point_counts��ֵ  
			CV_MAT_ELEM(*point_counts, int, successes, 0) = board_n;
			successes++;
		}
		//�ͷŸýǵ�ͼ��  
		cvReleaseImage(&img);
	}
	//�ر��ļ�  
	fclose(f);

	//��ʼ������ڲξ���  
	CV_MAT_ELEM(*intrinsic_matrix, float, 0, 0) = 1.0f;
	CV_MAT_ELEM(*intrinsic_matrix, float, 1, 1) = 1.0f;
	//�궨������ڲξ���ͻ���ϵ������  
	cvCalibrateCamera2(object_points, image_points, point_counts, *imgSize, intrinsic_matrix, distortion_coeffs, NULL, NULL, 0);

}

static void StereoCalib(CvMat* _left_object_points, CvMat* _left_image_points, CvMat* _right_image_points, CvMat* left_intrinsic, CvMat* right_intrinsic,
	CvMat* left_distortion, CvMat* right_distortion, CvMat* _R, CvMat* _T, CvMat* _E, CvMat* _F, int _n_boards, int _board_w, int _board_h, CvSize _imgSize)
{
	int board_n = _board_w*_board_h;//ÿ��ͼƬ�нǵ�������  

									//��ʼ��Number_perImg  
	CvMat* Number_perImg = cvCreateMat(1, _n_boards, CV_32SC1);
	int* pInt;
	pInt = (int*)(Number_perImg->data.ptr);
	for (int i = 0; i < _n_boards; ++i)
	{
		*pInt = board_n;
		pInt++;
	}
	//Show(���ڵ���)  
	/*for (int i = 0; i < _n_boards; i++)
	{
	cout << CV_MAT_ELEM(*Number_perImg, int, 0, i) << endl;
	}*/


	//����궨.���� _R, _T, _E, _F,����ͬʱ����Number_perImg, left_intrinsic, D_left,right_intrinsic, D_right  
	cvStereoCalibrate(_left_object_points, _left_image_points, _right_image_points, Number_perImg, left_intrinsic, left_distortion,
		right_intrinsic, right_distortion, _imgSize, _R, _T, _E, _F,
		cvTermCriteria(CV_TERMCRIT_ITER + CV_TERMCRIT_EPS, 100, 1e-5), CV_CALIB_FIX_ASPECT_RATIO + CV_CALIB_ZERO_TANGENT_DIST + CV_CALIB_SAME_FOCAL_LENGTH);
}

/*����궨����
CvMat* left_image_points     IN�����ͼ������ϵ��ľ���
CvMat* right_image_points    IN�����ͼ������ϵ��ľ���
CvMat* left_intrinsic        INandOUT��������ڲξ��󣬾�����궨���������
CvMat* right_intrinsic       INandOUT��������ڲξ��󣬾�����궨���������
CvMat* left_distortion       INandOUT������Ļ���������������궨���������
CvMat* right_distortion      INandOUT������Ļ���������������궨���������
CvMat* _F                    OUT��������
int n_boards                 INͼƬ������
int board_w                  INÿ��ͼƬx����ǵ�����
int board_h                  INÿ��ͼƬy����ǵ�����
*/
double Calib_Quality_Check(CvMat* _left_image_points, CvMat* _right_image_points, CvMat* left_intrinsic, CvMat* right_intrinsic,
	CvMat* left_distortion, CvMat* right_distortion, CvMat* _F, int _n_boards, int _board_w, int _board_h)
{
	int board_n = _board_w*_board_h;//ÿ��ͼƬ�нǵ�������  
									//����LeftPoints����Ϊ_left_image_points����������_n_boards*board_n*3*CV_32FC1��left����������1*_n_boards*board_n*CV_32FC2���Բ�������ת��  
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
	//����LeftPoints����Ϊ_left_image_points����������_n_boards*board_n*3*CV_32FC1��left����������1*_n_boards*board_n*CV_32FC2���Բ�������ת��  
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
	//����������������Ϊleft_distortion����������5*1��������D_left����������1*5�������Բ�������ת��  
	CvMat* D_left = cvCreateMat(1, 5, CV_32FC1);
	for (int i = 0; i < 5; ++i)
	{
		CV_MAT_ELEM(*D_left, float, 0, i) = CV_MAT_ELEM(*left_distortion, float, i, 0);
	}
	//����������������Ϊright_distortion����������5*1��������D_right����������1*5�������Բ�������ת��  
	CvMat* D_right = cvCreateMat(1, 5, CV_32FC1);
	for (int i = 0; i < 5; ++i)
	{
		CV_MAT_ELEM(*D_right, float, 0, i) = CV_MAT_ELEM(*right_distortion, float, i, 0);
	}
	//���弫��  
	CvMat* Lines1 = cvCreateMat(1, _n_boards*board_n, CV_32FC3);
	CvMat* Lines2 = cvCreateMat(1, _n_boards*board_n, CV_32FC3);
	//У����Ļ���  
	cvUndistortPoints(_left_image_points, _left_image_points, left_intrinsic, D_left, 0, left_intrinsic);
	cvUndistortPoints(_right_image_points, _right_image_points, right_intrinsic, D_right, 0, right_intrinsic);
	//���㼫��  
	cvComputeCorrespondEpilines(_left_image_points, 1, _F, Lines1);
	cvComputeCorrespondEpilines(_right_image_points, 2, _F, Lines2);
	//���������ͼ���ͼ��ߵľ���  
	float a, b, c;//����ϵ��ax+by+c=0  
	double err, avgErr = 0;//������  
	float* p_temp = (float*)cvPtr2D(_left_image_points, 0, 0);//������ʱ����ĵ�  
	float* l_temp = (float*)cvPtr2D(Lines2, 0, 0, 0);//������ʱ����ļ���  
	for (int i = 0; i < _n_boards*board_n; i++)
	{
		//��ȡ������  
		x = *p_temp;
		y = *(p_temp + 1);
		p_temp += 2;
		//��ȡ����ϵ��  
		a = *l_temp;
		b = *(l_temp + 1);
		c = *(l_temp + 2);
		l_temp += 3;
		//����㵽ֱ�ߵľ���  
		err = fabs(a*x + b*y + c) / sqrt(a*a + b*b);
		//�ۼ����  
		avgErr += err;
	}
	//���������ͼ���ͼ��ߵľ���  
	p_temp = (float*)cvPtr2D(_right_image_points, 0, 0);//������ʱ����ĵ�  
	l_temp = (float*)cvPtr2D(Lines1, 0, 0);//������ʱ����ļ���  
	for (int i = 0; i < _n_boards*board_n; i++)
	{
		//��ȡ������  
		x = *p_temp;
		y = *(p_temp + 1);
		p_temp += 2;
		//��ȡ����ϵ��  
		a = *l_temp;
		b = *(l_temp + 1);
		c = *(l_temp + 2);
		l_temp += 3;
		//����㵽ֱ�ߵľ���  
		err = fabs(a*x + b*y + c) / sqrt(a*a + b*b);
		//�ۼ����  
		avgErr += err;
	}
	//������ƽ��ֵ  
	avgErr /= (_n_boards*board_n);
	return avgErr;
}


/*����У����ӳ�����mapx��mapy����
CvMat* left_image_points     IN�����ͼ������ϵ��ľ���
CvMat* right_image_points    IN�����ͼ������ϵ��ľ���
CvMat* left_intrinsic        INandOUT��������ڲξ��󣬾�����궨���������
CvMat* right_intrinsic       INandOUT��������ڲξ��󣬾�����궨���������
CvMat* left_distortion       INandOUT������Ļ���������������궨���������
CvMat* right_distortion      INandOUT������Ļ���������������궨���������
CvMat* R                     OUT��������ת����
CvMat* T                     OUT������ƽ�ƾ���
CvMat* F                     OUT��������
CvSize* imgSize              INÿ��ͼƬ�����سߴ�
int n_boards                 INͼƬ������
int board_w                  INÿ��ͼƬx����ǵ�����
int board_h                  INÿ��ͼƬy����ǵ�����
bool useUncalibrated         IN�Ƿ���÷Ǳ궨�㷨Hartly��������У����ӳ�����mapx��mapy
CvMat* mxl                   OUT���������У����ӳ�����mapx
CvMat* myl                   OUT���������У����ӳ�����mapy
CvMat* mxr                   OUT���������У����ӳ�����mapx
CvMat* myr                   OUT���������У����ӳ�����mapy
*/
/*mapx��mapy�Ǻ������صĲ���ӳ�������Ŀ��ͼ���ϵ�ÿһ��������˵��ӳ���˵��Ӧ�ô�ʲôλ�ò�ֵԴ���أ�
����ӳ������ֱ�ӱ�Remap()Ƕ��ʹ�á�����cvInitUndistortRectifyMap����������������ֱ���ã��������Ի�����ǵĸ��Բ�ͬ����ӳ�����mapx��mapy��
ÿ���������µ���������ͼ����ҪУ��ʱ����ʹ������ӳ���*/
static void ComputeRectification(CvMat* imagePoints_left, CvMat* imagePoints_right, CvMat* _left_intrinsic, CvMat* _right_intrinsic, CvMat* _left_distortion, CvMat* _right_distortion,
	CvMat* R, CvMat* T, CvMat* F, CvSize imageSize, int _n_boards, int _board_w, int _board_h, CvMat* mxl, CvMat* myl, CvMat* mxr, CvMat* myr)
{
	/*�������������������CV_64F���͵ģ���Ҳ��֪��Ϊʲô��cvStereoRectify�������CV_32FC1�Ĳ������жϣ���ΪCV_64FC1��ɹ�*/
	//����������ڲξ�����Ϊ_left_intrinsic��������ΪCV_32FC1��M_leftΪCV_64FC1  
	CvMat* M_left = cvCreateMat(3, 3, CV_64FC1);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			CV_MAT_ELEM(*M_left, double, i, j) = CV_MAT_ELEM(*_left_intrinsic, float, i, j);
	//����������ڲξ�����Ϊ_left_intrinsic��������ΪCV_32FC1��M_leftΪCV_64FC1  
	CvMat* M_right = cvCreateMat(3, 3, CV_64FC1);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			CV_MAT_ELEM(*M_right, double, i, j) = CV_MAT_ELEM(*_right_intrinsic, float, i, j);
	//����������������Ϊleft_distortion����������5*1��������D_left����������1*5������ͬʱת��ΪCV_64FC1��  
	CvMat* D_left = cvCreateMat(1, 5, CV_64FC1);
	for (int i = 0; i < 5; ++i)
	{
		CV_MAT_ELEM(*D_left, double, 0, i) = CV_MAT_ELEM(*_left_distortion, float, i, 0);
	}
	//����������������Ϊright_distortion����������5*1��������D_right����������1*5������ͬʱת��ΪCV_64FC1��  
	CvMat* D_right = cvCreateMat(1, 5, CV_64FC1);
	for (int i = 0; i < 5; ++i)
	{
		CV_MAT_ELEM(*D_right, double, 0, i) = CV_MAT_ELEM(*_right_distortion, float, i, 0);
	}
	//ת����ת����CV_32FC1ת��ΪCV_64FC1  
	CvMat* RR = cvCreateMat(3, 3, CV_64FC1);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			CV_MAT_ELEM(*RR, double, i, j) = CV_MAT_ELEM(*R, float, i, j);
	//����ƽ��������CV_32FC1ת��ΪCV_64FC1  
	CvMat* TT = cvCreateMat(3, 1, CV_64FC1);
	for (int i = 0; i < 3; ++i)
	{
		CV_MAT_ELEM(*TT, double, i, 0) = CV_MAT_ELEM(*T, float, i, 0);
	}
	//����ͳ�ʼ������������ж�׼��ת����  
	CvMat* Rl = cvCreateMat(3, 3, CV_64FC1);
	cvZero(Rl);
	CvMat* Rr = cvCreateMat(3, 3, CV_64FC1);
	cvZero(Rr);
	//�������������ж�׼��У����  
	if (useUncalibrated == false)//useUncalibrated == 0��ʹ��Bouguet��������У����  
	{
		//����ͳ�ʼ���������ͶӰ���󣬸þ�����Լ���У���������ڲ�������  
		CvMat* Pl = cvCreateMat(3, 4, CV_64FC1);
		cvZero(Pl);
		CvMat* Pr = cvCreateMat(3, 4, CV_64FC1);
		cvZero(Pr);
		//�������У�������õ� Rl, Rr, Pl, Pr  
		cvStereoRectify(M_left, M_right, D_left, D_right, imageSize,
			RR, TT, Rl, Rr, Pl, Pr, 0, 0, -1.0, cvSize(0, 0), (CvRect*)0, (CvRect*)0);
		//Show�����ڵ��ԣ�  
		/*cout << endl << "Pr:" << endl;
		for (int i = 0; i < 3; i++)
		{
		const double *p = (const double*)(Pr->data.ptr + i*Pr->step);
		cout << *p << "    " << *(p + 1) << "       " << *(p + 2)<<"        "<<*(p+3) << endl;
		}*/
		/*ͨ���Ƚ�ת������Pλ��[0][3]��λ��[1][3]�Ĵ�С �ж������ˮƽ���û��Ǵ�ֱ���á�
		��Ϊ������λ�õ�Ԫ�طֱ��ʾ���������ϵԭ���ƶ������������ϵԭ���x��y�����ƶ����롣*/
		isVerticalStereo = fabs(CV_MAT_ELEM(*Pr, double, 1, 3)) > fabs(CV_MAT_ELEM(*Pr, double, 0, 3));
		//Precompute maps for cvRemap()  
		//����ͨ������cvInitUndistortRectifyMap()��������У����ӳ�����mx��my  
		cvInitUndistortRectifyMap(_left_intrinsic, D_left, Rl, Pl, mxl, myl);
		cvInitUndistortRectifyMap(_right_intrinsic, D_right, Rr, Pr, mxr, myr);
	}
	else//useUncalibrated�� = 0��ʹ��Hartley��������У����  
	{
		int board_n = _board_w*_board_h;//ÿ��ͼƬ�нǵ�������  
		CvMat* Hl = cvCreateMat(3, 3, CV_64FC1);//Hartley��������õ���У����  
		CvMat* Hr = cvCreateMat(3, 3, CV_64FC1);//Hartley��������õ���У����  
		CvMat *invM = cvCreateMat(3, 3, CV_64FC1);//������ʱ�������ڲ�������  
												  //����imagePoints_left��imagePoints_right  
		float x, y;//���ڴ�ֵ����ʱ����  
		double* p;//���ڴ�ֵ����ʱ����  
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
		//ͨ�����������ͼ�������������F  
		cvFindFundamentalMat(left, right, F);
		//ʹ�÷Ǳ궨�㷨Hartley����У����������������У����Ӧ�Ծ���H1,H2  
		cvStereoRectifyUncalibrated(left, right, F, imageSize, Hl, Hr, 3);
		/*ʹ��cvStereoRectifyUncalibrated()�Ľ��H1,H2��У�������������������Ԥ����һ�µ�Ӧ�Ծ���Hl����Hr��
		���û����ǰ�궨�л��Mrect����Mrect=M,�ֱ��������У����ת����Rl=inv(Mrect_l)Hl(Ml)��Rr=inv(Mrect_r)Hr(Mr)��
		�����Ҫ�ֱ������������Ļ���ϵ������*/
		//Show�����ڵ��ԣ�  
		/*cout << endl << "M1:";
		for (int i = 0; i < 3; i++)
		{
		cout << endl;
		for (int j = 0; j < 3; j++)
		{
		cout << CV_MAT_ELEM(*_left_intrinsic, float, i, j) << "     ";
		}
		}*/
		//Ԥ����Rl=inv(Mrect_l)Hl(Ml)  
		cvInvert(M_left, invM);
		cvMatMul(Hl, M_left, Rl);
		cvMatMul(invM, Rl, Rl);
		//Ԥ����Rr=inv(Mrect_r)Hr(Mr)  
		cvInvert(M_right, invM);//_iMΪ_M1�������  
		cvMatMul(Hr, M_right, Rr);//_R1=_H1*_M1  
		cvMatMul(invM, Rr, Rr);//_R2=_iM*_R1  
							   //Precompute map for cvRemap()����Ԥ������_R1��_R2���Լ�ǰ�����������������cvInitUndistortRectifyMap����mapx��mapy  
		cvInitUndistortRectifyMap(M_left, D_left, Rl, M_left, mxl, myl);
		cvInitUndistortRectifyMap(M_right, D_right, Rr, M_right, mxr, myr);

	}
}

/*У��ͼƬ�ͼ����Ӳ�ӳ��ĺ���
const char* imageList_left       IN���������ͼƬ����txt�ļ�
const char* imageList_right      IN���������ͼƬ����txt�ļ�
CvSize* imgSize                  INÿ��ͼƬ�����سߴ�
CvMat* mxl                       IN���������У����ӳ�����mapx
CvMat* myl                       IN���������У����ӳ�����mapy
CvMat* mxr                       IN���������У����ӳ�����mapx
CvMat* myr                       IN���������У����ӳ�����mapy
*/
static void Rectify_and_ComputeDisparity(const char* imageList_left, const char* imageList_right, CvSize imSize, CvMat* mxl, CvMat* myl, CvMat* mxr, CvMat* myr)
{
	FILE *f_left, *f_right;//���������ļ���������  
	fopen_s(&f_left, imageList_left, "rt");//�������ͼ���ͼ�����ļ�  
	fopen_s(&f_right, imageList_right, "rt");//�������ͼ���ͼ�����ļ�  
	IplImage *img_left, *img_right;//��������ͼ�����  
								   //���������������ڴ��У�����ͼ��  
	CvMat* img1r = cvCreateMat(imSize.height, imSize.width, CV_8U);
	CvMat* img2r = cvCreateMat(imSize.height, imSize.width, CV_8U);

	CvMat* disp = cvCreateMat(imSize.height, imSize.width, CV_16S);//�Ӳ�ӳ��  
	CvMat* vdisp = cvCreateMat(imSize.height, imSize.width, CV_8U);//��һ������Ӳ�ӳ��  
																   /*�����ͳ�ʼ��һ���ṹ��CvStereoBMState����cvFindStereoCorrespondenceBM()��ƥ����̵Ĳ�������Ϊ���ÿ�ε���cvFindStereoCorrespondenceBM()
																   ����ɢ�����ڲ���ɢ�����Ӱ�����������ٶȡ�*/
	CvStereoBMState *BMState = cvCreateStereoBMState();
	assert(BMState != 0);//�������ʧ������ֹ����  
	BMState->preFilterSize = 41;//Ԥ���˴��ڵĳߴ磬ʹͼ�����ȹ�һ������ǿ����  
	BMState->preFilterCap = 31;//Ԥ���˵���ֵ  
	BMState->SADWindowSize = 41;//"��������ۼ�"С����(SAD)�ĳߴ�  
	BMState->minDisparity = -64;//ƥ��������ʼ��λ��  
	BMState->numberOfDisparities = 128;//ƥ������ ���ӲΧ  
	BMState->textureThreshold = 10;//�������ֵ����ӦС�ڸ�ֵ��ƥ�䲻�迼��  
	BMState->uniquenessRatio = 15;//����˱�����ֵ��uniquenessRatio>(match_val-min_match)/min_match�ĵط����˵�ƥ��ֵ  

	cvNamedWindow("disparity");//����һ������������ʾ�Ӳ�ӳ��  
	cvNamedWindow("rectified", 1);//�������ڣ�������ʾУ�����ͼ��  
	CvMat* pair;//����һ���������ڴ��У���������ͼƬ  
	CvMat part;//����һ���������ڴ��У�����һ��ͼƬ  
	if (!isVerticalStereo)//���ˮƽ�����򴴽��ߴ�Ϊ��imageSize.height, imageSize.width * 2����ͼƬ������ʾ�����������У�����ͼƬ  
		pair = cvCreateMat(imSize.height, imSize.width * 2, CV_8UC3);
	else//�����ֱ�����򴴽��ߴ�Ϊ��imageSize.height * 2, imageSize.width����ͼƬ������ʾ�����������У�����ͼƬ  
		pair = cvCreateMat(imSize.height * 2, imSize.width, CV_8UC3);

	//ѭ����ȡ���������ͼ��  
	for (int i = 0;; i++)
	{
		//��ȡͼƬ  
		char buf_left[1024], buf_right[1024];//��Ŷ�ȡ���ַ�����  
		if ((!fgets(buf_left, sizeof(buf_left) - 3, f_left)) || (!fgets(buf_right, sizeof(buf_right) - 3, f_right)))//��ȡ�ļ����ַ���ŵ�buf  
			break;
		size_t len = strlen(buf_left);//lenΪ�ַ�����ĳ���  
		while (len > 0 && isspace(buf_left[len - 1]))//int isspace(int c)������c�Ƿ�Ϊ�ո��ַ���Ҳ�����ж��Ƿ�Ϊ�ո�(' ')����λ�ַ�(' \t ')��CR(' \r ')������(' \n ')  
			buf_left[--len] = '\0';                 //����ֱ��λ�ַ�(' \v ')��ҳ(' \f ')������������ڷǿհ��ַ��ĺ�����Ϊ��ӡ�\0��  
		len = strlen(buf_right);//lenΪ�ַ�����ĳ���  
		while (len > 0 && isspace(buf_right[len - 1]))//int isspace(int c)������c�Ƿ�Ϊ�ո��ַ���Ҳ�����ж��Ƿ�Ϊ�ո�(' ')����λ�ַ�(' \t ')��CR(' \r ')������(' \n ')  
			buf_right[--len] = '\0';                    //����ֱ��λ�ַ�(' \v ')��ҳ(' \f ')������������ڷǿհ��ַ��ĺ�����Ϊ��ӡ�\0��  
		if ((buf_left[0] == '#') || (buf_right[0] == '#'))//��ͷΪ'#',�����˴�ѭ��  
			continue;
		img_left = cvLoadImage(buf_left, 0);//��ȡͼƬ  
		img_right = cvLoadImage(buf_right, 0);//��ȡͼƬ  
		if ((!img_left) || (!img_right))
			break;
		//ʹ����ǰ������mx��my����У��ͼƬ�����������ͼƬ����ͬһ�����ƽ�����ж�׼  
		cvRemap(img_left, img1r, mxl, myl);
		cvRemap(img_right, img2r, mxr, myr);
		/*������������У׼���ͼƬ����˫Ŀ�Ӿ�ģ�͵��Ӳ�ӳ�䡣��Ϊ���ˮƽ���úʹ�ֱ�������������
		�����������ֱ���õ����Σ����û���Լ����ͼ��ת�õĴ��룬����cvFindStereoCorrespondenceBM()
		ֻ�ܼ���Ǳ궨У�����Ӳ죬���������ˮƽ��������Σ�����cvFindStereoCorrespondenceBM()���Լ���
		�궨�ͷǱ궨��У��������Ե��Ӳ*/
		if (!isVerticalStereo || useUncalibrated != 0)
		{
			cvFindStereoCorrespondenceBM(img1r, img2r, disp, BMState);//�����Ӳ�ӳ��  
			cvNormalize(disp, vdisp, 0, 256, CV_MINMAX);//��һ���Ӳ�ӳ��  
														//��ʾ�Ӳ�ӳ��  
			cvShowImage("disparity", vdisp);
		}
		//������ʾ����У��������������Ӧ��ͼƬ  
		if (!isVerticalStereo)//ˮƽ�������������  
		{
			//pair����ิ�����������У�������Ƭ  
			cvGetCols(pair, &part, 0, imSize.width);/*CvMat* cvGetCols( const CvArr* arr, CvMat* submat, int start_col, int end_col );
													���������һ������ڵ���
													arr-��������
													submat-ָ����������ͷָ��.
													start_col-��ȵĿ�ʼ�У��������У������±꣬���±���0Ϊ��׼��
													end_col-��ȵĽ����У����������У������±꣬���±���0Ϊ��׼��*/
			cvCvtColor(img1r, &part, CV_GRAY2BGR);//GRAY->RGB  
												  //pair�Ҳิ�����������У�������Ƭ  
			cvGetCols(pair, &part, imSize.width, imSize.width * 2);
			cvCvtColor(img2r, &part, CV_GRAY2BGR);//GRAY->RGB  
												  //��һ��������������ˮƽֱ��  
			for (int j = 0; j < imSize.height; j += 16)
				cvLine(pair, cvPoint(0, j), cvPoint(imSize.width * 2, j), CV_RGB(0, 255, 0));
		}
		else//��ֱ������������  
		{
			//pair���ϲิ�����������У�������Ƭ  
			cvGetRows(pair, &part, 0, imSize.height);
			cvCvtColor(img1r, &part, CV_GRAY2BGR);//GRAY->RGB  
												  //pair���²ิ�����������У�������Ƭ  
			cvGetRows(pair, &part, imSize.height, imSize.height * 2);
			cvCvtColor(img2r, &part, CV_GRAY2BGR);//GRAY->RGB  
												  //��һ����������������ֱ��ֱ��  
			for (int j = 0; j < imSize.width; j += 16)
				cvLine(pair, cvPoint(j, 0), cvPoint(j, imSize.height * 2), CV_RGB(0, 255, 0));
		}
		//��ʾpair  
		cvShowImage("rectified", pair);
		//�ȴ��û����롣����  
		if (cvWaitKey() == 27)
			break;
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n_boards = 13;//�궨��ͼƬ����  
	int board_w = 8;//x����ǵ���Ŀ  
	int board_h = 6;//y����ǵ���Ŀ  
	int board_n = board_w*board_h;//n_boards��ͼƬ�нǵ�������  

								  //�궨�����  
	const char* left_imageList = "left_ImgName.txt";
	CvMat* left_object_points = cvCreateMat(n_boards*board_n, 3, CV_32FC1);//�������Ӧ�����������еĽǵ�  
	CvMat* left_image_points = cvCreateMat(n_boards*board_n, 2, CV_32FC1);//�����ͼƬ�����еĽǵ�  
	CvMat* left_intrinsic_matrix = cvCreateMat(3, 3, CV_32FC1);//������ڲ�������  
	CvMat* left_distortion_coeffs = cvCreateMat(5, 1, CV_32FC1);//���������ϵ������  
	CvSize* left_imgSize = new CvSize;//ͼƬ�ߴ�  
									  //�궨������ڲ����ͻ�������  
	SingleCalib(left_imageList, left_object_points, left_image_points, left_intrinsic_matrix, left_distortion_coeffs,
		n_boards, board_w, board_h, left_imgSize);
	//Show(���ڵ���)  
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

	//�궨�����  
	const char* right_imageList = "right_ImgName.txt";
	CvMat* right_object_points = cvCreateMat(n_boards*board_n, 3, CV_32FC1);//�������Ӧ�����������еĽǵ�  
	CvMat* right_image_points = cvCreateMat(n_boards*board_n, 2, CV_32FC1);//�����ͼƬ�����еĽǵ�  
	CvMat* right_intrinsic_matrix = cvCreateMat(3, 3, CV_32FC1);//������ڲ�������  
	CvMat* right_distortion_coeffs = cvCreateMat(5, 1, CV_32FC1);//���������ϵ������  
	CvSize* right_imgSize = new CvSize;//ͼƬ�ߴ�  
									   //�궨������ڲ����ͻ�������  
	SingleCalib(right_imageList, right_object_points, right_image_points, right_intrinsic_matrix, right_distortion_coeffs,
		n_boards, board_w, board_h, right_imgSize);
	//����궨���ڲ�������  
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
	//ͼƬ�ߴ磨ȡ��������Ķ�һ����  
	CvSize imgSize;
	imgSize = *right_imgSize;
	CvMat* R = cvCreateMat(3, 3, CV_32FC1);//��ת����  
	CvMat* T = cvCreateMat(3, 1, CV_32FC1);//ƽ�ƾ���  
	CvMat* E = cvCreateMat(3, 3, CV_32FC1);//��������  
	CvMat* F = cvCreateMat(3, 3, CV_32FC1);//��������  
										   /*����궨���õ�R, T, E, F��ͬʱ����У��right_image_points, left_intrinsic_matrix, right_intrinsic_matrix,
										   left_distortion_coeffs*/
	StereoCalib(left_object_points, left_image_points, right_image_points, left_intrinsic_matrix, right_intrinsic_matrix,
		left_distortion_coeffs, right_distortion_coeffs, R, T, E, F, n_boards, board_w, board_h, imgSize);

	//�������궨���  
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

	//����궨��������  
	double Error = 0;
	Error = Calib_Quality_Check(left_image_points, right_image_points, left_intrinsic_matrix, right_intrinsic_matrix,
		left_distortion_coeffs, right_distortion_coeffs, F, n_boards, board_w, board_h);
	//����궨���  
	cout << endl << "Error:" << Error << endl;

	//����궨���  
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

	//��������У����ӳ�����mapx��mapy  
	CvMat* mxl = cvCreateMat(imgSize.height, imgSize.width, CV_32F);
	CvMat* myl = cvCreateMat(imgSize.height, imgSize.width, CV_32F);
	CvMat* mxr = cvCreateMat(imgSize.height, imgSize.width, CV_32F);
	CvMat* myr = cvCreateMat(imgSize.height, imgSize.width, CV_32F);
	ComputeRectification(left_image_points, right_image_points, left_intrinsic_matrix, right_intrinsic_matrix, left_distortion_coeffs, right_distortion_coeffs,
		R, T, F, imgSize, n_boards, board_w, board_h, mxl, myl, mxr, myr);

	//У��ͼƬ�ͼ����Ӳ�ӳ�䲢��ʾ  
	Rectify_and_ComputeDisparity(left_imageList, right_imageList, imgSize, mxl, myl, mxr, myr);


	return 0;
}
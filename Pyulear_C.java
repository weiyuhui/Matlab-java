package pyulear;

public class Pyulear_C {

	private double[] data;
	private static int rst_len;
	//private static double[] prstdata; 返回值
	private static int AR_LEAST_ORDER = 60;

	public Pyulear_C(double[] data, int rst_len) {
		this.data = data;
		//this.prstdata = prstdata;
		//this.AR_LEAST_ORDER = AR_LEAST_ORDER;
		this.rst_len = rst_len;
	}

	public static double[] freqMe(double[] data, int rst_len, boolean bRemoveDC) throws Exception {
		// TODO Auto-generated method stub
		double[] prstdata = new double[rst_len];
		if (data == null | data.length < AR_LEAST_ORDER) {
			return null;
		}
		double[] pdata_buf = null;
		if (bRemoveDC == true) {
			pdata_buf = new double[data.length];
			double dc = 0;
			for (int j = 0; j < data.length; j++) {
				dc += data[j];
			}
			dc = dc / data.length;
			for (int j = 0; j < data.length; j++) {
				pdata_buf[j] = data[j] - dc;
			}
		}
		
		//C
		
		int degree = 8;
		int length = data.length;
	    double[][] ar = new double[degree+1][degree+1];
	    double[] coef = new double[degree];
	    double[] per = new double[length+1];  
	    double[] pef = new double[length+1];   
	    double[] h = new double[degree+1]; 
	    double[] g = new double[degree+2];
//	    double t1, t2;
//	    int n;
//
//	    for (n = 1; n <= degree; n++)
//	    {
//	        double sn = 0.0;
//	        double sd = 0.0;
//	        int j;
//	        int jj = length - n;
//
//	        for (j = 0; j < jj; j++)
//	        {
//	            t1 = data [j + n] + pef [j];
//	            t2 = data [j] + per [j];
//	            sn -= 2.0 * t1 * t2;
//	            sd += (t1 * t1) + (t2 * t2);
//	        }
//
//	        t1 = g [n] = sn / sd;
//	        if (n != 1)
//	        {
//	            for (j = 1; j < n; j++) 
//	                h [j] = g [j] + t1 * g [n - j];
//	            for (j = 1; j < n; j++)
//	                g [j] = h [j];
//	            jj--;
//	        }
//
//	        for (j = 0; j < jj; j++)
//	        {
//	            per [j] += t1 * pef [j] + t1 * data [j + n];
//	            pef [j] = pef [j + 1] + t1 * per [j + 1] + t1 * data [j + 1];
//	        }
//
//	        for (j = 0; j < n; j++)
//	           ar [n][j] = g [j + 1];
//
//	    for (n = 0; n < degree; n++)
//	        coef [n] = g [n + 1];
//	    }
	    
	    
	    int j,n,nn,jj;
	    double sn,sd;
	    double t1,t2;

	    for (j=1;j<=length;j++) {
	       pef[j] = 0;
	       per[j] = 0;
	    }
	       
	    for (nn=2;nn<=degree+1;nn++) {
	       n  = nn - 2;
	       sn = 0.0;
	       sd = 0.0;
	       jj = length - n - 1;
	       for (j=1;j<=jj;j++) {
	          t1 = data[j+n] + pef[j];
	          t2 = data[j-1] + per[j];
	          sn -= 2.0 * t1 * t2;
	          sd += (t1 * t1) + (t2 * t2);
	       }
	       g[nn] = sn / sd;
	       t1 = g[nn];
	       if (n != 0) {
	          for (j=2;j<nn;j++) 
	             h[j] = g[j] + (t1 * g[n - j + 3]);
	          for (j=2;j<nn;j++)
	             g[j] = h[j];
	          jj--;
	       }
	       for (j=1;j<=jj;j++) {
	          per[j] += (t1 * pef[j]) + (t1 * data[j+nn-2]);
	          pef[j]  = pef[j+1] + (t1 * per[j+1]) + (t1 * data[j]);
	       }

	       for (j=2;j<=nn;j++)
	          ar[nn-1][j-1] = g[j];
	    }
	    
	    for (int i=1;i<=degree;i++)
	         coef[i-1] = -ar[degree][i];
	    
	    


		/*
		 * 0-0.5之间选择频率成分，所以实际可能需要运行三次，得到三个不懂的频段数据。
		 * 1.看看Matlab的频段划分是如何实现的
		 * 2.缩减选择阶数的那个环节；
		 * 3.整理返回数据的长度与AR阶数的关系
		 */
		double[] coe = AutoRegression.calculateARCoefficients(data, 8, true);
		// reserve
		double[] coer = new double[coe.length];
		for(int i =0;i<coe.length;i++){
			coer[i] = coe[coe.length-1-i];
		}
		/*
		 * 1.coef:老师给的c语言网站上转换而来的函数；
		 * 2.coe:AutoRegression函数得来的
		 * 3.coer:coe经过反序得到的
		 * 1、2、3的阶数都是8
		 */
		prstdata = CalData(coef, 0.01, 0.5);

		double sum1, sum2;

		if (bRemoveDC == true) {
			sum1 = 0.f;
			for (int j3 = 0; j3 < data.length; j3++) {
				sum1 += pdata_buf[j3] * pdata_buf[j3];
			}
		} else {
			sum1 = 0;
			for (int j3 = 0; j3 < data.length; j3++) {
				sum1 += data[j3] * data[j3];
			}
		}

		sum2 = 0;
		for (int j3 = 0; j3 < rst_len; j3++) {
			sum2 += prstdata[j3];
		}

		double coeff;
		if (sum2 != 0.f)
			coeff = sum1 / sum2;
		else
			coeff = 1.f;

		for (int j3 = 0; j3 < rst_len; j3++) {
			prstdata[j3] = prstdata[j3] * coeff;
		}

		// if (pdata_buf!=NULL)
		// delete []pdata_buf;

		return prstdata;
	}

	private static double[] CalData(double[] porgpara, double f0, double f1) {
		// TODO Auto-generated method stub
		if (f0<0.) f0=0.;
		if (f0>0.49) f0=0.49;
		if (f1<0.01) f1=0.01;
		if (f1>0.5) f1=0.5;
		if (f1<f0) f1=f0+0.01;
		
		int datalen = rst_len;
		int paralen = porgpara.length;//返回数值的长度，自己定义
		double[] pdata= new double[datalen];

		double f=f0;
		double det_f=((f1-f0)/datalen);
		int	nf=datalen;
		double	sumr, sumi, tmpf;
		int j, k;
		for (j=1; j<=nf; j++)
		{
			tmpf=(2*Math.PI*f);
			sumr=1.; sumi=0.;
			for (k=1; k<paralen; k++)
			{
				sumr=sumr+(porgpara[k]*Math.cos(tmpf*k));
				sumi=sumi-(porgpara[k]*Math.sin(tmpf*k));
			}
			tmpf=(sumr*sumr+sumi*sumi);
			if (tmpf>0)
				pdata[j-1]=(porgpara[0]/(sumr*sumr+sumi*sumi));
			else
				pdata[j-1] = 0;
			f+=det_f;
		}

		return pdata;
		
		

	}


}
